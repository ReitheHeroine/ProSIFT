#!/usr/bin/env python3
# title: uniprot_mapping.py
# project: ProSIFT
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-03-26
# last modified: 2026-03-26
#
# purpose:
#   Processes 4.5-4.9 of Module 01 (Input Validation). Takes the
#   detection-filtered protein list, queries the UniProt bulk ID mapping
#   API to retrieve gene symbols, Entrez IDs, and Ensembl gene IDs, handles
#   edge cases (unmapped, isoform-level accessions, multiple mappings), and
#   caches results to avoid redundant API calls across re-runs.
#
# inputs:
#   - {run_id}.filtered_matrix.parquet  (from filter_proteins.py)
#   - {run_id}_params.yml               (ProSIFT run parameters)
#
# outputs:
#   - {run_id}.id_mapping.parquet       (one row per input protein ID)
#   - {run_id}.validation_report_part3.txt
#
# usage example:
#   python bin/uniprot_mapping.py \
#       --matrix CTXcyto_WT_vs_CTXcyto_KO.filtered_matrix.parquet \
#       --params CTXcyto_WT_vs_CTXcyto_KO_params.yml \
#       --run_id CTXcyto_WT_vs_CTXcyto_KO \
#       --cachedir results/cache/id_mapping \
#       --outdir .

import argparse
import hashlib
import json
import os
import re
import sys
import time
from datetime import datetime, timezone
from typing import Optional

try:
    import pandas as pd
    import requests
    import yaml
except ImportError as exc:
    print(f'Error: required package not available -- {exc}')
    print('Activate the prosift conda environment before running.')
    sys.exit(1)


# ============================================================
# Constants
# ============================================================

UNIPROT_BASE   = 'https://rest.uniprot.org'
IDMAP_RUN_URL  = f'{UNIPROT_BASE}/idmapping/run'
IDMAP_STATUS   = f'{UNIPROT_BASE}/idmapping/status'

FROM_TYPE   = 'UniProtKB_AC-ID'
TO_TYPE     = 'UniProtKB'

# Page size for paginated results fetch
RESULTS_PAGE_SIZE = 500

RETRY_DELAYS   = [5, 15, 45]   # seconds between job-submission retries
POLL_INTERVAL  = 5             # seconds between status polls
MAX_POLLS      = 120           # give up after 10 minutes (120 * 5s)

# Output column schema (Process 4.8)
SCHEMA_COLS = [
    'input_id',
    'uniprot_accession',
    'gene_symbol_mouse',
    'entrez_id_mouse',
    'ensembl_gene_mouse',
    'human_ortholog_symbol',       # null: Phase 3
    'human_ortholog_entrez',       # null: Phase 3
    'ortholog_mapping_status',     # null: Phase 3
    'mapping_status',
    'mapping_notes',
]

# Regex pattern for isoform-level UniProt accessions (e.g. Q8K1M6-14)
ISOFORM_RE = re.compile(r'^([A-Z0-9]+)-\d+$', re.IGNORECASE)


# ============================================================
# CLI
# ============================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description='Processes 4.5-4.9: UniProt bulk ID mapping for ProSIFT.'
    )
    parser.add_argument(
        '--matrix', '-m', required=True,
        help='Path to {run_id}.filtered_matrix.parquet'
    )
    parser.add_argument(
        '--params', '-p', required=True,
        help='Path to {run_id}_params.yml'
    )
    parser.add_argument(
        '--run_id', '-r', required=True,
        help='Run identifier (used to name output files)'
    )
    parser.add_argument(
        '--cachedir', '-d', default='cache/id_mapping',
        help='Directory for UniProt mapping cache (default: cache/id_mapping)'
    )
    parser.add_argument(
        '--outdir', '-o', default='.',
        help='Output directory (default: current directory)'
    )
    return parser.parse_args()


# ============================================================
# Load inputs
# ============================================================

def load_params(params_path: str) -> dict:
    '''Load and return the ProSIFT params YAML.'''
    with open(params_path, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)


def extract_protein_ids(matrix_path: str, id_col: str) -> list[str]:
    '''
    Read filtered_matrix.parquet and return a list of protein IDs in
    the order they appear. The matrix already has first-accession extraction
    applied (from validate_inputs.py), so no further parsing is needed.
    '''
    df = pd.read_parquet(matrix_path, columns=[id_col])
    return df[id_col].tolist()


# ============================================================
# Cache management
# ============================================================

def _cache_paths(cachedir: str, cache_key: str) -> tuple[str, str]:
    parquet = os.path.join(cachedir, f'{cache_key}.parquet')
    meta    = os.path.join(cachedir, f'{cache_key}.json')
    return parquet, meta


def compute_cache_key(protein_ids: list[str], organism: str) -> str:
    '''MD5 of sorted IDs + source type + organism.'''
    payload = '|'.join(sorted(protein_ids)) + f'|{FROM_TYPE}|{organism}'
    return hashlib.md5(payload.encode('utf-8')).hexdigest()


def load_cache(cachedir: str, cache_key: str, cache_days: int) -> Optional[pd.DataFrame]:
    '''Return cached DataFrame if valid (exists and not expired), else None.'''
    parquet_path, meta_path = _cache_paths(cachedir, cache_key)
    if not os.path.isfile(parquet_path) or not os.path.isfile(meta_path):
        return None
    with open(meta_path, 'r') as f:
        meta = json.load(f)
    age_days = (datetime.now(timezone.utc).timestamp() - meta['timestamp']) / 86400
    if age_days > cache_days:
        print(f'  Cache expired ({age_days:.1f} days old, limit {cache_days} days).')
        return None
    print(f'  Cache hit ({age_days:.1f} days old). Loading from {parquet_path}')
    return pd.read_parquet(parquet_path)


def save_cache(cachedir: str, cache_key: str, df: pd.DataFrame) -> None:
    '''Write mapping table and metadata to cache directory.'''
    os.makedirs(cachedir, exist_ok=True)
    parquet_path, meta_path = _cache_paths(cachedir, cache_key)
    df.to_parquet(parquet_path, index=False)
    meta = {
        'timestamp': datetime.now(timezone.utc).timestamp(),
        'n_proteins': len(df),
        'cache_key': cache_key,
    }
    with open(meta_path, 'w') as f:
        json.dump(meta, f, indent=2)
    print(f'  Cached to {parquet_path}')


# ============================================================
# UniProt ID mapping API
# ============================================================

def _submit_job(ids: list[str]) -> str:
    '''Submit bulk ID mapping job. Returns job ID.'''
    resp = requests.post(
        IDMAP_RUN_URL,
        data={'ids': ','.join(ids), 'from': FROM_TYPE, 'to': TO_TYPE},
        timeout=60,
    )
    resp.raise_for_status()
    return resp.json()['jobId']


def _get_results_url(job_id: str) -> str:
    '''
    Poll the status endpoint until the job finishes, then return the
    results URL.

    The UniProt API v2 redirects the status endpoint to the results endpoint
    once a job is complete. When `to=UniProtKB`, the redirect goes to:
      .../idmapping/uniprotkb/results/{jobId}
    We follow the redirect and detect completion by URL pattern or inline
    results in the response body.
    '''
    for _ in range(MAX_POLLS):
        resp = requests.get(
            f'{IDMAP_STATUS}/{job_id}',
            timeout=30,
            allow_redirects=True,
        )
        resp.raise_for_status()
        # If we were redirected to a results endpoint, the job is done.
        # The redirect URL IS the results URL.
        if '/results/' in resp.url:
            return resp.url
        data = resp.json()
        status = data.get('jobStatus', '')
        if status == 'FINISHED':
            # No redirect yet -- construct the results URL from convention
            return f'{UNIPROT_BASE}/idmapping/uniprotkb/results/{job_id}'
        if status in ('ERROR', 'FAILED'):
            raise RuntimeError(f'UniProt ID mapping job failed: {data}')
        time.sleep(POLL_INTERVAL)
    raise TimeoutError(
        f'UniProt ID mapping job {job_id} did not finish after '
        f'{MAX_POLLS * POLL_INTERVAL}s'
    )


def _fetch_all_results(results_url: str) -> list[dict]:
    '''
    Paginate through all results at the given URL and return a flat list
    of {"from": ..., "to": {...}} dicts. Pagination is driven by the
    Link: <url>; rel="next" header.
    '''
    all_results: list[dict] = []
    next_url: Optional[str] = results_url
    page = 0

    while next_url:
        page += 1
        sep = '&' if '?' in next_url else '?'
        paged_url = f'{next_url}{sep}size={RESULTS_PAGE_SIZE}'
        resp = requests.get(paged_url, timeout=120)
        resp.raise_for_status()
        data = resp.json()
        batch = data.get('results', [])
        all_results.extend(batch)
        print(f'    Page {page}: {len(batch)} results (total so far: {len(all_results)})')

        # Follow Link: <...>; rel="next" header
        link_header = resp.headers.get('Link', '')
        next_url = _parse_next_link(link_header)

    return all_results


def _parse_next_link(link_header: str) -> Optional[str]:
    '''Extract the "next" URL from a Link header, or None.'''
    # Format: <https://...>; rel="next"
    for part in link_header.split(','):
        part = part.strip()
        if 'rel="next"' in part:
            # extract URL between < and >
            start = part.find('<')
            end   = part.find('>')
            if start >= 0 and end > start:
                return part[start + 1:end]
    return None


def query_uniprot(ids: list[str]) -> list[dict]:
    '''
    Submit a bulk ID mapping job, poll until done, paginate through all
    results, and return a flat list of raw result dicts.
    Retries the full sequence up to 3 times on network errors.
    '''
    last_exc: Optional[Exception] = None
    for attempt, delay in enumerate([0] + RETRY_DELAYS):
        if delay:
            print(f'  Retrying in {delay}s (attempt {attempt + 1})...')
            time.sleep(delay)
        try:
            print(f'  Submitting {len(ids)} IDs to UniProt ID mapping API...')
            job_id = _submit_job(ids)
            print(f'  Job ID: {job_id} -- polling for completion...')
            results_url = _get_results_url(job_id)
            print(f'  Job finished. Fetching results from: {results_url}')
            return _fetch_all_results(results_url)
        except (requests.RequestException, RuntimeError, TimeoutError) as exc:
            last_exc = exc
            print(f'  Warning: attempt {attempt + 1} failed -- {exc}')
    raise RuntimeError(
        f'UniProt ID mapping failed after {len(RETRY_DELAYS) + 1} attempts. '
        f'Last error: {last_exc}'
    )


# ============================================================
# Response parsing
# ============================================================

def _is_isoform(accession: str) -> bool:
    return bool(ISOFORM_RE.match(accession))


def _canonical(accession: str) -> str:
    '''Strip isoform suffix if present.'''
    m = ISOFORM_RE.match(accession)
    return m.group(1) if m else accession


def _extract_fields(to_entry: dict) -> tuple[Optional[str], Optional[str], Optional[str], Optional[str]]:
    '''
    Extract (canonical_accession, gene_symbol, entrez_id, ensembl_gene)
    from a UniProt KB entry JSON object (the "to" field in API results).

    Gene symbol: to.genes[0].geneName.value
    Entrez ID:   first cross-reference with database == "GeneID", .id field
    Ensembl gene: first cross-reference with database == "Ensembl", property
                  with key == "GeneId" (note: GeneId not GeneID)
    '''
    accession = to_entry.get('primaryAccession')

    # Gene symbol
    genes = to_entry.get('genes', [])
    gene_symbol: Optional[str] = None
    if genes:
        gene_name = genes[0].get('geneName', {})
        gene_symbol = gene_name.get('value') or None

    # Entrez ID and Ensembl gene from cross-references
    entrez_id:    Optional[str] = None
    ensembl_gene: Optional[str] = None
    for xref in to_entry.get('uniProtKBCrossReferences', []):
        db = xref.get('database', '')
        if db == 'GeneID' and entrez_id is None:
            entrez_id = xref.get('id') or None
        elif db == 'Ensembl' and ensembl_gene is None:
            for prop in xref.get('properties', []):
                if prop.get('key') == 'GeneId':
                    ensembl_gene = prop.get('value') or None
                    break
        if entrez_id and ensembl_gene:
            break  # both found, no need to scan further

    return accession, gene_symbol, entrez_id, ensembl_gene


def parse_json_results(api_results: list[dict], input_ids: list[str]) -> pd.DataFrame:
    '''
    Parse the list of {"from": ..., "to": {...}} dicts returned by the
    UniProt ID mapping API and build the mapping table.

    Steps:
    1. Group API results by input ID.
    2. Classify mapping_status: mapped, unmapped, isoform_collapsed,
       accession_redirected, multiple_mappings.
    3. Add null ortholog columns (Phase 3).
    4. Build one row per input ID preserving input order.
    '''
    # Build a lookup: input_id -> list of "to" entry dicts
    results_by_id: dict[str, list[dict]] = {}
    for item in api_results:
        from_id  = item.get('from', '')
        to_entry = item.get('to', {})
        results_by_id.setdefault(from_id, []).append(to_entry)

    # --- Build one row per input ID ---
    rows: list[dict] = []

    for input_id in input_ids:
        hits = results_by_id.get(input_id, [])

        if len(hits) == 0:
            # --- Unmapped ---
            rows.append({
                'input_id':                 input_id,
                'uniprot_accession':        None,
                'gene_symbol_mouse':        None,
                'entrez_id_mouse':          None,
                'ensembl_gene_mouse':       None,
                'human_ortholog_symbol':    None,
                'human_ortholog_entrez':    None,
                'ortholog_mapping_status':  None,
                'mapping_status':           'unmapped',
                'mapping_notes':            'ID not found in UniProt response.',
            })

        elif len(hits) == 1:
            accession, gene_sym, entrez, ensembl = _extract_fields(hits[0])
            is_isoform = _is_isoform(input_id) and accession == _canonical(input_id)
            is_redirected = not is_isoform and accession != input_id

            if is_isoform:
                status = 'isoform_collapsed'
                notes  = f'Isoform {input_id} collapsed to canonical {accession}.'
            elif is_redirected:
                # Input was a secondary/merged accession; UniProt returned a
                # different primary. The abundance data is still keyed on input_id
                # but all annotation columns use the canonical accession.
                status = 'accession_redirected'
                notes  = (
                    f'Input accession {input_id} is a secondary or merged '
                    f'accession; UniProt canonical is {accession}. '
                    f'Abundance data remains keyed on input_id.'
                )
            else:
                status = 'mapped'
                notes  = ''
            rows.append({
                'input_id':                 input_id,
                'uniprot_accession':        accession,
                'gene_symbol_mouse':        gene_sym,
                'entrez_id_mouse':          entrez,
                'ensembl_gene_mouse':       ensembl,
                'human_ortholog_symbol':    None,
                'human_ortholog_entrez':    None,
                'ortholog_mapping_status':  None,
                'mapping_status':           status,
                'mapping_notes':            notes,
            })

        else:
            # --- Multiple mappings: take first, log the rest ---
            accession, gene_sym, entrez, ensembl = _extract_fields(hits[0])
            other_accessions = '; '.join(
                h.get('primaryAccession', '') for h in hits[1:]
            )
            rows.append({
                'input_id':                 input_id,
                'uniprot_accession':        accession,
                'gene_symbol_mouse':        gene_sym,
                'entrez_id_mouse':          entrez,
                'ensembl_gene_mouse':       ensembl,
                'human_ortholog_symbol':    None,
                'human_ortholog_entrez':    None,
                'ortholog_mapping_status':  None,
                'mapping_status':           'multiple_mappings',
                'mapping_notes':            (
                    f'Input ID mapped to multiple UniProt entries. '
                    f'Using first ({accession}); '
                    f'others: {other_accessions}.'
                ),
            })

    df = pd.DataFrame(rows, columns=SCHEMA_COLS)
    return df


# ============================================================
# Validation report (part 3)
# ============================================================

def write_report(
    filepath: str,
    run_id: str,
    n_input: int,
    df: pd.DataFrame,
    cache_hit: bool,
    cache_key: str,
    cache_days: int,
) -> None:
    '''Write validation_report_part3.txt summarizing ID mapping results.'''

    status_counts = df['mapping_status'].value_counts()

    def _pct(n: int) -> str:
        return f'{100 * n / n_input:.1f}%'

    n_mapped      = int(status_counts.get('mapped', 0))
    n_isoform     = int(status_counts.get('isoform_collapsed', 0))
    n_redirected  = int(status_counts.get('accession_redirected', 0))
    n_multi       = int(status_counts.get('multiple_mappings', 0))
    n_unmapped    = int(status_counts.get('unmapped', 0))

    # Coverage of individual annotation fields (excluding unmapped rows)
    mapped_df = df[df['mapping_status'] != 'unmapped']
    n_gene_sym  = int(mapped_df['gene_symbol_mouse'].notna().sum())
    n_entrez    = int(mapped_df['entrez_id_mouse'].notna().sum())
    n_ensembl   = int(mapped_df['ensembl_gene_mouse'].notna().sum())

    unmapped_ids = df.loc[df['mapping_status'] == 'unmapped', 'input_id'].tolist()
    unmapped_preview = ', '.join(unmapped_ids[:10])
    if len(unmapped_ids) > 10:
        unmapped_preview += f' ... ({len(unmapped_ids) - 10} more)'

    lines = [
        '=' * 70,
        f'  ProSIFT Validation Report -- Part 3: UniProt ID Mapping',
        f'  Run: {run_id}',
        f'  Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}',
        '=' * 70,
        '',
        '--- Mapping Summary ---',
        '',
        f'  Proteins submitted to UniProt:    {n_input}',
        f'  Cache:                            {"HIT" if cache_hit else "MISS"} '
                                               f'(key: {cache_key[:12]}...)',
        '',
        '  Mapping status:',
        f'    Mapped:                         {n_mapped} ({_pct(n_mapped)})',
        f'    Isoform collapsed:              {n_isoform} ({_pct(n_isoform)})',
        f'    Accession redirected:           {n_redirected} ({_pct(n_redirected)})',
        f'    Multiple mappings:              {n_multi} ({_pct(n_multi)})',
        f'    Unmapped:                       {n_unmapped} ({_pct(n_unmapped)})',
        '',
        '--- Annotation Field Coverage ---',
        '  (excludes unmapped proteins)',
        '',
        f'  Gene symbol (mouse):              {n_gene_sym} / {n_input - n_unmapped}',
        f'  Entrez ID (mouse):                {n_entrez} / {n_input - n_unmapped}',
        f'  Ensembl gene (mouse):             {n_ensembl} / {n_input - n_unmapped}',
        f'  (Denominator excludes {n_unmapped} unmapped protein(s) with no annotation.)',
        '',
        '  Ortholog columns:                 null (Phase 3 -- not yet implemented)',
        '',
    ]

    if n_unmapped > 0:
        lines += [
            '--- Unmapped Protein IDs ---',
            '',
            f'  {n_unmapped} protein(s) could not be mapped. These are retained in',
            '  the mapping table with null annotation fields.',
            '  Common causes: obsolete accessions, contaminant prefixes,',
            '  non-UniProt IDs.',
            '',
            f'  Unmapped IDs: {unmapped_preview}',
            '',
        ]

    if n_redirected > 0:
        redirected = df.loc[df['mapping_status'] == 'accession_redirected',
                            ['input_id', 'uniprot_accession']]
        redirect_lines = [
            f'    {row.input_id}  ->  {row.uniprot_accession}'
            for row in redirected.itertuples()
        ]
        lines += [
            '--- Accession-Redirected Protein IDs ---',
            '',
            f'  {n_redirected} protein(s) submitted as secondary or merged accessions.',
            '  UniProt returned a different canonical accession for each. Abundance',
            '  data remains keyed on input_id; downstream annotation uses uniprot_accession.',
            '  This typically indicates the search database FASTA was built from an',
            '  older UniProt release.',
            '',
        ] + redirect_lines + ['']

    if n_multi > 0:
        multi_ids = df.loc[
            df['mapping_status'] == 'multiple_mappings', 'input_id'
        ].tolist()
        lines += [
            '--- Multiple-Mapping Protein IDs ---',
            '',
            f'  {n_multi} protein(s) returned more than one UniProt entry.',
            '  The first entry was used; see mapping_notes for alternatives.',
            f'  IDs: {", ".join(multi_ids[:10])}',
            '',
        ]

    lines += [
        '--- Cache Configuration ---',
        '',
        f'  Cache lifetime:                   {cache_days} days',
        f'  Cache key:                        {cache_key}',
        '',
        '--- Notes ---',
        '',
        '  The ortholog mapping columns (human_ortholog_symbol,',
        '  human_ortholog_entrez, ortholog_mapping_status) are reserved',
        '  for Phase 3. All values are null in this run.',
        '',
        '  Downstream modules that use human-centric databases (DGIdb,',
        '  DisGeNET) will require ortholog mapping to be implemented.',
        '=' * 70,
    ]

    with open(filepath, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines) + '\n')


# ============================================================
# Main
# ============================================================

def main() -> None:
    args = parse_args()

    run_id    = args.run_id
    outdir    = args.outdir
    cachedir  = args.cachedir

    os.makedirs(outdir, exist_ok=True)

    # ----------------------------------------------------------
    # 1. Load params and extract settings
    # ----------------------------------------------------------
    print(f'[{run_id}] Loading params...')
    params = load_params(args.params)

    id_col      = params['input']['protein_id_column']
    organism    = params['project']['organism']
    cache_days  = params.get('databases', {}).get('cache_days', 30)

    # ----------------------------------------------------------
    # 2. Extract protein IDs from filtered matrix
    # ----------------------------------------------------------
    print(f'[{run_id}] Reading filtered matrix: {args.matrix}')
    protein_ids = extract_protein_ids(args.matrix, id_col)
    print(f'  {len(protein_ids)} proteins')

    if len(protein_ids) == 0:
        print('Error: filtered matrix contains no proteins.')
        sys.exit(1)

    # ----------------------------------------------------------
    # 3. Check cache
    # ----------------------------------------------------------
    cache_key = compute_cache_key(protein_ids, organism)
    print(f'[{run_id}] Cache key: {cache_key}')

    mapping_df = load_cache(cachedir, cache_key, cache_days)
    cache_hit  = mapping_df is not None

    # ----------------------------------------------------------
    # 4. Query UniProt if cache miss
    # ----------------------------------------------------------
    if not cache_hit:
        print(f'[{run_id}] Cache miss -- querying UniProt ID mapping API...')
        api_results = query_uniprot(protein_ids)
        print(f'  Total results returned: {len(api_results)}')
        mapping_df = parse_json_results(api_results, protein_ids)

        # ----------------------------------------------------------
        # 5. Cache the result
        # ----------------------------------------------------------
        print(f'[{run_id}] Saving mapping table to cache...')
        save_cache(cachedir, cache_key, mapping_df)
    else:
        # Re-order cached table to match current protein_id order (may differ
        # if proteins were re-ordered between runs with the same set of IDs)
        id_order = {pid: i for i, pid in enumerate(protein_ids)}
        mapping_df = mapping_df.copy()
        mapping_df['_order'] = mapping_df['input_id'].map(id_order)
        mapping_df = mapping_df.sort_values('_order').drop(columns='_order')

    # ----------------------------------------------------------
    # 6. Write outputs
    # ----------------------------------------------------------
    mapping_out = os.path.join(outdir, f'{run_id}.id_mapping.parquet')
    report_out  = os.path.join(outdir, f'{run_id}.validation_report_part3.txt')

    mapping_df.to_parquet(mapping_out, index=False)
    print(f'[{run_id}] Written: {mapping_out}')

    write_report(
        report_out,
        run_id=run_id,
        n_input=len(protein_ids),
        df=mapping_df,
        cache_hit=cache_hit,
        cache_key=cache_key,
        cache_days=cache_days,
    )
    print(f'[{run_id}] Written: {report_out}')

    # Summary
    status_counts = mapping_df['mapping_status'].value_counts()
    print(f'\n[{run_id}] Mapping complete:')
    for status, count in status_counts.items():
        print(f'  {status}: {count} ({100 * count / len(protein_ids):.1f}%)')

    print(f'\n[{run_id}] Done.')


if __name__ == '__main__':
    main()

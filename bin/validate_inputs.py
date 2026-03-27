#!/usr/bin/env python3
# title: validate_inputs.py
# project: ProSIFT
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-03-26
# last modified: 2026-03-26
#
# purpose:
#   Module 01, Processes 4.1-4.3. Validates the abundance matrix and metadata
#   file for a single ProSIFT run, then cross-validates sample IDs between them.
#   Produces validated Parquet outputs containing only the matched working set of
#   samples, and a partial validation report (to be concatenated with reports from
#   subsequent processes).
#
# inputs:
#   --abundance  : master abundance CSV (protein_id + abundance + peptide_count cols)
#   --metadata   : per-run metadata CSV (sample_id + group column)
#   --params     : per-run params.yml
#   --run-id     : string identifier used to name output files
#   --outdir     : output directory (default: current working directory)
#
# outputs:
#   {outdir}/{run_id}.validated_matrix.parquet
#   {outdir}/{run_id}.validated_metadata.parquet
#   {outdir}/{run_id}.validation_report_part1.txt
#
# usage example:
#   python bin/validate_inputs.py \
#       --abundance prosift_inputs/CTX_synaptosome_abundance.csv \
#       --metadata  prosift_inputs/CTXcyto_WT_vs_CTXcyto_KO/CTXcyto_WT_vs_CTXcyto_KO_metadata.csv \
#       --params    prosift_inputs/CTXcyto_WT_vs_CTXcyto_KO/CTXcyto_WT_vs_CTXcyto_KO_params.yml \
#       --run-id    CTXcyto_WT_vs_CTXcyto_KO \
#       --outdir    results/CTXcyto_WT_vs_CTXcyto_KO/validation

import argparse
import sys
import os
from datetime import datetime
from io import StringIO

import pandas as pd
import yaml


# ============================================================
# Argument parsing
# ============================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description='ProSIFT Module 01 -- validate abundance matrix and metadata.'
    )
    parser.add_argument('--abundance', required=True,
                        help='Path to master abundance CSV/TSV')
    parser.add_argument('--metadata', required=True,
                        help='Path to per-run metadata CSV/TSV')
    parser.add_argument('--params', required=True,
                        help='Path to per-run params.yml')
    parser.add_argument('--run-id', required=True, dest='run_id',
                        help='Run identifier for output file naming')
    parser.add_argument('--outdir', default='.',
                        help='Output directory (default: current directory)')
    return parser.parse_args()


# ============================================================
# Config loading
# ============================================================

def load_params(params_path: str) -> dict:
    '''Load and return params.yml as a dict.'''
    with open(params_path, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)


def infer_separator(filepath: str) -> str:
    '''Return "," for .csv, "\t" for .tsv/.txt, else try to sniff.'''
    ext = os.path.splitext(filepath)[1].lower()
    if ext == '.csv':
        return ','
    if ext in ('.tsv', '.txt'):
        return '\t'
    # sniff first line
    with open(filepath, 'r', encoding='utf-8-sig') as f:
        first = f.readline()
    return '\t' if first.count('\t') > first.count(',') else ','


# ============================================================
# Report helpers
# ============================================================

class Report:
    '''Accumulates validation report lines; writes to file at the end.'''

    def __init__(self, run_id: str):
        self.run_id = run_id
        self._lines: list[str] = []
        self._warnings: list[str] = []

    def section(self, title: str) -> None:
        self._lines.append('')
        self._lines.append('=' * 64)
        self._lines.append(title)
        self._lines.append('=' * 64)

    def line(self, text: str = '') -> None:
        self._lines.append(text)

    def warn(self, text: str) -> None:
        msg = f'  WARNING: {text}'
        self._lines.append(msg)
        self._warnings.append(text)

    def ok(self, text: str) -> None:
        self._lines.append(f'  OK: {text}')

    def info(self, text: str) -> None:
        self._lines.append(f'  {text}')

    def error_exit(self, text: str) -> None:
        '''Log the error to the report then exit non-zero.'''
        self._lines.append(f'  ERROR: {text}')
        self._lines.append('')
        self._lines.append('Validation failed. Pipeline stopped.')
        # write partial report before exiting so the user has context
        print(f'ERROR [{self.run_id}]: {text}', file=sys.stderr)
        raise SystemExit(1)

    def write(self, outpath: str) -> None:
        with open(outpath, 'w', encoding='utf-8') as f:
            f.write('\n'.join(self._lines))
            f.write('\n')

    @property
    def warning_count(self) -> int:
        return len(self._warnings)


# ============================================================
# Process 4.1 -- Validate Matrix
# ============================================================

def validate_matrix(
    abundance_path: str,
    params: dict,
    report: Report
) -> pd.DataFrame:
    '''
    Read and validate the abundance matrix.

    Returns a DataFrame with the protein_id column plus all abundance and
    peptide_count columns. The protein_id column contains the representative
    accession (first token of semicolon-delimited groups).

    Raises SystemExit on any hard validation failure.
    '''
    report.section('PROCESS 4.1 -- VALIDATE MATRIX')

    # --- Read file ---
    sep = infer_separator(abundance_path)
    try:
        df = pd.read_csv(abundance_path, sep=sep, encoding='utf-8-sig',
                         low_memory=False)
    except Exception as exc:
        report.error_exit(f'Cannot parse abundance file: {exc}')

    report.ok(f'File parsed: {os.path.basename(abundance_path)}')
    report.info(f'{len(df)} rows, {len(df.columns)} columns')

    # --- Protein ID column ---
    id_col: str = params['input']['protein_id_column']
    if id_col not in df.columns:
        report.error_exit(
            f'protein_id_column "{id_col}" not found in abundance file. '
            f'Columns present: {list(df.columns)}'
        )
    report.ok(f'Protein ID column found: "{id_col}"')

    # --- Extract representative accession from semicolon-delimited groups ---
    original_ids = df[id_col].astype(str)
    df[id_col] = original_ids.str.split(';').str[0].str.strip()
    multi_member = (original_ids.str.contains(';', na=False)).sum()
    if multi_member > 0:
        report.info(
            f'{multi_member} protein group rows with semicolon-delimited '
            f'accessions -- first accession used as representative ID.'
        )

    # --- Duplicate protein IDs ---
    dupes = df[id_col].duplicated(keep=False)
    if dupes.any():
        dupe_ids = df.loc[dupes, id_col].unique().tolist()[:10]
        report.error_exit(
            f'{dupes.sum()} duplicate protein IDs after representative '
            f'accession extraction. Examples: {dupe_ids}'
        )
    report.ok(f'{len(df)} unique protein IDs')

    # --- Identify abundance and peptide count columns ---
    abund_prefix: str = params['input'].get('abundance_prefix', '')
    pep_prefix: str | None = params['input'].get('peptide_count_prefix') or None

    non_id_cols = [c for c in df.columns if c != id_col]

    if pep_prefix:
        pep_cols = [c for c in non_id_cols if c.startswith(pep_prefix)]
        abund_cols = [c for c in non_id_cols if c.startswith(abund_prefix)
                      and not c.startswith(pep_prefix)]
    else:
        pep_cols = []
        abund_cols = [c for c in non_id_cols if c.startswith(abund_prefix)]

    if not abund_cols:
        report.error_exit(
            f'No abundance columns found with prefix "{abund_prefix}". '
            f'Non-ID columns present: {non_id_cols[:10]}'
        )
    report.ok(f'{len(abund_cols)} abundance columns identified '
              f'(prefix: "{abund_prefix}")')

    # --- Abundance columns: numeric or NA ---
    bad_abund: list[str] = []
    for col in abund_cols:
        coerced = pd.to_numeric(df[col], errors='coerce')
        # if coercion introduced NAs that weren't already NA, the column has
        # non-numeric, non-empty values
        original_na = df[col].isna()
        new_na = coerced.isna()
        spurious_na = new_na & ~original_na
        if spurious_na.any():
            bad_abund.append(col)
    if bad_abund:
        report.error_exit(
            f'{len(bad_abund)} abundance column(s) contain non-numeric values: '
            f'{bad_abund[:5]}'
        )
    report.ok('All abundance columns are numeric or NA')

    # --- Convert abundance columns to float ---
    df[abund_cols] = df[abund_cols].apply(pd.to_numeric, errors='coerce')

    # --- Completely empty rows (all abundance values NA) ---
    all_na_rows = df[abund_cols].isna().all(axis=1)
    if all_na_rows.any():
        report.warn(
            f'{all_na_rows.sum()} protein row(s) have no abundance values '
            f'across any sample -- these rows carry no information.'
        )

    # --- Completely empty columns ---
    all_na_cols = [c for c in abund_cols if df[c].isna().all()]
    if all_na_cols:
        report.warn(
            f'{len(all_na_cols)} abundance column(s) are entirely NA: '
            f'{all_na_cols}'
        )

    # --- Minimum protein count ---
    if len(df) < 2:
        report.error_exit(
            f'Abundance matrix contains fewer than 2 proteins ({len(df)} rows). '
            f'Cannot proceed.'
        )

    # --- Peptide count columns (if specified) ---
    if pep_prefix:
        if not pep_cols:
            report.error_exit(
                f'peptide_count_prefix "{pep_prefix}" is set but no matching '
                f'columns found in abundance file.'
            )
        report.ok(f'{len(pep_cols)} peptide count columns identified '
                  f'(prefix: "{pep_prefix}")')

        # sample IDs embedded in peptide count column names must match
        # those embedded in abundance column names
        abund_sample_ids = {c[len(abund_prefix):] for c in abund_cols}
        pep_sample_ids = {c[len(pep_prefix):] for c in pep_cols}
        missing_in_pep = abund_sample_ids - pep_sample_ids
        extra_in_pep = pep_sample_ids - abund_sample_ids
        if missing_in_pep:
            report.error_exit(
                f'Peptide count columns missing for sample(s) present in '
                f'abundance columns: {sorted(missing_in_pep)[:10]}'
            )
        if extra_in_pep:
            report.warn(
                f'Peptide count columns exist for sample(s) not present in '
                f'abundance columns (will be ignored): {sorted(extra_in_pep)[:10]}'
            )

        # peptide count columns: numeric (int or float), zero, or NA
        bad_pep: list[str] = []
        for col in pep_cols:
            coerced = pd.to_numeric(df[col], errors='coerce')
            original_na = df[col].isna()
            new_na = coerced.isna()
            spurious_na = new_na & ~original_na
            if spurious_na.any():
                bad_pep.append(col)
        if bad_pep:
            report.error_exit(
                f'{len(bad_pep)} peptide count column(s) contain non-numeric '
                f'values: {bad_pep[:5]}'
            )
        df[pep_cols] = df[pep_cols].apply(pd.to_numeric, errors='coerce')
        report.ok('All peptide count columns are numeric or NA')
    else:
        report.info('No peptide_count_prefix specified -- DEqMS weighting '
                    'will not be available; falling back to limma.')

    # --- Return subset of columns only ---
    keep_cols = [id_col] + abund_cols + pep_cols
    df = df[keep_cols].copy()

    report.line()
    report.info(f'Matrix validation complete: {len(df)} proteins, '
                f'{len(abund_cols)} abundance columns'
                + (f', {len(pep_cols)} peptide count columns' if pep_cols else ''))
    return df


# ============================================================
# Process 4.2 -- Validate Metadata
# ============================================================

def validate_metadata(
    metadata_path: str,
    params: dict,
    report: Report
) -> pd.DataFrame:
    '''
    Read and validate the metadata file.

    Returns a DataFrame with at least sample_id and the group column.
    Raises SystemExit on any hard validation failure.
    '''
    report.section('PROCESS 4.2 -- VALIDATE METADATA')

    sep = infer_separator(metadata_path)
    try:
        meta = pd.read_csv(metadata_path, sep=sep, encoding='utf-8-sig')
    except Exception as exc:
        report.error_exit(f'Cannot parse metadata file: {exc}')

    report.ok(f'File parsed: {os.path.basename(metadata_path)}')
    report.info(f'{len(meta)} rows, {len(meta.columns)} columns')

    # --- sample_id column ---
    # The metadata must have a sample_id column. Its name is always "sample_id"
    # (prepare_prosift_input.py writes it that way).
    if 'sample_id' not in meta.columns:
        report.error_exit(
            f'"sample_id" column not found in metadata. '
            f'Columns present: {list(meta.columns)}'
        )

    # --- Duplicate sample IDs ---
    dupes = meta['sample_id'].duplicated(keep=False)
    if dupes.any():
        report.error_exit(
            f'{dupes.sum()} duplicate sample_id values in metadata: '
            f'{meta.loc[dupes, "sample_id"].tolist()}'
        )
    report.ok(f'{len(meta)} unique sample IDs')

    # --- Group column ---
    group_col: str = params['design']['group_column']
    if group_col not in meta.columns:
        report.error_exit(
            f'group_column "{group_col}" not found in metadata. '
            f'Columns present: {list(meta.columns)}'
        )
    groups = meta[group_col].dropna().unique()
    if len(groups) < 2:
        report.error_exit(
            f'group_column "{group_col}" contains fewer than 2 distinct groups '
            f'(found: {list(groups)}). Cannot run a comparison.'
        )
    report.ok(f'Group column "{group_col}" found: {sorted(str(g) for g in groups)}')

    # --- Covariate columns ---
    covariates: list[str] = params['design'].get('covariates') or []
    for cov in covariates:
        if cov not in meta.columns:
            report.error_exit(
                f'Covariate column "{cov}" specified in params but not found '
                f'in metadata. Columns present: {list(meta.columns)}'
            )
    if covariates:
        report.ok(f'Covariate column(s) found: {covariates}')

    # --- Batch column ---
    batch_col: str | None = params['design'].get('batch_column') or None
    if batch_col:
        if batch_col not in meta.columns:
            report.error_exit(
                f'batch_column "{batch_col}" specified in params but not found '
                f'in metadata. Columns present: {list(meta.columns)}'
            )
        report.ok(f'Batch column "{batch_col}" found')

    report.line()
    report.info(f'Metadata validation complete: {len(meta)} samples, '
                f'{len(groups)} groups')
    return meta


# ============================================================
# Process 4.3 -- Cross-Validate Samples
# ============================================================

def cross_validate(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    params: dict,
    report: Report
) -> tuple[pd.DataFrame, pd.DataFrame]:
    '''
    Match sample IDs between the abundance matrix and the metadata.

    Returns (matrix_subset, metadata_subset) containing only the working set
    of matched samples. Raises SystemExit if no samples match or if any group
    falls below qc.min_samples_per_group.
    '''
    report.section('PROCESS 4.3 -- CROSS-VALIDATE SAMPLES')

    id_col: str = params['input']['protein_id_column']
    abund_prefix: str = params['input'].get('abundance_prefix', '')
    pep_prefix: str | None = params['input'].get('peptide_count_prefix') or None
    group_col: str = params['design']['group_column']
    min_per_group: int = params['qc'].get('min_samples_per_group', 2)

    # Derive sample IDs from abundance column names by stripping the prefix
    abund_cols = [c for c in matrix.columns if c != id_col
                  and (not pep_prefix or not c.startswith(pep_prefix))]
    matrix_sample_ids: set[str] = {c[len(abund_prefix):] for c in abund_cols}
    meta_sample_ids: set[str] = set(metadata['sample_id'].astype(str))

    matched = matrix_sample_ids & meta_sample_ids
    only_in_matrix = matrix_sample_ids - meta_sample_ids
    only_in_meta = meta_sample_ids - matrix_sample_ids

    report.info(f'Abundance matrix samples: {len(matrix_sample_ids)}')
    report.info(f'Metadata samples:         {len(meta_sample_ids)}')
    report.info(f'Matched (working set):    {len(matched)}')

    if only_in_matrix:
        report.warn(
            f'{len(only_in_matrix)} sample(s) in abundance matrix but not '
            f'in metadata -- dropped from analysis: '
            f'{sorted(only_in_matrix)}'
        )
    if only_in_meta:
        report.warn(
            f'{len(only_in_meta)} sample(s) in metadata but not in '
            f'abundance matrix -- check for typos or missing data: '
            f'{sorted(only_in_meta)}'
        )

    if not matched:
        report.error_exit(
            'No sample IDs match between the abundance matrix and metadata. '
            'Check that sample IDs in the metadata match the column names in '
            'the abundance matrix (after stripping abundance_prefix).'
        )

    # --- Subset matrix to matched samples ---
    keep_abund_cols = [c for c in abund_cols
                       if c[len(abund_prefix):] in matched]
    if pep_prefix:
        pep_cols = [c for c in matrix.columns if c.startswith(pep_prefix)]
        keep_pep_cols = [c for c in pep_cols
                         if c[len(pep_prefix):] in matched]
    else:
        keep_pep_cols = []

    matrix_out = matrix[[id_col] + keep_abund_cols + keep_pep_cols].copy()

    # --- Subset metadata to matched samples ---
    meta_out = metadata[metadata['sample_id'].astype(str).isin(matched)].copy()

    # --- Check group sizes ---
    group_counts = meta_out.groupby(group_col)['sample_id'].count()
    report.line()
    report.info(f'Samples per group (working set):')
    for grp, n in group_counts.items():
        report.info(f'  {grp}: {n}')

    undersized = group_counts[group_counts < min_per_group]
    if not undersized.empty:
        report.error_exit(
            f'The following group(s) have fewer than {min_per_group} samples '
            f'after cross-validation: '
            f'{dict(undersized)}. '
            f'Adjust qc.min_samples_per_group or provide more samples.'
        )

    report.line()
    report.ok(f'Cross-validation complete. Working set: {len(matched)} samples, '
              f'{matrix_out.shape[0]} proteins.')
    return matrix_out, meta_out


# ============================================================
# Output writing
# ============================================================

def write_outputs(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    report: Report,
    run_id: str,
    outdir: str
) -> None:
    os.makedirs(outdir, exist_ok=True)

    matrix_path = os.path.join(outdir, f'{run_id}.validated_matrix.parquet')
    meta_path = os.path.join(outdir, f'{run_id}.validated_metadata.parquet')
    report_path = os.path.join(outdir, f'{run_id}.validation_report_part1.txt')

    matrix.to_parquet(matrix_path, index=False)
    metadata.to_parquet(meta_path, index=False)

    # finalise report header before writing
    report.write(report_path)

    print(f'  Validated matrix:   {matrix_path}')
    print(f'  Validated metadata: {meta_path}')
    print(f'  Validation report:  {report_path}')


# ============================================================
# Main
# ============================================================

def main() -> None:
    args = parse_args()
    params = load_params(args.params)
    run_id = args.run_id

    report = Report(run_id)

    # --- Report header ---
    report.line(f'ProSIFT Validation Report -- {run_id}')
    report.line(f'Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
    report.line(f'Part 1 of 2: Input Validation (Processes 4.1-4.3)')
    report.line(f'  (Part 2: Detection Filter -- see validation_report_part2.txt)')

    print(f'[{run_id}] Validating inputs...')

    # --- Process 4.1: Validate matrix ---
    matrix = validate_matrix(args.abundance, params, report)

    # --- Process 4.2: Validate metadata ---
    metadata = validate_metadata(args.metadata, params, report)

    # --- Process 4.3: Cross-validate ---
    matrix, metadata = cross_validate(matrix, metadata, params, report)

    # --- Summary ---
    report.section('VALIDATION SUMMARY')
    report.info(f'Warnings: {report.warning_count}')
    if report.warning_count == 0:
        report.info('All checks passed with no warnings.')
    report.info(f'Output: {matrix.shape[0]} proteins x '
                f'{matrix.shape[1] - 1} columns (working set)')
    report.info('Status: PASSED -- proceed to detection filter (Part 2)')

    # --- Write outputs ---
    write_outputs(matrix, metadata, report, run_id, args.outdir)
    print(f'[{run_id}] Validation complete. '
          f'{report.warning_count} warning(s). '
          f'{matrix.shape[0]} proteins, {metadata.shape[0]} samples.')


if __name__ == '__main__':
    main()

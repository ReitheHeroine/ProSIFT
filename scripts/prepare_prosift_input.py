# title: prepare_prosift_input.py
# project: ProSIFT
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-03-23
# last modified: 2026-03-23
#
# purpose:
#   Reads a run_config.yml (produced by generate_run_config.py) and the
#   master abundance/metadata files it references. Splits the master data
#   into per-run input files that ProSIFT Module 01 expects:
#     - {run}_abundance.csv   (protein_id + abundance + peptide count columns)
#     - {run}_metadata.csv    (sample_id + group column)
#     - {run}_params.yml      (draft ProSIFT parameters)
#
# inputs:
#   - run_config.yml (references master abundance and metadata CSVs)
#
# outputs:
#   - Per-run subdirectories with the 3 files above
#
# usage example:
#   python scripts/prepare_prosift_input.py --config prosift_inputs/run_config.yml
#   python scripts/prepare_prosift_input.py --config prosift_inputs/run_config.yml --outdir prosift_inputs

import argparse
import csv
import os
import shutil
import sys
from collections import OrderedDict
from datetime import date

try:
    import yaml
except ImportError:
    print('Error: PyYAML is required. Install with: pip install pyyaml')
    sys.exit(1)


# ============================================================
# Loading and validation
# ============================================================

def load_config(config_path):
    '''Load and return the run config YAML.'''
    with open(config_path, 'r', encoding='utf-8') as f:
        config = yaml.safe_load(f)
    return config


def resolve_path(base_dir, filename):
    '''Resolve a filename relative to a base directory.'''
    return os.path.join(base_dir, filename)


def load_csv(filepath):
    '''Load a CSV file and return (headers, rows) where rows is a list of dicts.'''
    with open(filepath, 'r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        headers = reader.fieldnames
        rows = list(reader)
    return headers, rows


def validate_config(config, headers, config_dir):
    '''Validate that the config references valid columns and files.'''
    source = config['source']
    runs = config['runs']

    # check source files exist
    abund_path = resolve_path(config_dir, source['abundance_file'])
    if not os.path.isfile(abund_path):
        print(f'Error: abundance file not found: {abund_path}')
        sys.exit(1)

    # check protein_id column exists
    protein_id_col = source['protein_id_column']
    if protein_id_col not in headers:
        print(f'Error: protein ID column "{protein_id_col}" not found in abundance file.')
        sys.exit(1)

    # check all sample columns exist
    abund_prefix = source.get('abundance_prefix', '')
    pep_prefix = source.get('peptide_count_prefix', None)

    for run_name, run_def in runs.items():
        for sample_id in run_def['samples']:
            abund_col = abund_prefix + sample_id
            if abund_col not in headers:
                print(f'Error: run "{run_name}": abundance column "{abund_col}" not found.')
                sys.exit(1)
            if pep_prefix:
                pep_col = pep_prefix + sample_id
                if pep_col not in headers:
                    print(f'Error: run "{run_name}": peptide count column "{pep_col}" not found.')
                    sys.exit(1)


# ============================================================
# Per-run output writers
# ============================================================

def write_run_abundance(filepath, rows, sample_ids, protein_id_col,
                        abund_prefix, pep_prefix):
    '''
    Write the per-run abundance CSV.
    Columns: protein_id, then bare sample IDs, then peptide_count_ columns.
    '''
    # source column names in the master file
    abund_src_cols = [abund_prefix + sid for sid in sample_ids]
    pep_src_cols = [pep_prefix + sid for sid in sample_ids] if pep_prefix else []

    # output column names
    out_headers = (
        ['protein_id']
        + list(sample_ids)
        + (['peptide_count_' + sid for sid in sample_ids] if pep_prefix else [])
    )

    with open(filepath, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(out_headers)
        for row in rows:
            protein_id = row[protein_id_col]
            abund_vals = [row.get(col, '') for col in abund_src_cols]
            pep_vals = [row.get(col, '') for col in pep_src_cols]
            writer.writerow([protein_id] + abund_vals + pep_vals)


def write_run_metadata(filepath, sample_group_map, group_column):
    '''Write the per-run metadata CSV.'''
    with open(filepath, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['sample_id', group_column])
        for sample_id, group in sample_group_map.items():
            writer.writerow([sample_id, group])


def write_run_params(filepath, run_name, group_column, pep_prefix):
    '''Write a draft params.yml for one run.'''
    today = date.today().isoformat()
    pep_line = f'  peptide_count_prefix: "peptide_count_"' if pep_prefix else f'  peptide_count_prefix: null'

    content = f'''# ProSIFT params -- {run_name}
# Generated by prepare_prosift_input.py on {today}

project:
  name: "{run_name}"
  organism: "mouse"

input:
  abundance_matrix: "{run_name}_abundance.csv"
  metadata: "{run_name}_metadata.csv"
  format: "csv"
  protein_id_column: "protein_id"
  abundance_type: "raw"
  abundance_prefix: ""
{pep_line}

design:
  group_column: "{group_column}"
  covariates: []
  batch_column: null

qc:
  min_samples_per_group: 2
  min_detection_fraction: 0.67

databases:
  cache_days: 30

# --- Modules below not yet implemented ---
# normalization, imputation, differential_abundance,
# enrichment, and database query parameters will be
# added as those modules are designed.
'''
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(content)


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description='Split master ProSIFT input files into per-run directories.'
    )
    parser.add_argument(
        '--config', '-c', required=True,
        help='Path to run_config.yml'
    )
    parser.add_argument(
        '--outdir', '-o', default=None,
        help='Output directory (default: same directory as config file)'
    )
    args = parser.parse_args()

    config_path = args.config
    if not os.path.isfile(config_path):
        print(f'Error: config file not found: {config_path}')
        sys.exit(1)

    config_dir = os.path.dirname(os.path.abspath(config_path))
    outdir = args.outdir if args.outdir else config_dir

    # ----------------------------------------------------------
    # 1. Load config
    # ----------------------------------------------------------
    config = load_config(config_path)
    source = config['source']
    runs = config['runs']

    protein_id_col = source['protein_id_column']
    abund_prefix = source.get('abundance_prefix', '')
    pep_prefix = source.get('peptide_count_prefix', None)

    # ----------------------------------------------------------
    # 2. Load and validate master abundance data
    # ----------------------------------------------------------
    abund_path = resolve_path(config_dir, source['abundance_file'])
    print(f'Loading {os.path.basename(abund_path)}...')
    headers, rows = load_csv(abund_path)
    print(f'  {len(rows)} proteins, {len(headers)} columns')

    validate_config(config, headers, config_dir)
    print(f'  Validation passed.')

    # ----------------------------------------------------------
    # 3. Process each run
    # ----------------------------------------------------------
    print(f'\nGenerating {len(runs)} runs in {outdir}/\n')

    for run_name, run_def in runs.items():
        sample_group_map = OrderedDict(run_def['samples'])
        sample_ids = list(sample_group_map.keys())
        group_column = run_def['group_column']

        # create run directory
        run_dir = os.path.join(outdir, run_name)
        os.makedirs(run_dir, exist_ok=True)

        # write the 3 output files
        write_run_abundance(
            os.path.join(run_dir, f'{run_name}_abundance.csv'),
            rows, sample_ids, protein_id_col, abund_prefix, pep_prefix
        )

        write_run_metadata(
            os.path.join(run_dir, f'{run_name}_metadata.csv'),
            sample_group_map, group_column
        )

        write_run_params(
            os.path.join(run_dir, f'{run_name}_params.yml'),
            run_name, group_column, pep_prefix
        )

        # summary
        groups = OrderedDict()
        for sid, grp in sample_group_map.items():
            groups.setdefault(grp, []).append(sid)
        group_summary = ', '.join(f'{len(v)} {k}' for k, v in groups.items())
        print(f'  {run_name}: {group_summary} ({group_column})')

    # ----------------------------------------------------------
    # 4. Copy run_config.yml for provenance
    # ----------------------------------------------------------
    config_dest = os.path.join(outdir, 'run_config.yml')
    if os.path.abspath(config_path) != os.path.abspath(config_dest):
        shutil.copy2(config_path, config_dest)
        print(f'\nCopied run_config.yml to {outdir}/')

    print(f'\nDone. {len(runs)} runs written to {outdir}/')


if __name__ == '__main__':
    main()

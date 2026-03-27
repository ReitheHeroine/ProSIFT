# title: prepare_prosift_input.py
# project: ProSIFT
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-03-23
# last modified: 2026-03-26  (added samplesheet.csv generation)
#
# purpose:
#   Reads a run_config.yml (produced by generate_run_config.py) and the
#   master abundance file it references. For each run, writes two files into
#   a per-run subdirectory:
#     - {run}_metadata.csv    (sample_id + group column for this run's samples)
#     - {run}_params.yml      (draft ProSIFT parameters pointing to master abundance)
#
#   All runs share the master abundance file. Module 01's cross-validation step
#   subsets the abundance matrix to the samples listed in the metadata file, so
#   per-run abundance copies are unnecessary.
#
# inputs:
#   - run_config.yml (references master abundance CSV)
#
# outputs:
#   - Per-run subdirectories with {run}_metadata.csv and {run}_params.yml
#   - run_config.yml copied to outdir for provenance
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

def write_samplesheet(filepath, samplesheet_rows):
    '''
    Write samplesheet.csv for Nextflow input.

    Each row has: run_id, abundance, metadata, params -- all absolute paths.
    Absolute paths are used so Nextflow can be run from any working directory.
    Regenerate this file after rsyncing to a new machine (absolute paths
    will differ on the cluster).
    '''
    with open(filepath, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=['run_id', 'abundance', 'metadata', 'params'])
        writer.writeheader()
        writer.writerows(samplesheet_rows)


def write_run_metadata(filepath, sample_group_map, group_column):
    '''Write the per-run metadata CSV.'''
    with open(filepath, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['sample_id', group_column])
        for sample_id, group in sample_group_map.items():
            writer.writerow([sample_id, group])


def write_run_params(filepath, run_name, group_column, abund_prefix, pep_prefix, master_abund_relpath):
    '''Write a draft params.yml for one run.'''
    today = date.today().isoformat()
    pep_line = f'  peptide_count_prefix: "{pep_prefix}"' if pep_prefix else f'  peptide_count_prefix: null'

    content = f'''# ProSIFT params -- {run_name}
# Generated by prepare_prosift_input.py on {today}

project:
  name: "{run_name}"
  organism: "mouse"

input:
  abundance_matrix: "{master_abund_relpath}"  # shared master file
  metadata: "{run_name}_metadata.csv"
  format: "csv"
  protein_id_column: "protein_id"
  abundance_type: "raw"
  abundance_prefix: "{abund_prefix}"
{pep_line}

design:
  group_column: "{group_column}"
  covariates: []
  batch_column: null

qc:
  min_samples_per_group: 2
  min_detections_per_group: 2

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
    # 2. Validate master abundance file headers
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

    samplesheet_rows = []

    for run_name, run_def in runs.items():
        sample_group_map = OrderedDict(run_def['samples'])
        group_column = run_def['group_column']

        # create run directory
        run_dir = os.path.join(outdir, run_name)
        os.makedirs(run_dir, exist_ok=True)

        # relative path from run directory to master abundance file
        master_abund_relpath = os.path.relpath(
            os.path.abspath(abund_path), os.path.abspath(run_dir)
        )

        meta_path   = os.path.join(run_dir, f'{run_name}_metadata.csv')
        params_path = os.path.join(run_dir, f'{run_name}_params.yml')

        write_run_metadata(meta_path, sample_group_map, group_column)
        write_run_params(params_path, run_name, group_column, abund_prefix, pep_prefix, master_abund_relpath)

        # collect absolute paths for samplesheet
        samplesheet_rows.append({
            'run_id':   run_name,
            'abundance': os.path.abspath(abund_path),
            'metadata':  os.path.abspath(meta_path),
            'params':    os.path.abspath(params_path),
        })

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

    # ----------------------------------------------------------
    # 5. Write samplesheet.csv for Nextflow
    # ----------------------------------------------------------
    samplesheet_path = os.path.join(outdir, 'samplesheet.csv')
    write_samplesheet(samplesheet_path, samplesheet_rows)
    print(f'Written samplesheet.csv ({len(samplesheet_rows)} runs) to {outdir}/')
    print(f'  Note: samplesheet uses absolute paths -- regenerate after rsyncing to a new machine.')

    print(f'\nDone. {len(runs)} runs written to {outdir}/')


if __name__ == '__main__':
    main()

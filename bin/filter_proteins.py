#!/usr/bin/env python3
# title: filter_proteins.py
# project: ProSIFT
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-03-26
# last modified: 2026-04-23
#
# purpose:
#   Module 01, Process 4.4. Applies the per-run detection filter to the
#   validated abundance matrix. Classifies each protein as PASSED, PARTIAL,
#   SINGLE-GROUP, SPARSE, or ABSENT based on per-group detection counts,
#   removes proteins with insufficient signal, and produces a missingness
#   overview and detection filter summary for the validation report.
#
# inputs:
#   --matrix    : {run_id}.validated_matrix.parquet (from validate_inputs.py)
#   --metadata  : {run_id}.validated_metadata.parquet (from validate_inputs.py)
#   --params    : per-run params.yml
#   --run-id    : string identifier used to name output files
#   --outdir    : output directory (default: current working directory)
#
# outputs:
#   {outdir}/{run_id}.filtered_matrix.parquet
#   {outdir}/{run_id}.detection_filter_table.csv
#   {outdir}/{run_id}.validation_report_part2.txt
#
# usage example:
#   python bin/filter_proteins.py \
#       --matrix   results/CTXcyto_WT_vs_CTXcyto_KO/CTXcyto_WT_vs_CTXcyto_KO.validated_matrix.parquet \
#       --metadata results/CTXcyto_WT_vs_CTXcyto_KO/CTXcyto_WT_vs_CTXcyto_KO.validated_metadata.parquet \
#       --params   prosift_inputs/CTXcyto_WT_vs_CTXcyto_KO/CTXcyto_WT_vs_CTXcyto_KO_params.yml \
#       --run-id   CTXcyto_WT_vs_CTXcyto_KO \
#       --outdir   results/CTXcyto_WT_vs_CTXcyto_KO

import argparse
import os
import sys
from datetime import datetime

import pandas as pd
import yaml


# ============================================================
# Argument parsing
# ============================================================

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description='ProSIFT Module 01 -- per-run detection filter.'
    )
    parser.add_argument('--matrix', required=True,
                        help='Path to validated matrix Parquet')
    parser.add_argument('--metadata', required=True,
                        help='Path to validated metadata Parquet')
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
    with open(params_path, 'r', encoding='utf-8') as f:
        return yaml.safe_load(f)


# ============================================================
# Report helper
# ============================================================

class Report:
    '''Accumulates validation report lines; writes to file at the end.

    If `report_path` is supplied, `error_exit` will flush the accumulated
    lines to that path before raising SystemExit so the user gets a partial
    report on failure.
    '''

    def __init__(self, run_id: str, report_path: str | None = None):
        self.run_id = run_id
        self._report_path = report_path
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
        '''Log the error to the report, flush it to disk, then exit non-zero.'''
        self._lines.append(f'  ERROR: {text}')
        self._lines.append('')
        self._lines.append('Validation failed. Pipeline stopped.')
        # Flush partial report so the user has diagnostic context on failure.
        if self._report_path is not None:
            try:
                parent = os.path.dirname(self._report_path)
                if parent:
                    os.makedirs(parent, exist_ok=True)
                self.write(self._report_path)
            except OSError as exc:
                print(
                    f'WARNING: could not write partial report to '
                    f'{self._report_path}: {exc}',
                    file=sys.stderr,
                )
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
# Detection counting and classification
# ============================================================

def get_abundance_cols(matrix: pd.DataFrame, id_col: str,
                       pep_prefix: str | None) -> list[str]:
    '''Return the abundance column names from the matrix.'''
    return [c for c in matrix.columns
            if c != id_col
            and (not pep_prefix or not c.startswith(pep_prefix))]


def count_detections_per_group(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    id_col: str,
    abund_cols: list[str],
    abund_prefix: str,
    group_col: str
) -> pd.DataFrame:
    '''
    For each protein, count non-missing abundance values per group.

    Returns a DataFrame indexed by protein_id with one column per group
    containing the integer detection count for that group.
    '''
    # map sample_id -> group
    sample_to_group: dict[str, str] = dict(
        zip(metadata['sample_id'].astype(str),
            metadata[group_col].astype(str))
    )

    # map abundance column -> group
    col_to_group: dict[str, str] = {
        col: sample_to_group[col[len(abund_prefix):]]
        for col in abund_cols
    }

    groups = metadata[group_col].astype(str).unique()
    abund = matrix.set_index(id_col)[abund_cols]

    detection_counts = pd.DataFrame(index=abund.index)
    for grp in groups:
        grp_cols = [c for c, g in col_to_group.items() if g == grp]
        detection_counts[grp] = abund[grp_cols].notna().sum(axis=1)

    return detection_counts


def classify_proteins(
    detection_counts: pd.DataFrame,
    min_detections: int
) -> pd.Series:
    '''
    Classify each protein into PASSED, SINGLE-GROUP, PARTIAL, SPARSE, or ABSENT.

    Rules (applied in order):
      ABSENT       -- zero detections in ALL groups
      PASSED       -- >= min_detections in ALL groups
      PARTIAL      -- >= min_detections in at least one group AND
                      sub-threshold but non-zero detections in at least
                      one other group (e.g., 2/3 WT, 1/3 KO with min=2)
      SINGLE-GROUP -- >= min_detections in exactly one group AND
                      zero detections in all other groups
      SPARSE       -- everything else (detections exist but no group meets
                      the threshold)
    '''
    total_detections = detection_counts.sum(axis=1)
    meets_threshold = detection_counts >= min_detections
    zero_in_group = detection_counts == 0

    n_groups = detection_counts.shape[1]
    n_meeting = meets_threshold.sum(axis=1)
    n_zero = zero_in_group.sum(axis=1)

    status = pd.Series('SPARSE', index=detection_counts.index, dtype=str)

    # ABSENT: no detections anywhere
    status[total_detections == 0] = 'ABSENT'

    # PASSED: meets threshold in all groups
    status[n_meeting == n_groups] = 'PASSED'

    # SINGLE-GROUP: meets threshold in exactly one group,
    #               and all other groups have zero detections
    single_group = (n_meeting == 1) & (n_zero == n_groups - 1)
    status[single_group] = 'SINGLE-GROUP'

    # PARTIAL: meets threshold in at least one group, but at least one
    #          other group has sub-threshold non-zero detections.
    #          This covers the gap between SINGLE-GROUP (zero in other
    #          groups) and PASSED (all groups meet threshold).
    partial = (n_meeting >= 1) & (n_meeting < n_groups) & (n_zero < n_groups - n_meeting)
    # Only apply to proteins still labeled SPARSE (not already SINGLE-GROUP)
    status[(status == 'SPARSE') & partial] = 'PARTIAL'

    # SPARSE: anything remaining (already set as default)

    return status


# ============================================================
# Report blocks
# ============================================================

def format_pct(n: int, total: int) -> str:
    return f'{n / total * 100:5.1f}%'


def write_detection_filter_block(
    report: Report,
    status: pd.Series,
    min_detections: int,
    n_input: int
) -> None:
    counts = status.value_counts()
    n_passed = counts.get('PASSED', 0)
    n_single = counts.get('SINGLE-GROUP', 0)
    n_partial = counts.get('PARTIAL', 0)
    n_sparse = counts.get('SPARSE', 0)
    n_absent = counts.get('ABSENT', 0)
    n_retained = n_passed + n_single + n_partial
    n_removed = n_sparse + n_absent

    report.section('DETECTION FILTER SUMMARY')
    report.line()
    report.line(f'Filter threshold: min_detections_per_group = {min_detections}')
    report.line(f'  A protein passes if detected (non-missing) in at least '
                f'{min_detections} replicates')
    report.line(f'  in at least one group.')
    report.line()
    report.line('Category definitions:')
    report.line(f'  PASSED           - Detected in >= {min_detections} replicates '
                f'in both groups.')
    report.line(f'                     Sufficient data for reliable statistical testing.')
    report.line(f'  PARTIAL          - Detected in >= {min_detections} replicates '
                f'in one group')
    report.line(f'                     but below threshold (1 to {min_detections - 1} '
                f'detections) in the other.')
    report.line(f'                     Retained: the sub-threshold group has some data')
    report.line(f'                     but not enough for reliable group-level estimates.')
    report.line(f'                     Missing values imputed as a mix of MNAR and MAR.')
    report.line(f'  SINGLE-GROUP     - Detected in >= {min_detections} replicates '
                f'in one group')
    report.line(f'                     but completely absent (0 detections) in the other.')
    report.line(f'                     Retained as potential presence/absence candidates.')
    report.line(f'                     Flagged for careful interpretation: the absent')
    report.line(f'                     group will be entirely imputed.')
    report.line(f'  SPARSE (removed) - Detected in < {min_detections} replicates '
                f'in ALL groups.')
    report.line(f'                     Insufficient data for reliable statistical testing')
    report.line(f'                     in any group.')
    report.line(f'  ABSENT (removed) - Zero detections across all samples in this run.')
    report.line(f'                     Protein was retained in the source dataset only')
    report.line(f'                     because it was detected in conditions not included')
    report.line(f'                     in this run.')
    report.line()
    report.line('Results:')
    report.line(f'  Input proteins:    {n_input:>6,}')
    report.line(f'  Passed:            {n_passed:>6,}  ({format_pct(n_passed, n_input)})')
    report.line(f'  Partial:           {n_partial:>6,}  ({format_pct(n_partial, n_input)})  [flagged]')
    report.line(f'  Single-group:      {n_single:>6,}  ({format_pct(n_single, n_input)})  [flagged]')
    report.line(f'  Sparse (removed):  {n_sparse:>6,}  ({format_pct(n_sparse, n_input)})')
    report.line(f'  Absent (removed):  {n_absent:>6,}  ({format_pct(n_absent, n_input)})')
    report.line(f'  {"─" * 34}')
    report.line(f'  Retained:          {n_retained:>6,}  ({format_pct(n_retained, n_input)})')
    report.line(f'  Removed:           {n_removed:>6,}  ({format_pct(n_removed, n_input)})')


def write_missingness_overview_block(
    report: Report,
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    id_col: str,
    abund_cols: list[str],
    abund_prefix: str,
    group_col: str
) -> None:
    n_proteins = len(matrix)
    n_samples = len(abund_cols)
    total_values = n_proteins * n_samples
    abund = matrix.set_index(id_col)[abund_cols]
    total_missing = abund.isna().sum().sum()

    sample_to_group: dict[str, str] = dict(
        zip(metadata['sample_id'].astype(str),
            metadata[group_col].astype(str))
    )
    col_to_group: dict[str, str] = {
        col: sample_to_group[col[len(abund_prefix):]]
        for col in abund_cols
    }
    groups = metadata[group_col].astype(str).unique()

    report.section('PRE-FILTER MISSINGNESS OVERVIEW')
    report.line()
    report.line(f'Total abundance values:  {total_values:>8,}  '
                f'({n_proteins:,} proteins x {n_samples} samples)')
    report.line(f'Missing values:          {total_missing:>8,}  '
                f'({total_missing / total_values * 100:.1f}%)')
    report.line()
    report.line('Per-group detection rates:')

    for grp in sorted(groups):
        grp_cols = [c for c, g in col_to_group.items() if g == grp]
        # per-sample detection counts (axis=0: sum across proteins for each sample)
        per_sample = abund[grp_cols].notna().sum(axis=0)
        grp_mean_detected = per_sample.mean()
        report.line(
            f'  {grp}:  mean {grp_mean_detected:,.0f} detected '
            f'({grp_mean_detected / n_proteins * 100:.1f}%), '
            f'range {per_sample.min():,}-{per_sample.max():,}'
        )

    report.line()
    report.line('Per-sample detection:')
    for col in abund_cols:
        sample_id = col[len(abund_prefix):]
        n_detected = abund[col].notna().sum()
        n_missing = abund[col].isna().sum()
        report.line(
            f'  {sample_id:<20}  {n_detected:>6,} detected '
            f'({n_detected / n_proteins * 100:.1f}%),  '
            f'{n_missing:>5,} missing '
            f'({n_missing / n_proteins * 100:.1f}%)'
        )


# ============================================================
# Main
# ============================================================

def main() -> None:
    args = parse_args()
    params = load_params(args.params)
    run_id = args.run_id

    # Compute report path up front so Report can flush on error_exit.
    report_path = os.path.join(
        args.outdir, f'{run_id}.validation_report_part2.txt'
    )
    report = Report(run_id, report_path=report_path)
    report.line(f'ProSIFT Validation Report -- {run_id}')
    report.line(f'Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}')
    report.line(f'Part 2 of 2: Detection Filter (Process 4.4)')
    report.line(f'  (Part 1: Input Validation -- see validation_report_part1.txt)')

    print(f'[{run_id}] Applying detection filter...')

    # --- Load inputs ---
    matrix = pd.read_parquet(args.matrix)
    metadata = pd.read_parquet(args.metadata)

    id_col: str = params['input']['protein_id_column']
    abund_prefix: str = params['input'].get('abundance_prefix', '')
    pep_prefix: str | None = params['input'].get('peptide_count_prefix') or None
    group_col: str = params['design']['group_column']
    min_detections: int = params['qc']['min_detections_per_group']

    abund_cols = get_abundance_cols(matrix, id_col, pep_prefix)
    n_input = len(matrix)

    # --- Pre-filter missingness overview (captures state before any removal) ---
    write_missingness_overview_block(
        report, matrix, metadata, id_col, abund_cols, abund_prefix, group_col
    )

    # --- Classify proteins ---
    detection_counts = count_detections_per_group(
        matrix, metadata, id_col, abund_cols, abund_prefix, group_col
    )
    status = classify_proteins(detection_counts, min_detections)

    # --- Detection filter summary block ---
    write_detection_filter_block(report, status, min_detections, n_input)

    # --- Warnings and hard stops ---
    # NOTE: PARTIAL is retained (one group passes + non-meeting groups have
    # some non-zero detections); see Section 4.4 of the module spec.
    retain_mask = status.isin(['PASSED', 'PARTIAL', 'SINGLE-GROUP'])
    n_retained = int(retain_mask.sum())
    n_removed = n_input - n_retained
    pct_removed = n_removed / n_input * 100

    if pct_removed > 20:
        report.warn(
            f'{pct_removed:.1f}% of proteins removed by detection filter '
            f'(>{20}% threshold). This may indicate the input dataset is not '
            f'appropriate for this run configuration.'
        )

    if n_retained < 50:
        report.error_exit(
            f'Only {n_retained} proteins remain after detection filtering '
            f'(minimum required: 50). Insufficient data for meaningful analysis.'
        )

    # --- Filter matrix ---
    filtered_matrix = matrix[matrix[id_col].isin(status[retain_mask].index)].copy()

    # --- Build detection filter summary table ---
    filter_table = detection_counts.copy()
    filter_table.columns = [f'detections_{g}' for g in filter_table.columns]
    filter_table['filter_status'] = status
    filter_table = filter_table.reset_index()

    # --- Summary ---
    report.section('FILTER SUMMARY')
    report.info(f'Warnings: {report.warning_count}')
    report.info(f'Retained: {n_retained:,} proteins '
                f'({n_retained / n_input * 100:.1f}%)')
    report.info(f'Removed:  {n_removed:,} proteins '
                f'({pct_removed:.1f}%)')
    report.info('Status: PASSED -- proceed to ID mapping')

    # --- Write outputs ---
    os.makedirs(args.outdir, exist_ok=True)

    matrix_path = os.path.join(args.outdir, f'{run_id}.filtered_matrix.parquet')
    table_path = os.path.join(args.outdir, f'{run_id}.detection_filter_table.csv')

    filtered_matrix.to_parquet(matrix_path, index=False)
    filter_table.to_csv(table_path, index=False)
    report.write(report_path)

    print(f'  Filtered matrix:      {matrix_path}')
    print(f'  Detection filter table: {table_path}')
    print(f'  Validation report:    {report_path}')
    print(f'[{run_id}] Detection filter complete. '
          f'{report.warning_count} warning(s). '
          f'{n_retained:,} proteins retained, {n_removed:,} removed.')


if __name__ == '__main__':
    main()

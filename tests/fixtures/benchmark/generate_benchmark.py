#!/usr/bin/env python3
# title: generate_benchmark.py
# project: ProSIFT (PROtein Statistical Integration and Filtering Tool)
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-07-08
# last modified: 2026-07-08
#
# purpose:
#   Deterministically generate the ProSIFT minimal regression benchmark: a
#   tiny synthetic 2-group (WT vs KO) proteomics dataset in which every protein's
#   true behavior is fixed BY CONSTRUCTION. Because the ground truth is set here
#   -- not read back from any pipeline output -- the generated ground-truth table
#   is a trustworthy oracle for testing Modules 01 (validate + filter), 02
#   (pre-norm QC), and 03 (normalize + impute).
#
#   Soundness by design (see the benchmark discussion): the answer key is the
#   INPUT to the data-generating process, so no module under test can influence
#   it. We can therefore assert directional / structural facts (spiked-up
#   proteins are higher in KO, MNAR proteins are absent in a whole group, a
#   zero-variance protein has CV 0) independently of the code.
#
#   Scope note: ID mapping (network), DEqMS (R), enrichment (gene sets), and the
#   database queries (network) are intentionally OUT of scope for this fixture.
#   They need different oracle strategies and are tracked separately.
#
# inputs:
#   --seed     RNG seed (default 42; fixed seed = byte-identical output)
#   --outdir   directory to write the benchmark into (default: this file's dir)
#
# outputs:
#   <outdir>/benchmark_abundance.csv     raw-intensity matrix (Module 01 input)
#   <outdir>/benchmark_metadata.csv      sample_id, genotype (Module 01 input)
#   <outdir>/benchmark_params.yml        params.yml for the run
#   <outdir>/benchmark_ground_truth.csv  the answer key (one row per protein)
#
# usage example:
#   python tests/fixtures/benchmark/generate_benchmark.py --outdir tests/fixtures/benchmark
#
#   copy/paste: python tests/fixtures/benchmark/generate_benchmark.py

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

# ============================================================
# BENCHMARK DESIGN CONSTANTS
# ============================================================
# The whole design lives here so the answer key is auditable at a glance.

# --- Sample layout: 3 WT + 3 KO (matches the real min_samples_per_group=2
#     with a margin, and the CTXcyto benchmark's 3v3 shape). ---
SAMPLES = {
    'WT': ['WT-1', 'WT-2', 'WT-3'],
    'KO': ['KO-1', 'KO-2', 'KO-3'],
}

# --- Protein cohort sizes. Total = 10 + 10 + 24 + 6 edge = 50. ---
N_TRUE_UP = 10      # higher in KO
N_TRUE_DOWN = 10    # lower in KO
N_TRUE_NULL = 24    # no genotype effect

# --- Effect sizes and noise, on log2 scale. ---
SPIKE_LOG2FC = 1.5      # |true log2 fold change| for spiked proteins (KO - WT)
BIO_NOISE_SD = 0.35     # within-group biological + technical noise (log2 units)
BASE_LOG2_LOW = 17.0    # base abundance range (log2); 2^17 ~ 1.3e5
BASE_LOG2_HIGH = 22.0   # 2^22 ~ 4.2e6 -- realistic DIA-NN intensity span

# --- Peptide counts drawn uniformly in this inclusive range. ---
PEPTIDE_MIN, PEPTIDE_MAX = 1, 10


def _all_sample_ids() -> list[str]:
    '''WT sample ids followed by KO sample ids, in fixed order.'''
    return SAMPLES['WT'] + SAMPLES['KO']


# ============================================================
# Section 1: core (non-edge) proteins
# ============================================================

def _make_regular_proteins(rng: np.random.Generator) -> tuple[list[dict], list[dict]]:
    '''
    Build the true-up / true-down / true-null proteins.

    Returns (abundance_rows, truth_rows). Each abundance row is a dict of
    raw intensities keyed by 'abundance_<sample>'; each truth row records the
    designed behavior.
    '''
    abundance_rows: list[dict] = []
    truth_rows: list[dict] = []

    cohort = (
        [('true_up', +SPIKE_LOG2FC)] * N_TRUE_UP
        + [('true_down', -SPIKE_LOG2FC)] * N_TRUE_DOWN
        + [('true_null', 0.0)] * N_TRUE_NULL
    )

    for i, (klass, log2fc) in enumerate(cohort):
        protein_id = f'BENCH{i:05d}'
        # Base (WT) mean abundance on log2 scale, unique per protein.
        base_log2 = rng.uniform(BASE_LOG2_LOW, BASE_LOG2_HIGH)

        row = {'protein_id': protein_id}
        for sid in SAMPLES['WT']:
            val_log2 = base_log2 + rng.normal(0, BIO_NOISE_SD)
            row[f'abundance_{sid}'] = 2.0 ** val_log2
        for sid in SAMPLES['KO']:
            # KO mean is shifted by the designed log2 fold change.
            val_log2 = base_log2 + log2fc + rng.normal(0, BIO_NOISE_SD)
            row[f'abundance_{sid}'] = 2.0 ** val_log2

        abundance_rows.append(row)
        truth_rows.append({
            'protein_id': protein_id,
            'class': klass,
            'true_log2fc': log2fc,
            'missingness_class': 'complete',
            # Complete data -> 3 detections in each group -> PASSED for any
            # min_detections <= 3. Set here by design, not read from any module.
            'expected_filter_status': 'PASSED',
            'notes': '',
        })

    return abundance_rows, truth_rows


# ============================================================
# Section 2: edge-case proteins
# ============================================================
# Six proteins that deliberately trip specific code paths. IDs are descriptive
# (not BENCH#####) so a failure names the scenario.

def _make_edge_proteins(rng: np.random.Generator) -> tuple[list[dict], list[dict]]:
    abundance_rows: list[dict] = []
    truth_rows: list[dict] = []
    all_ids = _all_sample_ids()

    def blank_row(pid: str) -> dict:
        return {'protein_id': pid, **{f'abundance_{s}': np.nan for s in all_ids}}

    # --- EDGE_CONST: identical value in every sample -> zero within-group
    #     variance PRE-normalization. Exercises the pre-norm zero-variance
    #     guards in prenorm_qc (KDE / Q-Q / correlation). NOTE: median
    #     normalization applies a different per-sample shift to every sample,
    #     so this protein is NO LONGER constant after normalization -- do not
    #     assert CV=0 on the normalized matrix. Complete data -> PASSED. ---
    r = blank_row('EDGE_CONST')
    const_val = 2.0 ** 20.0
    for s in all_ids:
        r[f'abundance_{s}'] = const_val
    abundance_rows.append(r)
    truth_rows.append({'protein_id': 'EDGE_CONST', 'class': 'edge',
                       'true_log2fc': 0.0, 'missingness_class': 'complete',
                       'expected_filter_status': 'PASSED',
                       'notes': 'zero within-group variance PRE-norm only; '
                                'exercises prenorm zero-variance guards'})

    # --- EDGE_MNAR_WT: present in KO, entirely missing in WT. 0 detections in
    #     WT -> filter_proteins should flag/drop; single-observation group in
    #     CV; MNAR imputation target. ---
    r = blank_row('EDGE_MNAR_WT')
    for s in SAMPLES['KO']:
        r[f'abundance_{s}'] = 2.0 ** rng.uniform(18, 21)
    abundance_rows.append(r)
    # 0 detections in WT, 3 in KO -> meets threshold in exactly one group,
    # zero in the other -> SINGLE-GROUP.
    truth_rows.append({'protein_id': 'EDGE_MNAR_WT', 'class': 'edge',
                       'true_log2fc': np.nan, 'missingness_class': 'mnar_absent_in_WT',
                       'expected_filter_status': 'SINGLE-GROUP',
                       'notes': 'missing entirely in WT (0 detections); MNAR'})

    # --- EDGE_MNAR_KO: mirror image -- absent in KO. ---
    r = blank_row('EDGE_MNAR_KO')
    for s in SAMPLES['WT']:
        r[f'abundance_{s}'] = 2.0 ** rng.uniform(18, 21)
    abundance_rows.append(r)
    truth_rows.append({'protein_id': 'EDGE_MNAR_KO', 'class': 'edge',
                       'true_log2fc': np.nan, 'missingness_class': 'mnar_absent_in_KO',
                       'expected_filter_status': 'SINGLE-GROUP',
                       'notes': 'missing entirely in KO (0 detections); MNAR'})

    # --- EDGE_PARTIAL: meets threshold in KO (3 obs) but only 1 non-zero
    #     sub-threshold detection in WT. This is the ONLY protein that lands in
    #     the PARTIAL class -- the path the 2026-04 filter audit fixed (PARTIAL
    #     proteins must be RETAINED, not silently dropped). Without it the
    #     benchmark cannot detect a regression of that fix. ---
    r = blank_row('EDGE_PARTIAL')
    base = rng.uniform(18, 21)
    for s in SAMPLES['KO']:
        r[f'abundance_{s}'] = 2.0 ** (base + rng.normal(0, BIO_NOISE_SD))
    r['abundance_WT-1'] = 2.0 ** (base + rng.normal(0, BIO_NOISE_SD))  # lone WT obs
    abundance_rows.append(r)
    truth_rows.append({'protein_id': 'EDGE_PARTIAL', 'class': 'edge',
                       'true_log2fc': np.nan, 'missingness_class': 'partial',
                       'expected_filter_status': 'PARTIAL',
                       'notes': 'KO=3 obs, WT=1 obs (non-zero, sub-threshold); '
                                'exercises PARTIAL retention path'})

    # --- EDGE_MAR: complete-ish but with two random missing values (one per
    #     group), so each group keeps >=2 detections. Exercises MAR/KNN path
    #     without dropping the protein at filtering. ---
    r = blank_row('EDGE_MAR')
    base = rng.uniform(18, 21)
    for s in all_ids:
        r[f'abundance_{s}'] = 2.0 ** (base + rng.normal(0, BIO_NOISE_SD))
    r['abundance_WT-2'] = np.nan   # one MAR dropout in WT
    r['abundance_KO-3'] = np.nan   # one MAR dropout in KO
    abundance_rows.append(r)
    # 2 detections per group == min_detections -> PASSED. This is the exact
    # lower-boundary case for the threshold (2 passes, 1 would not).
    truth_rows.append({'protein_id': 'EDGE_MAR', 'class': 'edge',
                       'true_log2fc': 0.0, 'missingness_class': 'mar_random',
                       'expected_filter_status': 'PASSED',
                       'notes': '1 random dropout per group; 2 obs per group '
                                '== min_detections boundary (PASSED)'})

    # --- EDGE_ZERO: contains a literal 0 intensity. apply_log2 must convert it
    #     to NaN and record a warning (log2(0) = -inf otherwise). ---
    r = blank_row('EDGE_ZERO')
    base = rng.uniform(18, 21)
    for s in all_ids:
        r[f'abundance_{s}'] = 2.0 ** (base + rng.normal(0, BIO_NOISE_SD))
    r['abundance_KO-1'] = 0.0
    abundance_rows.append(r)
    # Literal 0 intensity. Under the 2026-07-08 non-positive policy, validation
    # (Module 01) converts raw <= 0 to NaN and zeros the paired peptide count, so
    # KO-1 becomes missing: KO keeps 2 detections -> PASSED. Exercises the
    # validation conversion + peptide-count zeroing end to end.
    truth_rows.append({'protein_id': 'EDGE_ZERO', 'class': 'edge',
                       'true_log2fc': 0.0, 'missingness_class': 'contains_zero',
                       'expected_filter_status': 'PASSED',
                       'notes': 'literal 0 in KO-1; validation converts to NaN '
                                '(peptide zeroed); KO keeps 2 obs -> PASSED'})

    # --- EDGE_NEG: a negative raw value (e.g. background-subtracted). Validation
    #     converts <= 0 to NaN, so WT-2 becomes missing: WT keeps 2 detections,
    #     KO keeps 3 -> PASSED. Exercises the <= 0 path, not just == 0. ---
    r = blank_row('EDGE_NEG')
    base = rng.uniform(18, 21)
    for s in all_ids:
        r[f'abundance_{s}'] = 2.0 ** (base + rng.normal(0, BIO_NOISE_SD))
    r['abundance_WT-2'] = -5.0
    abundance_rows.append(r)
    truth_rows.append({'protein_id': 'EDGE_NEG', 'class': 'edge',
                       'true_log2fc': 0.0, 'missingness_class': 'contains_negative',
                       'expected_filter_status': 'PASSED',
                       'notes': 'negative in WT-2; validation converts <=0 to NaN '
                                '-> WT=2 obs -> PASSED; exercises the <=0 path'})

    # --- EDGE_LOWDET: 1 observed value per group. Both groups fall BELOW
    #     min_detections (1 < 2), so no group meets threshold but detections
    #     exist -> SPARSE (dropped). This is a below-threshold case, NOT the
    #     boundary (the boundary == 2 is covered by EDGE_MAR). ---
    r = blank_row('EDGE_LOWDET')
    r['abundance_WT-1'] = 2.0 ** rng.uniform(17, 18)
    r['abundance_KO-1'] = 2.0 ** rng.uniform(17, 18)
    abundance_rows.append(r)
    truth_rows.append({'protein_id': 'EDGE_LOWDET', 'class': 'edge',
                       'true_log2fc': 0.0, 'missingness_class': 'low_detection',
                       'expected_filter_status': 'SPARSE',
                       'notes': '1 obs per group (< min_detections=2); '
                                'below threshold -> SPARSE (dropped)'})

    return abundance_rows, truth_rows


# ============================================================
# Section 3: peptide counts
# ============================================================

def _add_peptide_counts(abundance_df: pd.DataFrame, rng: np.random.Generator) -> pd.DataFrame:
    '''
    Append peptide_count_<sample> columns. A count is positive only where an
    abundance value is quantified (> 0); a missing OR non-positive abundance
    gets 0 peptides. This matches the real convention (peptide count 0 <-> no
    detection) that DEqMS relies on, and it means the committed benchmark is
    already consistent with the Module 01 non-positive conversion (a raw 0 or
    negative is not-detected, so it carries no supporting peptides).
    '''
    out = abundance_df.copy()
    for sid in _all_sample_ids():
        quantified = out[f'abundance_{sid}'] > 0   # NaN > 0 and (<=0) both False
        counts = rng.integers(PEPTIDE_MIN, PEPTIDE_MAX + 1, size=len(out))
        out[f'peptide_count_{sid}'] = np.where(quantified, counts, 0).astype(int)
    return out


# ============================================================
# Section 4: params.yml
# ============================================================

def _params_yaml_text() -> str:
    '''
    Minimal params.yml covering the fields Modules 01-03 read. Mirrors the real
    prosift_inputs/*/params.yml structure so the benchmark exercises the same
    code paths. Written as literal text (not yaml.dump) to keep comments.
    '''
    return '''# ProSIFT params -- minimal regression benchmark
# Generated by tests/fixtures/benchmark/generate_benchmark.py

project:
  name: "benchmark_WT_vs_KO"
  organism: "mouse"

input:
  abundance_matrix: "benchmark_abundance.csv"
  metadata: "benchmark_metadata.csv"
  format: "csv"
  protein_id_column: "protein_id"
  abundance_type: "raw"
  abundance_prefix: "abundance_"
  peptide_count_prefix: "peptide_count_"

design:
  group_column: "genotype"
  covariates: []
  batch_column: null
  contrasts:
    - "KO_vs_WT"

qc:
  min_samples_per_group: 2
  min_detections_per_group: 2

normalization:
  method: "median"

imputation:
  mode: "mixed"
  mnar_method: "minprob"
  mar_method: "knn"
  minprob_quantile: 0.01
  minprob_scale: 0.3
  knn_k: 3
  random_seed: 42
'''


# ============================================================
# Section 5: assembly + write
# ============================================================

def build_benchmark(seed: int) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    '''
    Assemble the full benchmark. A single seeded generator is threaded through
    every step so the whole dataset is reproducible from `seed` alone.
    '''
    rng = np.random.default_rng(seed)

    reg_abund, reg_truth = _make_regular_proteins(rng)
    edge_abund, edge_truth = _make_edge_proteins(rng)

    abundance_df = pd.DataFrame(reg_abund + edge_abund)
    truth_df = pd.DataFrame(reg_truth + edge_truth)

    # Column order: protein_id, all abundance cols, then peptide cols.
    abund_cols = [f'abundance_{s}' for s in _all_sample_ids()]
    abundance_df = abundance_df[['protein_id', *abund_cols]]
    abundance_df = _add_peptide_counts(abundance_df, rng)

    metadata_df = pd.DataFrame({
        'sample_id': _all_sample_ids(),
        'genotype': (['WT'] * len(SAMPLES['WT'])) + (['KO'] * len(SAMPLES['KO'])),
    })

    return abundance_df, metadata_df, truth_df


def main() -> None:
    parser = argparse.ArgumentParser(
        prog='generate_benchmark.py',
        description='Generate the ProSIFT minimal regression benchmark (deterministic).',
    )
    parser.add_argument('--seed', type=int, default=42,
                        help='RNG seed; fixed seed gives byte-identical output (default 42)')
    parser.add_argument('--outdir', type=Path, default=Path(__file__).resolve().parent,
                        help='Output directory (default: this script\'s directory)')
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    abundance_df, metadata_df, truth_df = build_benchmark(args.seed)

    # Round abundances to whole numbers -- real DIA-NN intensities are large,
    # and rounding keeps the committed CSV small and diff-friendly.
    abund_cols = [c for c in abundance_df.columns if c.startswith('abundance_')]
    abundance_df[abund_cols] = abundance_df[abund_cols].round(1)

    abundance_df.to_csv(args.outdir / 'benchmark_abundance.csv', index=False)
    metadata_df.to_csv(args.outdir / 'benchmark_metadata.csv', index=False)
    truth_df.to_csv(args.outdir / 'benchmark_ground_truth.csv', index=False)
    (args.outdir / 'benchmark_params.yml').write_text(_params_yaml_text(), encoding='utf-8')

    # Console summary for the operator's visual inspection.
    n_prot = len(abundance_df)
    n_missing = int(abundance_df[abund_cols].isna().sum().sum())
    print(f'Benchmark written to {args.outdir}')
    print(f'  proteins:        {n_prot}')
    print(f'  samples:         {len(metadata_df)} ({metadata_df.genotype.value_counts().to_dict()})')
    print(f'  missing cells:   {n_missing} / {n_prot * len(abund_cols)}')
    print(f'  class breakdown: {truth_df["class"].value_counts().to_dict()}')
    print('  files: benchmark_abundance.csv, benchmark_metadata.csv, '
          'benchmark_params.yml, benchmark_ground_truth.csv')


if __name__ == '__main__':
    main()

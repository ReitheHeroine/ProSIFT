#!/usr/bin/env python3
# title: test_benchmark_fixture.py
# project: ProSIFT (PROtein Statistical Integration and Filtering Tool)
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-07-08
# last modified: 2026-07-08
#
# purpose:
#   Validate the minimal regression benchmark BEFORE any module is tested
#   against it. A regression benchmark is only as trustworthy as its answer
#   key, so this file is the oracle's own guard. It checks two things:
#
#     1. DETERMINISM -- regenerating from the fixed seed reproduces the
#        committed files byte-for-byte. If this fails, the benchmark is not a
#        stable reference and every downstream regression test is meaningless.
#
#     2. INTERNAL CONSISTENCY -- the committed abundance data actually exhibits
#        the behavior the ground-truth table claims (spiked-up proteins really
#        are higher in KO, MNAR proteins really are absent in a whole group,
#        the zero-variance protein really is constant, etc.). This is the
#        "how do we know the benchmark is sound" check, made executable:
#        it compares two INDEPENDENT products of the generator (the data and
#        the answer key) and confirms they agree.
#
#   Note the direction of proof: these tests do NOT run any pipeline module.
#   They validate the fixture in isolation, so a failure here is unambiguous --
#   it means the benchmark generator is wrong, never that a module is wrong.
#
# inputs:
#   tests/fixtures/benchmark/benchmark_abundance.csv
#   tests/fixtures/benchmark/benchmark_metadata.csv
#   tests/fixtures/benchmark/benchmark_ground_truth.csv
#   tests/fixtures/benchmark/generate_benchmark.py (imported for determinism check)
#
# outputs:
#   Test results (stdout via pytest)
#
# usage example:
#   pytest tests/test_benchmark_fixture.py -v
#
#   copy/paste: pytest tests/test_benchmark_fixture.py -v

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# Locate the benchmark fixture directory and make its generator importable.
BENCH_DIR = Path(__file__).resolve().parent / 'fixtures' / 'benchmark'
sys.path.insert(0, str(BENCH_DIR))

import generate_benchmark as gen  # noqa: E402

WT = ['abundance_WT-1', 'abundance_WT-2', 'abundance_WT-3']
KO = ['abundance_KO-1', 'abundance_KO-2', 'abundance_KO-3']


# ============================================================
# FIXTURES: load the committed benchmark once per module
# ============================================================

@pytest.fixture(scope='module')
def abundance() -> pd.DataFrame:
    return pd.read_csv(BENCH_DIR / 'benchmark_abundance.csv').set_index('protein_id')


@pytest.fixture(scope='module')
def truth() -> pd.DataFrame:
    return pd.read_csv(BENCH_DIR / 'benchmark_ground_truth.csv').set_index('protein_id')


# ============================================================
# Section 1: determinism guard
# ============================================================

class TestDeterminism:

    def test_regeneration_is_byte_identical(self):
        '''
        Regenerate the benchmark from the default seed into a temp dir and
        compare against the committed files. Byte-identical output is what
        lets a frozen reference detect real drift rather than RNG noise.
        '''
        abund, _meta, _truth = gen.build_benchmark(seed=42)

        # Apply the same rounding main() applies before writing, then compare
        # the CSV serialization to the committed file.
        abund_cols = [c for c in abund.columns if c.startswith('abundance_')]
        abund[abund_cols] = abund[abund_cols].round(1)

        regenerated = abund.to_csv(index=False)
        committed = (BENCH_DIR / 'benchmark_abundance.csv').read_text()
        assert regenerated == committed, (
            'Regenerated abundance CSV differs from the committed file. '
            'The generator is non-deterministic or has changed.'
        )

    def test_different_seed_changes_data(self):
        '''Sanity check on the seed: a different seed must change the data,
        otherwise the "seeded" claim is vacuous.'''
        a42, _, _ = gen.build_benchmark(seed=42)
        a7, _, _ = gen.build_benchmark(seed=7)
        # Same shape, but the random draws differ.
        assert not np.allclose(
            a42[WT].values, a7[WT].values, equal_nan=True
        )


# ============================================================
# Section 2: structural consistency (shape + bookkeeping)
# ============================================================

class TestStructure:

    def test_ground_truth_covers_every_protein(self, abundance, truth):
        '''Every protein in the matrix has exactly one answer-key row.'''
        assert set(abundance.index) == set(truth.index)
        assert not truth.index.duplicated().any()

    def test_expected_cohort_sizes(self, truth):
        counts = truth['class'].value_counts().to_dict()
        assert counts['true_up'] == gen.N_TRUE_UP
        assert counts['true_down'] == gen.N_TRUE_DOWN
        assert counts['true_null'] == gen.N_TRUE_NULL
        # CONST, MNAR_WT, MNAR_KO, PARTIAL, MAR, ZERO, NEG, LOWDET
        assert counts['edge'] == 8

    def test_expected_filter_status_present_and_valid(self, truth):
        '''
        Every protein carries a designed expected_filter_status drawn from the
        classifier's vocabulary, and the four non-trivial classes the Module 01
        filter can emit are all represented (so the benchmark can regress each).
        '''
        valid = {'PASSED', 'PARTIAL', 'SINGLE-GROUP', 'SPARSE', 'ABSENT'}
        assert 'expected_filter_status' in truth.columns
        assert set(truth['expected_filter_status']).issubset(valid)
        present = set(truth['expected_filter_status'])
        for required in ('PASSED', 'PARTIAL', 'SINGLE-GROUP', 'SPARSE'):
            assert required in present, f'no benchmark protein exercises {required}'

    def test_peptide_counts_zero_iff_not_quantified(self, abundance):
        '''
        Generator contract, matching the pipeline invariant peptide count
        0 <-> no detection: a count is positive exactly where abundance is
        quantified (> 0); a missing OR non-positive abundance gets 0 peptides.
        '''
        for ab_col in WT + KO:
            sid = ab_col.replace('abundance_', '')
            pep_col = f'peptide_count_{sid}'
            quantified = abundance[ab_col] > 0   # NaN and <= 0 both False
            assert (abundance.loc[quantified, pep_col] > 0).all()
            assert (abundance.loc[~quantified, pep_col] == 0).all()


# ============================================================
# Section 3: the answer key matches the data (the real soundness check)
# ============================================================
# Each test confirms that an INDEPENDENT reading of the abundance matrix
# reproduces what the ground-truth table asserts. Agreement here is the
# evidence that the oracle is trustworthy.

class TestGroundTruthConsistency:

    def test_spiked_up_proteins_are_higher_in_ko(self, abundance, truth):
        '''true_up proteins: mean log2 abundance must be higher in KO than WT.'''
        up_ids = truth.index[truth['class'] == 'true_up']
        for pid in up_ids:
            wt_mean = np.log2(abundance.loc[pid, WT]).mean()
            ko_mean = np.log2(abundance.loc[pid, KO]).mean()
            assert ko_mean > wt_mean, f'{pid} labeled true_up but not higher in KO'

    def test_spiked_down_proteins_are_lower_in_ko(self, abundance, truth):
        down_ids = truth.index[truth['class'] == 'true_down']
        for pid in down_ids:
            wt_mean = np.log2(abundance.loc[pid, WT]).mean()
            ko_mean = np.log2(abundance.loc[pid, KO]).mean()
            assert ko_mean < wt_mean, f'{pid} labeled true_down but not lower in KO'

    def test_spike_direction_separates_from_null(self, abundance, truth):
        '''
        Aggregate check: the mean |KO-WT| log2 difference of spiked proteins
        should clearly exceed that of true_null proteins. Confirms the designed
        effect size is actually present and detectable above the noise.
        '''
        def mean_abs_diff(ids):
            diffs = [
                abs(np.log2(abundance.loc[p, KO]).mean()
                    - np.log2(abundance.loc[p, WT]).mean())
                for p in ids
            ]
            return float(np.mean(diffs))

        spiked = truth.index[truth['class'].isin(['true_up', 'true_down'])]
        null = truth.index[truth['class'] == 'true_null']
        assert mean_abs_diff(spiked) > 3 * mean_abs_diff(null)

    def test_edge_const_is_constant(self, abundance):
        '''EDGE_CONST must have zero variance across all samples.'''
        vals = abundance.loc['EDGE_CONST', WT + KO].astype(float)
        assert vals.std(ddof=1) == pytest.approx(0.0, abs=1e-9)

    def test_edge_mnar_wt_absent_in_wt_only(self, abundance):
        '''EDGE_MNAR_WT: fully missing in WT, fully observed in KO.'''
        assert abundance.loc['EDGE_MNAR_WT', WT].isna().all()
        assert abundance.loc['EDGE_MNAR_WT', KO].notna().all()

    def test_edge_mnar_ko_absent_in_ko_only(self, abundance):
        assert abundance.loc['EDGE_MNAR_KO', KO].isna().all()
        assert abundance.loc['EDGE_MNAR_KO', WT].notna().all()

    def test_edge_partial_detection_pattern(self, abundance, truth):
        '''
        EDGE_PARTIAL must have 3 observed in KO and exactly 1 in WT, the pattern
        that classifies as PARTIAL (meets threshold in KO, non-zero sub-threshold
        in WT). This is the only protein guarding the PARTIAL-retention fix.
        '''
        assert abundance.loc['EDGE_PARTIAL', KO].notna().sum() == 3
        assert abundance.loc['EDGE_PARTIAL', WT].notna().sum() == 1
        assert truth.loc['EDGE_PARTIAL', 'expected_filter_status'] == 'PARTIAL'

    def test_edge_mar_keeps_two_per_group(self, abundance):
        '''EDGE_MAR: one dropout per group, so each group retains >=2 obs.'''
        assert abundance.loc['EDGE_MAR', WT].notna().sum() == 2
        assert abundance.loc['EDGE_MAR', KO].notna().sum() == 2

    def test_edge_zero_contains_a_literal_zero(self, abundance):
        '''EDGE_ZERO must carry exactly one 0.0 intensity for the conversion path.'''
        vals = abundance.loc['EDGE_ZERO', WT + KO].astype(float)
        assert (vals == 0.0).sum() == 1

    def test_edge_neg_contains_a_negative(self, abundance):
        '''EDGE_NEG must carry exactly one negative intensity for the <= 0 path.'''
        vals = abundance.loc['EDGE_NEG', WT + KO].astype(float)
        assert (vals < 0).sum() == 1

    def test_edge_lowdet_one_observation_per_group(self, abundance):
        assert abundance.loc['EDGE_LOWDET', WT].notna().sum() == 1
        assert abundance.loc['EDGE_LOWDET', KO].notna().sum() == 1

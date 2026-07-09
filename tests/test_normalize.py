#!/usr/bin/env python3
# title: test_normalize.py
# project: ProSIFT (PROtein Statistical Integration and Filtering Tool)
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-07-08
# last modified: 2026-07-08
#
# purpose:
#   Unit tests for Module 03 NORMALIZE (bin/normalize.py). Scope of this file:
#   the pure numerical transforms and the method dispatcher -- the parts that
#   decide what numbers land in normalized_matrix.parquet and cv_summary.parquet.
#
#   Design philosophy demonstrated here (worked example):
#     1. GROUND-TRUTH tests -- tiny hand-computed inputs where the correct
#        output is known independently of the code under test. If the code and
#        the test agree, and the test's expected value was derived by hand,
#        that is real evidence of correctness (not just "the code did what the
#        code does").
#     2. PROPERTY / INVARIANT tests -- assert a mathematical guarantee that
#        must hold for ANY valid input (e.g. "after median normalization, every
#        sample shares the same median"). These catch bugs that ground-truth
#        cases with specific numbers might miss.
#     3. CONTRACT tests -- the function rejects bad input loudly (raises) rather
#        than producing silently-wrong output. ProSIFT's handoff QC protocol
#        explicitly lists "silent data loss" as a red flag; these enforce it.
#     4. EDGE cases -- NaN handling, single-observation groups, zero-variance.
#
#   VSN (normalize_vsn) is intentionally NOT unit-tested here: it requires R +
#   Bioconductor vsn via rpy2. Those belong in a separate module marked
#   @pytest.mark.requires_r so the default `pytest` run stays fast and offline.
#
# inputs:
#   None (tests build inputs in-memory; no files read)
#
# outputs:
#   Test results (stdout via pytest)
#
# usage example:
#   pytest tests/test_normalize.py -v
#   pytest tests/test_normalize.py::TestNormalizeMedian -v
#   pytest tests/test_normalize.py -k 'quantile'
#
#   copy/paste: pytest tests/test_normalize.py -v

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# Add bin/ to path so we can import normalize directly, no package install
# needed. Mirrors the mechanism in test_prenorm_qc.py / test_prosift_cache.py.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'bin'))

from normalize import (
    apply_log2,
    build_group_map,
    compute_cv_summary,
    normalize_median,
    normalize_quantile,
    prepare_abundance,
    run_normalization,
)

# ============================================================
# TEST HELPERS
# ============================================================
# Small in-memory builders. Each returns the minimal valid input for the
# function under test, so a test failure points at the function, not at
# unrelated setup.

def _log2_frame(data: dict[str, list[float]], protein_ids=None) -> pd.DataFrame:
    '''
    Build a log2-scale abundance frame: rows = proteins, columns = samples.
    `data` maps sample_id -> list of per-protein values.
    '''
    df = pd.DataFrame(data)
    if protein_ids is None:
        protein_ids = [f'P{i:05d}' for i in range(len(df))]
    df.index = pd.Index(protein_ids, name='protein_id')
    return df


def _params(method: str, abundance_type: str) -> dict:
    '''Minimal params dict for run_normalization: only the keys it reads.'''
    return {
        'input': {'abundance_type': abundance_type},
        'normalization': {'method': method},
    }


# ============================================================
# Section 1: normalize_median
# ============================================================
# Median normalization shifts each sample so all samples share the global
# median (median of the per-sample medians). NaN values pass through untouched.

class TestNormalizeMedian:

    def test_ground_truth_shifts(self):
        '''
        Hand-computed case. Per-sample medians: A=25, B=2.5, C=250.
        Global median = median(25, 2.5, 250) = 25.
        Shifts: A -> 0, B -> +22.5, C -> -225.
        '''
        df = _log2_frame({
            'A': [10.0, 20.0, 30.0, 40.0],   # median 25
            'B': [1.0, 2.0, 3.0, 4.0],       # median 2.5
            'C': [100.0, 200.0, 300.0, 400.0],  # median 250
        })
        out = normalize_median(df)

        # A was already at the global median -> unchanged
        assert out['A'].tolist() == [10.0, 20.0, 30.0, 40.0]
        # B shifted up by 22.5
        assert out['B'].tolist() == [23.5, 24.5, 25.5, 26.5]
        # C shifted down by 225: [100,200,300,400] - 225 -> [-125,-25,75,175]
        # (post-shift median = (-25 + 75) / 2 = 25 = global median, as required)
        assert out['C'].tolist() == [-125.0, -25.0, 75.0, 175.0]

    def test_invariant_medians_all_equal(self):
        '''
        Property: after median normalization every sample's median equals the
        global median, regardless of the input distribution.
        '''
        rng = np.random.default_rng(seed=7)
        df = _log2_frame({
            f'S{j}': rng.normal(loc=10 * j, scale=3, size=50).tolist()
            for j in range(1, 6)
        })
        out = normalize_median(df)

        sample_medians = out.median()
        # All post-norm medians should collapse to a single value.
        assert np.allclose(sample_medians.values, sample_medians.iloc[0])

    def test_nan_preserved(self):
        '''NaN positions must remain NaN (missing != imputed here).'''
        df = _log2_frame({
            'A': [10.0, np.nan, 30.0, 40.0],
            'B': [1.0, 2.0, np.nan, 4.0],
        })
        out = normalize_median(df)
        assert np.isnan(out.loc['P00001', 'A'])  # row index 1
        assert np.isnan(out.loc['P00002', 'B'])  # row index 2

    def test_shift_is_purely_additive(self):
        '''
        A median shift preserves within-sample differences: the gap between any
        two proteins in a sample is unchanged. This distinguishes a correct
        additive shift from an accidental scaling.
        '''
        df = _log2_frame({'A': [5.0, 8.0, 20.0], 'B': [1.0, 2.0, 3.0]})
        out = normalize_median(df)
        for col in df.columns:
            orig_gaps = df[col].diff().dropna().values
            new_gaps = out[col].diff().dropna().values
            assert np.allclose(orig_gaps, new_gaps)


# ============================================================
# Section 2: normalize_quantile
# ============================================================
# Quantile normalization forces every sample to share an identical value
# distribution (the mean of the sorted columns). After normalization, the
# sorted values of every sample are identical.

class TestNormalizeQuantile:

    def test_ground_truth_two_samples(self):
        '''
        Two samples, no ties, no NaN. Reference = row-wise mean of sorted cols.

          A sorted: [1, 3, 5, 7]
          B sorted: [2, 4, 6, 8]
          reference (mean at each rank): [1.5, 3.5, 5.5, 7.5]

        A's original order is [5, 7, 1, 3] -> ranks [3, 4, 1, 2]
          -> reference[rank-1] = [5.5, 7.5, 1.5, 3.5]
        B's original order is [8, 6, 4, 2] -> ranks [4, 3, 2, 1]
          -> reference[rank-1] = [7.5, 5.5, 3.5, 1.5]
        '''
        df = _log2_frame({
            'A': [5.0, 7.0, 1.0, 3.0],
            'B': [8.0, 6.0, 4.0, 2.0],
        })
        out = normalize_quantile(df)
        assert out['A'].tolist() == [5.5, 7.5, 1.5, 3.5]
        assert out['B'].tolist() == [7.5, 5.5, 3.5, 1.5]

    def test_invariant_identical_sorted_distributions(self):
        '''
        Property: after quantile normalization every sample has the same sorted
        value vector. This is the defining guarantee of the method.
        '''
        rng = np.random.default_rng(seed=11)
        df = _log2_frame({
            f'S{j}': rng.normal(loc=20, scale=4, size=30).tolist()
            for j in range(6)
        })
        out = normalize_quantile(df)

        sorted_first = np.sort(out.iloc[:, 0].values)
        for col in out.columns:
            assert np.allclose(np.sort(out[col].values), sorted_first)

    def test_rank_order_preserved(self):
        '''
        Quantile normalization changes values but must preserve within-sample
        rank order: the largest input stays the largest output.
        '''
        df = _log2_frame({
            'A': [5.0, 7.0, 1.0, 3.0],
            'B': [8.0, 6.0, 4.0, 2.0],
        })
        out = normalize_quantile(df)
        for col in df.columns:
            assert (df[col].rank().values == out[col].rank().values).all()

    def test_nan_preserved(self):
        '''NaN stays NaN; observed values are still mapped to the reference.'''
        df = _log2_frame({
            'A': [5.0, np.nan, 1.0, 3.0],
            'B': [8.0, 6.0, 4.0, 2.0],
        })
        out = normalize_quantile(df)
        assert np.isnan(out.loc['P00001', 'A'])
        # The observed cells remain finite.
        assert out['A'].notna().sum() == 3


# ============================================================
# Section 3: apply_log2
# ============================================================

class TestApplyLog2:

    def test_ground_truth_powers_of_two(self):
        '''log2 of exact powers of two gives exact integers.'''
        df = _log2_frame({'A': [1.0, 2.0, 4.0, 8.0]})
        out = apply_log2(df, warnings=[])
        assert out['A'].tolist() == [0.0, 1.0, 2.0, 3.0]

    def test_zeros_converted_to_nan_and_warned(self):
        '''
        Zeros must become NaN (log2(0) = -inf would poison downstream stats)
        and the event must be recorded in the warnings list for the summary.
        This is a defensive backstop -- validation should remove them first.
        '''
        df = _log2_frame({'A': [1.0, 0.0, 4.0]})
        warnings: list[str] = []
        out = apply_log2(df, warnings=warnings)
        assert np.isnan(out.loc['P00001', 'A'])   # the former zero
        assert len(warnings) == 1
        assert 'non-positive' in warnings[0].lower()

    def test_negatives_converted_to_nan_and_warned(self):
        '''
        Negative values must also become NaN (log2 of a negative is NaN with a
        RuntimeWarning, which the project's pytest config escalates to an error).
        The guard covers <= 0, not just == 0.
        '''
        df = _log2_frame({'A': [1.0, -3.0, 4.0]})
        warnings: list[str] = []
        out = apply_log2(df, warnings=warnings)
        assert np.isnan(out.loc['P00001', 'A'])   # the former negative
        assert len(warnings) == 1

    def test_no_nonpositive_no_warning(self):
        df = _log2_frame({'A': [1.0, 2.0, 4.0]})
        warnings: list[str] = []
        apply_log2(df, warnings=warnings)
        assert warnings == []


# ============================================================
# Section 4: compute_cv_summary
# ============================================================
# CV is computed on the back-transformed LINEAR scale (2**value), per protein
# per group, using observed values only. Groups with <2 observations -> NaN.

class TestComputeCvSummary:

    def test_ground_truth_cv(self):
        '''
        One group 'WT' with three samples. Protein P00000 log2 = [0, 1, 2]
        -> linear [1, 2, 4], mean = 7/3, sample SD (ddof=1) = sqrt(2.3333...).
        CV = SD / mean = 1.52753.../2.33333... = 0.65465...
        Protein P00001 log2 = [1, 1, 1] -> linear [2, 2, 2] -> CV = 0.
        '''
        norm_df = _log2_frame(
            {'WT-1': [0.0, 1.0], 'WT-2': [1.0, 1.0], 'WT-3': [2.0, 1.0]},
            protein_ids=['P00000', 'P00001'],
        )
        group_map = {'WT-1': 'WT', 'WT-2': 'WT', 'WT-3': 'WT'}
        cv = compute_cv_summary(norm_df, group_map, ['WT-1', 'WT-2', 'WT-3'])

        cv = cv.set_index('protein_id')
        assert cv.loc['P00000', 'cv_WT'] == pytest.approx(0.65465367, rel=1e-5)
        assert cv.loc['P00001', 'cv_WT'] == pytest.approx(0.0, abs=1e-12)
        assert cv.loc['P00000', 'n_observed_WT'] == 3

    def test_single_observation_group_is_nan(self):
        '''Fewer than 2 observed values -> CV undefined -> NaN (not 0, not error).'''
        norm_df = _log2_frame(
            {'WT-1': [1.0], 'WT-2': [np.nan], 'WT-3': [np.nan]},
            protein_ids=['P00000'],
        )
        group_map = {'WT-1': 'WT', 'WT-2': 'WT', 'WT-3': 'WT'}
        cv = compute_cv_summary(norm_df, group_map, ['WT-1', 'WT-2', 'WT-3'])
        assert np.isnan(cv.loc[0, 'cv_WT'])
        assert cv.loc[0, 'n_observed_WT'] == 1

    def test_two_groups_reported_separately(self):
        '''CV and n_observed columns exist per group; groups are independent.'''
        norm_df = _log2_frame(
            {'WT-1': [0.0], 'WT-2': [2.0], 'KO-1': [1.0], 'KO-2': [1.0]},
            protein_ids=['P00000'],
        )
        group_map = {'WT-1': 'WT', 'WT-2': 'WT', 'KO-1': 'KO', 'KO-2': 'KO'}
        cv = compute_cv_summary(
            norm_df, group_map, ['WT-1', 'WT-2', 'KO-1', 'KO-2']
        )
        assert set(['cv_WT', 'cv_KO', 'n_observed_WT', 'n_observed_KO']).issubset(
            cv.columns
        )
        # KO is constant on linear scale -> CV 0; WT is not.
        assert cv.loc[0, 'cv_KO'] == pytest.approx(0.0, abs=1e-12)
        assert cv.loc[0, 'cv_WT'] > 0


# ============================================================
# Section 5: run_normalization dispatcher (contract tests)
# ============================================================
# The dispatcher routes on (method, abundance_type). These tests assert it
# fails LOUDLY on invalid combinations and applies log2 only when it should.

class TestRunNormalizationDispatch:

    def test_unknown_method_raises(self):
        df = _log2_frame({'A': [1.0, 2.0], 'B': [3.0, 4.0]})
        with pytest.raises(ValueError, match='Unknown normalization'):
            run_normalization(df, _params('zscore', 'log2'), warnings=[])

    def test_vsn_requires_raw_data(self):
        '''VSN on log2 data is a user error and must be rejected, not silently run.'''
        df = _log2_frame({'A': [1.0, 2.0], 'B': [3.0, 4.0]})
        with pytest.raises(ValueError, match='VSN requires raw'):
            run_normalization(df, _params('vsn', 'log2'), warnings=[])

    def test_unknown_abundance_type_raises(self):
        df = _log2_frame({'A': [1.0, 2.0], 'B': [3.0, 4.0]})
        with pytest.raises(ValueError, match='Unknown abundance_type'):
            run_normalization(df, _params('median', 'bogus'), warnings=[])

    def test_normalized_passthrough(self):
        '''abundance_type 'normalized' -> data returned unchanged, no log2.'''
        df = _log2_frame({'A': [1.0, 2.0], 'B': [3.0, 4.0]})
        out, method_used, log2_applied = run_normalization(
            df, _params('median', 'normalized'), warnings=[]
        )
        assert log2_applied is False
        assert 'pass-through' in method_used
        pd.testing.assert_frame_equal(out, df)

    def test_log2_data_not_transformed_again(self):
        '''
        abundance_type 'log2' + method 'none' -> values pass straight through
        with no second log2. Guards against double-log bugs.
        '''
        df = _log2_frame({'A': [10.0, 20.0], 'B': [30.0, 40.0]})
        out, method_used, log2_applied = run_normalization(
            df, _params('none', 'log2'), warnings=[]
        )
        assert log2_applied is False
        assert method_used == 'none'
        pd.testing.assert_frame_equal(out, df)

    def test_raw_data_gets_log2(self):
        '''abundance_type 'raw' + method 'none' -> log2 applied, flag set True.'''
        df = _log2_frame({'A': [1.0, 4.0], 'B': [2.0, 8.0]})  # raw intensities
        out, _method_used, log2_applied = run_normalization(
            df, _params('none', 'raw'), warnings=[]
        )
        assert log2_applied is True
        # log2([1,4]) = [0,2]; log2([2,8]) = [1,3]
        assert out['A'].tolist() == [0.0, 2.0]
        assert out['B'].tolist() == [1.0, 3.0]


# ============================================================
# Section 6: prepare_abundance and build_group_map (I/O parsing contract)
# ============================================================

class TestPrepareAbundance:

    def test_strips_prefix_and_separates_peptides(self):
        matrix = pd.DataFrame({
            'protein_id': ['P1', 'P2'],
            'abundance_S1': [10.0, 20.0],
            'abundance_S2': [30.0, 40.0],
            'peptide_count_S1': [3, 5],
            'peptide_count_S2': [4, 6],
        })
        params = {
            'input': {
                'abundance_prefix': 'abundance_',
                'peptide_count_prefix': 'peptide_count_',
            }
        }
        raw_df, peptide_df, sample_ids = prepare_abundance(matrix, params)
        assert sample_ids == ['S1', 'S2']          # prefix stripped
        assert list(raw_df.columns) == ['S1', 'S2']
        assert list(peptide_df.columns) == ['peptide_count_S1', 'peptide_count_S2']
        assert raw_df.loc['P1', 'S1'] == 10.0       # indexed by protein_id

    def test_no_abundance_columns_raises(self):
        '''A matrix with only peptide columns is malformed -> hard stop.'''
        matrix = pd.DataFrame({
            'protein_id': ['P1'],
            'peptide_count_S1': [3],
        })
        params = {'input': {'abundance_prefix': '', 'peptide_count_prefix': 'peptide_count_'}}
        with pytest.raises(ValueError, match='No abundance columns'):
            prepare_abundance(matrix, params)


class TestBuildGroupMap:

    def test_maps_sample_to_group(self):
        meta = pd.DataFrame({'sample_id': ['S1', 'S2'], 'group': ['WT', 'KO']})
        params = {'design': {'group_column': 'group'}}
        assert build_group_map(meta, params) == {'S1': 'WT', 'S2': 'KO'}

    def test_missing_group_column_raises(self):
        meta = pd.DataFrame({'sample_id': ['S1'], 'condition': ['WT']})
        params = {'design': {'group_column': 'group'}}
        with pytest.raises(ValueError, match="Group column 'group' not found"):
            build_group_map(meta, params)


# ============================================================
# Section 7: parametrized invariant sweep
# ============================================================
# Runs the median-equalization invariant across several matrix shapes and seeds
# in one place, so adding a shape is a one-line change.

class TestParametrizedInvariants:

    @pytest.mark.parametrize('n_proteins,n_samples,seed', [
        (10, 2, 1),
        (50, 4, 2),
        (200, 6, 3),
        (5, 3, 4),
    ])
    def test_median_invariant_across_shapes(self, n_proteins, n_samples, seed):
        rng = np.random.default_rng(seed=seed)
        df = _log2_frame({
            f'S{j}': rng.normal(loc=5 * j, scale=2, size=n_proteins).tolist()
            for j in range(n_samples)
        })
        out = normalize_median(df)
        medians = out.median().values
        assert np.allclose(medians, medians[0])

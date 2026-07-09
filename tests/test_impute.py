#!/usr/bin/env python3
# title: test_impute.py
# project: ProSIFT (PROtein Statistical Integration and Filtering Tool)
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-07-09
# last modified: 2026-07-09
#
# purpose:
#   Unit tests for Module 03 IMPUTE (bin/impute.py). Scope: the MNAR/MAR
#   classification, the three imputation methods, and the dispatcher -- the
#   logic that decides which missing values get MinProb (below-detection) vs
#   KNN (technical dropout) imputation, and fills them.
#
#   classify_missingness is the scientific centerpiece (the analogue of Module
#   02's compute_sample_flags): a bug there silently changes which values are
#   treated as below-detection vs random, which changes every downstream
#   statistic. Its cases mirror the committed benchmark's edge proteins
#   (EDGE_MNAR_WT/KO -> SINGLE-GROUP, EDGE_PARTIAL -> PARTIAL, EDGE_MAR ->
#   PASSED), so the in-memory cases here match the validated oracle.
#
#   Test kinds (same taxonomy as test_normalize.py):
#     1. GROUND-TRUTH -- hand-built inputs with a known correct classification.
#     2. INVARIANT    -- guarantees that hold for any input (no NaN remains,
#                        observed values never change, masks disjoint).
#     3. CONTRACT     -- the dispatcher rejects bad config loudly.
#     4. EDGE         -- no-missing, protein absent from filter table, and an
#                        adversarial single-observation probe (xfail: see below).
#
#   Not tested here: Plotly plots, output writers, main() (I/O / integration).
#
# inputs:
#   None (tests build inputs in-memory).
#
# outputs:
#   Test results (stdout via pytest). Note: one xfail documents a real finding
#   (impute_minprob single-observation guard) rather than a test failure.
#
# usage example:
#   pytest tests/test_impute.py -v
#   pytest tests/test_impute.py::TestClassifyMissingness -v
#
#   copy/paste: pytest tests/test_impute.py -v

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# Add bin/ to path so we can import impute directly (mirrors other test files).
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'bin'))

from impute import (
    build_group_map,
    build_imputation_mask,
    build_imputation_summary,
    classify_missingness,
    impute_knn,
    impute_left_censored,
    impute_minprob,
    prepare_normalized_abundance,
    run_imputation,
)

# ============================================================
# SHARED LAYOUT + HELPERS
# ============================================================
# 3 WT + 3 KO, matching the benchmark's design so classification cases here
# read against the same edge scenarios the benchmark encodes.

SAMPLES = ['WT-1', 'WT-2', 'WT-3', 'KO-1', 'KO-2', 'KO-3']
GROUP_MAP = {s: ('WT' if s.startswith('WT') else 'KO') for s in SAMPLES}


def _norm(data: dict, protein_ids) -> pd.DataFrame:
    '''Build a normalized (log2) matrix: index=protein_id, columns=sample_ids.'''
    return pd.DataFrame(data, index=pd.Index(protein_ids, name='protein_id'))


def _filter(status_map: dict) -> pd.DataFrame:
    '''Build a detection filter table with protein_id + filter_status.'''
    return pd.DataFrame({
        'protein_id': list(status_map),
        'filter_status': list(status_map.values()),
    })


def _imp_params(**overrides) -> dict:
    '''Minimal params dict covering what the imputation functions read.'''
    imp = {
        'mode': 'mixed',
        'mnar_method': 'minprob',
        'mar_method': 'knn',
        'single_method': 'minprob',
        'minprob_quantile': 0.01,
        'minprob_scale': 0.3,
        'knn_k': 3,
        'left_censored_downshift': 1.8,
        'left_censored_width': 0.3,
        'random_seed': 42,
    }
    imp.update(overrides)
    return {
        'input': {'abundance_prefix': '', 'peptide_count_prefix': 'peptide_count_'},
        'design': {'group_column': 'group'},
        'imputation': imp,
    }


# ============================================================
# Section 1: classify_missingness (scientific centerpiece)
# ============================================================

class TestClassifyMissingness:

    def test_passed_sporadic_is_all_mar(self):
        '''PASSED protein: any missing value is sporadic technical dropout -> MAR.'''
        norm = _norm({
            'WT-1': [20.0], 'WT-2': [np.nan], 'WT-3': [20.0],
            'KO-1': [20.0], 'KO-2': [20.0], 'KO-3': [20.0],
        }, ['P1'])
        mnar, mar, cls = classify_missingness(norm, _filter({'P1': 'PASSED'}),
                                              GROUP_MAP, SAMPLES)
        assert not mnar.loc['P1'].any()
        assert bool(mar.loc['P1', 'WT-2']) is True
        assert int(mar.loc['P1'].sum()) == 1
        assert cls['P1'] == 'MAR'

    def test_single_group_absent_is_mnar_detected_sporadic_is_mar(self):
        '''
        SINGLE-GROUP: the fully-absent group's NaN are below-detection (MNAR);
        a sporadic NaN in the detected group is MAR. protein_class = MNAR.
        Mirrors benchmark EDGE_MNAR_WT.
        '''
        norm = _norm({
            'WT-1': [np.nan], 'WT-2': [np.nan], 'WT-3': [np.nan],   # absent group
            'KO-1': [20.0], 'KO-2': [np.nan], 'KO-3': [20.0],       # sporadic gap
        }, ['P1'])
        mnar, mar, cls = classify_missingness(norm, _filter({'P1': 'SINGLE-GROUP'}),
                                              GROUP_MAP, SAMPLES)
        assert mnar.loc['P1', ['WT-1', 'WT-2', 'WT-3']].all()
        assert bool(mar.loc['P1', 'KO-2']) is True
        assert int(mnar.loc['P1'].sum()) == 3
        assert int(mar.loc['P1'].sum()) == 1
        assert cls['P1'] == 'MNAR'

    def test_partial_subthreshold_mnar_passing_sporadic_mar_is_mixed(self):
        '''
        PARTIAL: the sub-threshold (minority-detected) group's NaN are MNAR; the
        passing (majority-detected) group's sporadic NaN are MAR. Both present ->
        protein_class = "mixed". Mirrors benchmark EDGE_PARTIAL.
          WT: 1 of 3 detected  -> 1 < 3/2 -> MNAR for WT-2, WT-3
          KO: 2 of 3 detected  -> 2 >= 3/2 -> MAR for KO-3
        '''
        norm = _norm({
            'WT-1': [20.0], 'WT-2': [np.nan], 'WT-3': [np.nan],
            'KO-1': [20.0], 'KO-2': [20.0], 'KO-3': [np.nan],
        }, ['P1'])
        mnar, mar, cls = classify_missingness(norm, _filter({'P1': 'PARTIAL'}),
                                              GROUP_MAP, SAMPLES)
        assert mnar.loc['P1', ['WT-2', 'WT-3']].all()
        assert bool(mar.loc['P1', 'KO-3']) is True
        assert cls['P1'] == 'mixed'

    def test_partial_heuristic_boundary_two_of_three_is_mar(self):
        '''
        Heuristic boundary (n_detected < n_total/2): a group with 2 of 3 detected
        is NOT sub-threshold (2 >= 1.5), so its single gap is MAR, not MNAR.
        '''
        norm = _norm({
            'WT-1': [20.0], 'WT-2': [20.0], 'WT-3': [np.nan],
            'KO-1': [20.0], 'KO-2': [20.0], 'KO-3': [20.0],
        }, ['P1'])
        mnar, mar, cls = classify_missingness(norm, _filter({'P1': 'PARTIAL'}),
                                              GROUP_MAP, SAMPLES)
        assert not mnar.loc['P1'].any()
        assert bool(mar.loc['P1', 'WT-3']) is True
        assert cls['P1'] == 'MAR'

    def test_no_missing_not_classified(self):
        '''A fully-observed protein needs no imputation and stays class "none".'''
        norm = _norm({s: [20.0] for s in SAMPLES}, ['P1'])
        mnar, mar, cls = classify_missingness(norm, _filter({'P1': 'PASSED'}),
                                              GROUP_MAP, SAMPLES)
        assert not mnar.loc['P1'].any()
        assert not mar.loc['P1'].any()
        assert cls['P1'] == 'none'

    def test_protein_absent_from_filter_table_is_skipped(self):
        '''A protein with missing values but no filter-table row is left
        unclassified (defensive: should not happen in the real pipeline).'''
        norm = _norm({
            'WT-1': [np.nan], 'WT-2': [20.0], 'WT-3': [20.0],
            'KO-1': [20.0], 'KO-2': [20.0], 'KO-3': [20.0],
        }, ['P1'])
        mnar, mar, cls = classify_missingness(norm, _filter({'OTHER': 'PASSED'}),
                                              GROUP_MAP, SAMPLES)
        assert not mnar.loc['P1'].any()
        assert not mar.loc['P1'].any()
        assert cls['P1'] == 'none'


# ============================================================
# Section 2: impute_minprob
# ============================================================

class TestImputeMinProb:

    def test_only_target_positions_filled_observed_untouched(self):
        '''
        MinProb fills only positions flagged in target_mask; NaN outside the
        mask remain NaN (KNN fills those later), and observed values are unchanged.
        '''
        norm = _norm({
            'S1': [10.0, np.nan, 12.0, np.nan, 14.0],
        }, ['P1', 'P2', 'P3', 'P4', 'P5'])
        target = pd.DataFrame(False, index=norm.index, columns=norm.columns)
        target.loc['P2', 'S1'] = True                 # this NaN is targeted
        # P4/S1 is NaN but NOT targeted -> must remain NaN

        out = impute_minprob(norm, target, _imp_params(), np.random.default_rng(0))
        assert not pd.isna(out.loc['P2', 'S1'])        # filled
        assert pd.isna(out.loc['P4', 'S1'])            # untargeted NaN preserved
        assert out.loc['P1', 'S1'] == 10.0             # observed unchanged
        assert out.loc['P5', 'S1'] == 14.0

    def test_determinism_same_seed(self):
        '''Same seed -> identical imputed values (seeded rng).'''
        norm = _norm({'S1': [10.0, np.nan, 12.0, np.nan, 14.0]},
                     ['P1', 'P2', 'P3', 'P4', 'P5'])
        target = norm.isna()
        out1 = impute_minprob(norm, target, _imp_params(), np.random.default_rng(7))
        out2 = impute_minprob(norm, target, _imp_params(), np.random.default_rng(7))
        pd.testing.assert_frame_equal(out1, out2)

    def test_imputed_values_are_left_shifted(self):
        '''
        Property: MinProb centers on a low quantile, so imputed MNAR values sit
        below the observed distribution -- their mean is below the observed
        median. This is the defining behavior of below-detection imputation.
        '''
        rng_data = np.random.default_rng(0)
        col = list(rng_data.normal(loc=20.0, scale=2.0, size=20))
        col[0] = col[1] = col[2] = np.nan
        norm = _norm({'S1': col}, [f'P{i}' for i in range(20)])
        target = norm.isna()

        out = impute_minprob(norm, target, _imp_params(), np.random.default_rng(1))
        imputed = out.loc[['P0', 'P1', 'P2'], 'S1'].to_numpy()
        observed = np.array([v for v in col if not np.isnan(v)])
        assert imputed.mean() < np.median(observed)

    @pytest.mark.xfail(
        reason=(
            'FINDING (2026-07-09): impute_minprob guards len(observed)==0 but not '
            'len<2. With exactly one observed value in a sample, np.std(ddof=1) '
            'divides by zero -> NaN width -> NaN imputed value (and a RuntimeWarning '
            'that this project escalates to an error). The guard should be len<2. '
            'Low real-data risk (a sample column spans all proteins, so it is never '
            'down to one observation), but reachable on small/degenerate inputs.'
        ),
        strict=False,
    )
    def test_single_observed_value_still_imputes_without_nan(self):
        '''Adversarial probe: one observed value in a sample. Correct behavior is
        no NaN left behind; current behavior raises / yields NaN.'''
        norm = _norm({'S1': [10.0, np.nan, np.nan]}, ['P1', 'P2', 'P3'])
        target = norm.isna()
        out = impute_minprob(norm, target, _imp_params(), np.random.default_rng(0))
        assert int(out.isna().sum().sum()) == 0


# ============================================================
# Section 3: impute_knn
# ============================================================

class TestImputeKNN:

    def _knn_frame(self):
        # 5 proteins x 3 samples, scattered NaN, no all-NaN row or column.
        return _norm({
            'S1': [10.0, 11.0, np.nan, 13.0, 14.0],
            'S2': [10.5, np.nan, 12.5, 13.5, 14.5],
            'S3': [np.nan, 11.2, 12.2, 13.2, 14.2],
        }, ['P1', 'P2', 'P3', 'P4', 'P5'])

    def test_fills_all_nan(self):
        norm = self._knn_frame()
        out = impute_knn(norm, norm.isna(), _imp_params())
        assert int(out.isna().sum().sum()) == 0

    def test_observed_values_unchanged(self):
        norm = self._knn_frame()
        out = impute_knn(norm, norm.isna(), _imp_params())
        assert out.loc['P1', 'S1'] == 10.0             # was observed
        assert out.loc['P4', 'S2'] == 13.5


# ============================================================
# Section 4: impute_left_censored (single mode)
# ============================================================

class TestImputeLeftCensored:

    def test_imputed_below_observed_mean(self):
        '''Left-censored centers at mean - downshift*SD, so imputed < observed mean.'''
        rng_data = np.random.default_rng(0)
        col = list(rng_data.normal(loc=20.0, scale=2.0, size=20))
        col[0] = col[1] = np.nan
        norm = _norm({'S1': col}, [f'P{i}' for i in range(20)])
        target = norm.isna()

        out = impute_left_censored(norm, target, _imp_params(), np.random.default_rng(1))
        imputed = out.loc[['P0', 'P1'], 'S1'].to_numpy()
        observed = np.array([v for v in col if not np.isnan(v)])
        assert imputed.mean() < observed.mean()


# ============================================================
# Section 5: run_imputation dispatcher (contract + completeness)
# ============================================================

class TestRunImputation:

    def _mixed_frame(self):
        # PASS proteins with sporadic gaps (MAR) + one SINGLE-GROUP (MNAR).
        norm = _norm({
            'WT-1': [20.0, np.nan, 20.0, 20.0, 20.0],
            'WT-2': [20.0, np.nan, 20.0, np.nan, 20.0],
            'WT-3': [20.0, np.nan, 20.0, 20.0, 20.0],
            'KO-1': [20.0, 20.0, 20.0, 20.0, 20.0],
            'KO-2': [20.0, 20.0, 20.0, 20.0, 20.0],
            'KO-3': [20.0, 20.0, 20.0, 20.0, 20.0],
        }, ['PASS1', 'SG1', 'PASS2', 'PASS3', 'PASS4'])
        filt = _filter({'PASS1': 'PASSED', 'SG1': 'SINGLE-GROUP',
                        'PASS2': 'PASSED', 'PASS3': 'PASSED', 'PASS4': 'PASSED'})
        return norm, filt

    def test_unknown_mode_raises(self):
        norm, filt = self._mixed_frame()
        with pytest.raises(ValueError, match='Unknown imputation'):
            run_imputation(norm, filt, GROUP_MAP, SAMPLES,
                           _imp_params(mode='bogus'), np.random.default_rng(0))

    def test_unknown_mnar_method_raises(self):
        norm, filt = self._mixed_frame()
        with pytest.raises(ValueError, match='mnar_method'):
            run_imputation(norm, filt, GROUP_MAP, SAMPLES,
                           _imp_params(mnar_method='bogus'), np.random.default_rng(0))

    def test_unknown_mar_method_raises(self):
        norm, filt = self._mixed_frame()
        with pytest.raises(ValueError, match='mar_method'):
            run_imputation(norm, filt, GROUP_MAP, SAMPLES,
                           _imp_params(mar_method='bogus'), np.random.default_rng(0))

    def test_unknown_single_method_raises(self):
        norm, filt = self._mixed_frame()
        with pytest.raises(ValueError, match='single_method'):
            run_imputation(norm, filt, GROUP_MAP, SAMPLES,
                           _imp_params(mode='single', single_method='bogus'),
                           np.random.default_rng(0))

    def test_mixed_leaves_no_nan(self):
        '''Core guarantee: after mixed imputation, no NaN remain.'''
        norm, filt = self._mixed_frame()
        imputed, _mnar, _mar, _cls, _m = run_imputation(
            norm, filt, GROUP_MAP, SAMPLES, _imp_params(), np.random.default_rng(0)
        )
        assert int(imputed.isna().sum().sum()) == 0

    def test_single_minprob_leaves_no_nan(self):
        norm, filt = self._mixed_frame()
        imputed, _mnar, _mar, _cls, _m = run_imputation(
            norm, filt, GROUP_MAP, SAMPLES,
            _imp_params(mode='single', single_method='minprob'),
            np.random.default_rng(0),
        )
        assert int(imputed.isna().sum().sum()) == 0

    def test_single_knn_leaves_no_nan(self):
        norm, filt = self._mixed_frame()
        imputed, _mnar, _mar, _cls, _m = run_imputation(
            norm, filt, GROUP_MAP, SAMPLES,
            _imp_params(mode='single', single_method='knn'),
            np.random.default_rng(0),
        )
        assert int(imputed.isna().sum().sum()) == 0


# ============================================================
# Section 6: cross-cutting invariants
# ============================================================

class TestInvariants:

    def _run(self):
        # 4 proteins with distinct values so each WT column keeps >= 2 observed
        # values even with SG1 absent across WT (avoids the single-observation
        # MinProb std issue that TestImputeMinProb documents separately).
        norm = _norm({
            'WT-1': [20.0, np.nan, 21.0, 22.0],
            'WT-2': [20.1, np.nan, np.nan, 22.1],   # PASS2 sporadic gap (MAR)
            'WT-3': [20.2, np.nan, 21.2, 22.2],
            'KO-1': [20.0, 20.5, 21.0, 22.0],
            'KO-2': [20.1, 20.6, 21.1, 22.1],
            'KO-3': [20.2, 20.7, 21.2, 22.2],
        }, ['PASS1', 'SG1', 'PASS2', 'PASS3'])
        filt = _filter({'PASS1': 'PASSED', 'SG1': 'SINGLE-GROUP',
                        'PASS2': 'PASSED', 'PASS3': 'PASSED'})
        return norm, run_imputation(norm, filt, GROUP_MAP, SAMPLES,
                                    _imp_params(), np.random.default_rng(0))

    def test_observed_values_unchanged(self):
        '''Imputation must never alter an observed value.'''
        norm, (imputed, _mnar, _mar, _cls, _m) = self._run()
        obs = (~norm.isna()).to_numpy()
        assert np.allclose(imputed.to_numpy()[obs], norm.to_numpy()[obs])

    def test_mnar_mar_masks_disjoint(self):
        '''A cell cannot be classified as both MNAR and MAR.'''
        _norm_df, (_imp, mnar, mar, _cls, _m) = self._run()
        assert not (mnar & mar).to_numpy().any()

    def test_build_imputation_mask_matches_masks(self):
        '''Categorical mask reflects the boolean masks: observed/mnar/mar.'''
        norm, (_imp, mnar, mar, _cls, _m) = self._run()
        mask = build_imputation_mask(norm, mnar, mar)
        assert (mask.to_numpy()[mnar.to_numpy()] == 'mnar').all()
        assert (mask.to_numpy()[mar.to_numpy()] == 'mar').all()
        observed = ~(mnar | mar)
        assert (mask.to_numpy()[observed.to_numpy()] == 'observed').all()


# ============================================================
# Section 7: prepare_normalized_abundance + build_group_map
# ============================================================

class TestPrepareAndGroupMap:

    def test_strips_prefix_and_separates_peptides(self):
        matrix = pd.DataFrame({
            'protein_id': ['P1', 'P2'],
            'abundance_S1': [10.0, 20.0],
            'abundance_S2': [30.0, 40.0],
            'peptide_count_S1': [3, 5],
            'peptide_count_S2': [4, 6],
        })
        params = {'input': {'abundance_prefix': 'abundance_',
                            'peptide_count_prefix': 'peptide_count_'}}
        norm_df, peptide_df, sample_ids, prefix = prepare_normalized_abundance(matrix, params)
        assert sample_ids == ['S1', 'S2']
        assert list(norm_df.columns) == ['S1', 'S2']
        assert list(peptide_df.columns) == ['peptide_count_S1', 'peptide_count_S2']
        assert prefix == 'abundance_'

    def test_no_abundance_columns_raises(self):
        matrix = pd.DataFrame({'protein_id': ['P1'], 'peptide_count_S1': [3]})
        params = {'input': {'abundance_prefix': '', 'peptide_count_prefix': 'peptide_count_'}}
        with pytest.raises(ValueError, match='No abundance columns'):
            prepare_normalized_abundance(matrix, params)

    def test_missing_group_column_raises(self):
        meta = pd.DataFrame({'sample_id': ['S1'], 'condition': ['WT']})
        params = {'design': {'group_column': 'group'}}
        with pytest.raises(ValueError, match="Group column 'group' not found"):
            build_group_map(meta, params)


# ============================================================
# Section 8: build_imputation_summary
# ============================================================

class TestImputationSummary:

    def test_per_protein_counts_match_masks(self):
        norm = _norm({s: [20.0, 20.0] for s in SAMPLES}, ['P1', 'P2'])
        filt = _filter({'P1': 'SINGLE-GROUP', 'P2': 'PASSED'})
        mnar = pd.DataFrame(False, index=norm.index, columns=norm.columns)
        mar = pd.DataFrame(False, index=norm.index, columns=norm.columns)
        mnar.loc['P1', ['WT-1', 'WT-2']] = True   # 2 MNAR
        mar.loc['P2', 'KO-3'] = True               # 1 MAR
        cls = pd.Series({'P1': 'MNAR', 'P2': 'MAR'})

        summary = build_imputation_summary(norm, filt, cls, mnar, mar).set_index('protein_id')
        assert summary.loc['P1', 'n_mnar_imputed'] == 2
        assert summary.loc['P1', 'n_imputed_total'] == 2
        assert summary.loc['P2', 'n_mar_imputed'] == 1
        assert summary.loc['P1', 'filter_status'] == 'SINGLE-GROUP'

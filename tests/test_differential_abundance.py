#!/usr/bin/env python3
# title: test_differential_abundance.py
# project: ProSIFT (PROtein Statistical Integration and Filtering Tool)
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-07-09
# last modified: 2026-07-09
#
# purpose:
#   Unit tests for Module 04 DIFFERENTIAL_ABUNDANCE (bin/differential_abundance.py).
#   The statistical fit itself (limma + DEqMS) runs in R via rpy2 and is an
#   integration concern for a cluster run (see NOTE below); it is not exercised
#   here. The bulk of the module's decision logic is pure Python and is fully
#   covered:
#     - parse_and_validate_contrasts : contrast-string parsing + validation
#     - summarize_peptide_counts      : per-protein min-nonzero peptide count
#                                       (the DEqMS variance predictor)
#     - assemble_results              : R-output -> ProSIFT schema, significance
#                                       call (FDR + FC), direction assignment
#     - extract_abundance_and_peptide_cols / build_group_map : I/O parsing
#
#   assemble_results is the centerpiece: it decides which proteins are called
#   significant and their direction. We fabricate an R-output-shaped frame (no R
#   needed) and check the calls exactly.
#
#   Not tested here: the plot_* Plotly functions, write_summary_txt, main, and
#   the R fit (_run_one_contrast_r / limma+DEqMS). NOTE: importing this module
#   starts embedded R via rpy2, which segfaults where R is not linked; the test
#   neutralizes that import (see below), so the R fit must be validated by a
#   cluster integration run, not this unit file.
#
# inputs:
#   None (tests build inputs in-memory).
#
# outputs:
#   Test results (stdout via pytest).
#
# usage example:
#   pytest tests/test_differential_abundance.py -v
#
#   copy/paste: pytest tests/test_differential_abundance.py -v

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# Add bin/ to path.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'bin'))

# differential_abundance imports rpy2.robjects at module load, which starts an
# embedded R. Where R is not properly linked (e.g. a dev laptop with a broken
# rpy2/R), that import SEGFAULTS -- and the module's `except ImportError` guard
# cannot catch a native crash. Neutralize the rpy2 import so the module loads
# with _HAVE_RPY2=False (the pure functions under test do not use rpy2), then
# restore sys.modules so other test files are unaffected. The real limma/DEqMS
# fit is exercised by cluster integration runs, not this unit file.
_RPY2_SENTINEL = object()
_rpy2_saved = sys.modules.get('rpy2', _RPY2_SENTINEL)
sys.modules['rpy2'] = None
try:
    from differential_abundance import (
        assemble_results,
        build_group_map,
        extract_abundance_and_peptide_cols,
        parse_and_validate_contrasts,
        summarize_peptide_counts,
    )
finally:
    if _rpy2_saved is _RPY2_SENTINEL:
        sys.modules.pop('rpy2', None)
    else:
        sys.modules['rpy2'] = _rpy2_saved


# ============================================================
# SHARED HELPERS
# ============================================================

def _metadata():
    return pd.DataFrame({
        'sample_id': ['WT-1', 'WT-2', 'WT-3', 'KO-1', 'KO-2', 'KO-3'],
        'genotype': ['WT', 'WT', 'WT', 'KO', 'KO', 'KO'],
    })


def _params(contrasts, fdr=0.05, fc=1.0):
    return {
        'design': {'group_column': 'genotype', 'contrasts': contrasts},
        'differential_abundance': {'significance': {'fdr_threshold': fdr,
                                                    'fc_threshold': fc}},
    }


# ============================================================
# Section 1: parse_and_validate_contrasts
# ============================================================

class TestParseContrasts:

    def test_valid_contrast_parsed(self):
        parsed = parse_and_validate_contrasts(_params(['KO_vs_WT']), _metadata())
        assert parsed == [('KO_vs_WT', 'KO', 'WT', 'KO - WT')]

    def test_multiple_contrasts(self):
        parsed = parse_and_validate_contrasts(
            _params(['KO_vs_WT', 'WT_vs_KO']), _metadata())
        assert len(parsed) == 2
        assert parsed[1] == ('WT_vs_KO', 'WT', 'KO', 'WT - KO')

    def test_missing_delimiter_raises(self):
        with pytest.raises(ValueError, match="_vs_"):
            parse_and_validate_contrasts(_params(['KOWT']), _metadata())

    def test_empty_denominator_raises(self):
        with pytest.raises(ValueError, match='denominator is empty'):
            parse_and_validate_contrasts(_params(['KO_vs_']), _metadata())

    def test_unknown_group_raises(self):
        with pytest.raises(ValueError, match='not found'):
            parse_and_validate_contrasts(_params(['KO_vs_XX']), _metadata())

    def test_no_contrasts_raises(self):
        with pytest.raises(ValueError, match='No contrasts defined'):
            parse_and_validate_contrasts(_params([]), _metadata())


# ============================================================
# Section 2: summarize_peptide_counts (DEqMS variance predictor)
# ============================================================

class TestSummarizePeptideCounts:

    def test_min_of_nonzero_excludes_zeros(self):
        '''Per protein: minimum of nonzero counts; zeros (imputed positions) skip.'''
        pep = pd.DataFrame(
            {'s1': [3, 0], 's2': [0, 7], 's3': [5, 2]},
            index=pd.Index(['P1', 'P2'], name='protein_id'),
        )
        out = summarize_peptide_counts(pep)
        assert out['P1'] == 3          # min(3, 5), the 0 excluded
        assert out['P2'] == 2          # min(7, 2)

    def test_all_positive(self):
        pep = pd.DataFrame({'s1': [4], 's2': [5], 's3': [6]},
                           index=pd.Index(['P1'], name='protein_id'))
        assert summarize_peptide_counts(pep)['P1'] == 4

    def test_returns_int64(self):
        pep = pd.DataFrame({'s1': [3], 's2': [5]},
                           index=pd.Index(['P1'], name='protein_id'))
        assert summarize_peptide_counts(pep).dtype == np.int64

    @pytest.mark.xfail(
        reason=(
            'FINDING (2026-07-09): summarize_peptide_counts assumes every protein '
            'has >= 1 nonzero peptide count (Module 01 removes ABSENT proteins). '
            'An all-zero row hits np.nanmin over an all-NaN slice -> NaN + '
            'RuntimeWarning (escalated to error under pytest); in production the '
            'NaN.astype(int64) yields a garbage sentinel (-2^63) silently. '
            'Low reachability (a retained protein has detections, hence peptides), '
            'but unguarded. A len/all-zero guard would make it explicit.'
        ),
        strict=False,
    )
    def test_all_zero_row_does_not_corrupt(self):
        '''Adversarial probe: an all-zero peptide row. Correct behavior is a
        sensible value or a clear error, not a silent garbage int.'''
        pep = pd.DataFrame({'s1': [0], 's2': [0]},
                           index=pd.Index(['P1'], name='protein_id'))
        out = summarize_peptide_counts(pep)
        assert out['P1'] >= 0          # not a -2^63 garbage sentinel


# ============================================================
# Section 3: assemble_results (significance + direction centerpiece)
# ============================================================

class TestAssembleResults:

    def _deqms_raw(self):
        '''An R-output-shaped frame as DEqMS would return (no R needed).'''
        return pd.DataFrame({
            'protein_id':   ['P1', 'P2', 'P3', 'P4'],
            'logFC':        [2.0, -1.5, 2.0, 0.5],
            'AveExpr':      [20.0, 21.0, 19.0, 22.0],
            't':            [5.0, -4.0, 3.0, 4.0],
            'P.Value':      [1e-4, 1e-3, 0.10, 1e-3],
            'adj.P.Val':    [1e-3, 5e-3, 0.20, 5e-3],
            'B':            [3.0, 2.0, -1.0, 2.0],
            'sca.t':        [5.5, -4.5, 3.2, 4.2],
            'sca.P.Value':  [8e-5, 8e-4, 0.09, 8e-4],
            'sca.adj.pval': [1e-2, 1e-3, 0.20, 1e-2],
            'count':        [3, 5, 2, 4],
        })

    def _id_mapping(self):
        return pd.DataFrame({'protein_id': ['P1', 'P2', 'P3', 'P4'],
                             'gene_symbol': ['Gene1', 'Gene2', 'Gene3', 'Gene4']})

    def test_schema_and_column_order(self):
        out = assemble_results(self._deqms_raw(), self._id_mapping(),
                               _params(['KO_vs_WT']), 'DEqMS', 'KO_vs_WT')
        assert list(out.columns) == [
            'protein_id', 'gene_symbol', 'log2_fc', 'avg_abundance',
            'limma_t', 'limma_pvalue', 'limma_adj_pvalue',
            'deqms_t', 'deqms_pvalue', 'deqms_adj_pvalue',
            'n_peptides', 'significant', 'direction', 'contrast',
        ]
        assert (out['contrast'] == 'KO_vs_WT').all()
        assert out.set_index('protein_id').loc['P1', 'n_peptides'] == 3

    def test_significance_and_direction_calls(self):
        '''
        With FDR<0.05 and |FC|>1, using the DEqMS adj p-value:
          P1: adj 0.01, FC +2   -> significant, up
          P2: adj 0.001, FC -1.5 -> significant, down
          P3: adj 0.20         -> ns (fails FDR)
          P4: adj 0.01, FC 0.5  -> ns (fails FC)
        '''
        out = assemble_results(self._deqms_raw(), self._id_mapping(),
                               _params(['KO_vs_WT']), 'DEqMS', 'KO_vs_WT'
                               ).set_index('protein_id')
        assert bool(out.loc['P1', 'significant']) is True
        assert out.loc['P1', 'direction'] == 'up'
        assert out.loc['P2', 'direction'] == 'down'
        assert bool(out.loc['P3', 'significant']) is False
        assert out.loc['P3', 'direction'] == 'ns'
        assert bool(out.loc['P4', 'significant']) is False      # fails FC
        assert out.loc['P4', 'direction'] == 'ns'

    def test_significance_uses_deqms_pvalue_when_deqms(self):
        '''
        A protein whose DEqMS adj p-value passes but whose limma adj p-value fails
        must be called significant (DEqMS is the primary p-value for method=DEqMS).
        '''
        raw = self._deqms_raw()
        raw.loc[raw['protein_id'] == 'P1', 'adj.P.Val'] = 0.30       # limma fails
        raw.loc[raw['protein_id'] == 'P1', 'sca.adj.pval'] = 0.01    # DEqMS passes
        out = assemble_results(raw, self._id_mapping(), _params(['KO_vs_WT']),
                               'DEqMS', 'KO_vs_WT').set_index('protein_id')
        assert bool(out.loc['P1', 'significant']) is True

    def test_limma_only_fills_deqms_columns_na(self):
        '''method=limma: DEqMS columns and n_peptides are NA; significance uses
        the limma adj p-value.'''
        raw = self._deqms_raw()[['protein_id', 'logFC', 'AveExpr', 't',
                                 'P.Value', 'adj.P.Val', 'B']]
        out = assemble_results(raw, self._id_mapping(), _params(['KO_vs_WT']),
                               'limma', 'KO_vs_WT').set_index('protein_id')
        assert pd.isna(out.loc['P1', 'deqms_t'])
        assert pd.isna(out.loc['P1', 'n_peptides'])
        assert bool(out.loc['P1', 'significant']) is True           # limma adj 1e-3

    def test_fc_threshold_zero_disables_fc_filter(self):
        '''fc_threshold=0 -> only FDR matters; a small-FC protein passing FDR is
        significant.'''
        raw = self._deqms_raw()
        out = assemble_results(raw, self._id_mapping(),
                               _params(['KO_vs_WT'], fc=0.0), 'DEqMS', 'KO_vs_WT'
                               ).set_index('protein_id')
        assert bool(out.loc['P4', 'significant']) is True          # FC no longer filters

    def test_gene_symbol_mouse_column_preferred(self):
        '''gene_symbol_mouse is used when present (ortholog-mapped runs).'''
        idmap = pd.DataFrame({'protein_id': ['P1', 'P2', 'P3', 'P4'],
                              'gene_symbol_mouse': ['M1', 'M2', 'M3', 'M4']})
        out = assemble_results(self._deqms_raw(), idmap, _params(['KO_vs_WT']),
                               'DEqMS', 'KO_vs_WT').set_index('protein_id')
        assert out.loc['P1', 'gene_symbol'] == 'M1'

    def test_missing_gene_column_yields_na(self):
        idmap = pd.DataFrame({'protein_id': ['P1', 'P2', 'P3', 'P4']})
        out = assemble_results(self._deqms_raw(), idmap, _params(['KO_vs_WT']),
                               'DEqMS', 'KO_vs_WT').set_index('protein_id')
        assert pd.isna(out.loc['P1', 'gene_symbol'])


# ============================================================
# Section 4: I/O parsing
# ============================================================

class TestIOParsing:

    def test_splits_abundance_and_peptides(self):
        matrix = pd.DataFrame({
            'protein_id': ['P1', 'P2'],
            'abundance_WT-1': [20.0, 21.0],
            'abundance_KO-1': [22.0, 23.0],
            'peptide_count_WT-1': [3, 5],
            'peptide_count_KO-1': [4, 6],
        })
        params = {'input': {'abundance_prefix': 'abundance_'}}
        abund, pep, sample_ids = extract_abundance_and_peptide_cols(matrix, params)
        assert sample_ids == ['WT-1', 'KO-1']
        assert list(abund.columns) == ['WT-1', 'KO-1']
        assert list(pep.columns) == ['peptide_count_WT-1', 'peptide_count_KO-1']

    def test_no_abundance_columns_raises(self):
        matrix = pd.DataFrame({'protein_id': ['P1'], 'peptide_count_WT-1': [3]})
        with pytest.raises(ValueError, match='No abundance columns'):
            extract_abundance_and_peptide_cols(matrix, {'input': {}})

    def test_build_group_map_missing_column_raises(self):
        meta = pd.DataFrame({'sample_id': ['WT-1'], 'condition': ['WT']})
        with pytest.raises(ValueError, match="group_column 'genotype' not found"):
            build_group_map(meta, {'design': {'group_column': 'genotype'}})

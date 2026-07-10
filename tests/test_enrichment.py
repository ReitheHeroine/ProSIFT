#!/usr/bin/env python3
# title: test_enrichment.py
# project: ProSIFT (PROtein Statistical Integration and Filtering Tool)
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-07-09
# last modified: 2026-07-09
#
# purpose:
#   Unit tests for Module 05 ENRICHMENT (bin/enrichment.py). The ORA/GSEA calls
#   (gseapy) and the GO-term redundancy reduction (rrvgo via R) are integration
#   concerns; the module's pure-Python data logic is covered here:
#     - _library_short_name             : GMT path -> short library id
#     - prepare_gene_symbols            : per-contrast filter, drop unmapped,
#                                         dedup by gene (keep most significant)
#     - build_ranked_series             : GSEA ranking metric (3 modes)
#     - build_protein_term_mapping      : GMT parse -> many-to-many protein/term
#     - _truncate_label / _msigdb_name_to_go_lookup_phrase : string helpers
#     - cluster_go_terms (no-rpy2 path) : graceful degradation to null columns
#
#   Not tested here: run_ora / run_gsea (gseapy), the plot_* functions, main,
#   and the real rrvgo clustering (needs R + rrvgo + GO.db + org db -> a cluster
#   integration test, like Module 04's R fit).
#
#   enrichment imports cleanly: rpy2 is loaded lazily inside cluster_go_terms, so
#   importing the module does not start embedded R. The cluster_go_terms no-rpy2
#   degradation test blocks rpy2 for the duration of the call (via monkeypatch)
#   so the lazy import fails with ImportError rather than starting R.
#
# inputs:
#   None (tests build inputs in-memory / in tmp_path).
#
# outputs:
#   Test results (stdout via pytest).
#
# usage example:
#   pytest tests/test_enrichment.py -v
#
#   copy/paste: pytest tests/test_enrichment.py -v

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'bin'))

from enrichment import (
    _library_short_name,
    _msigdb_name_to_go_lookup_phrase,
    _truncate_label,
    build_protein_term_mapping,
    build_ranked_series,
    cluster_go_terms,
    prepare_gene_symbols,
)

# ============================================================
# Section 1: _library_short_name
# ============================================================

class TestLibraryShortName:

    def test_go_bp(self):
        assert _library_short_name('gmt/m5.go.bp.v2026.1.Mm.symbols.gmt') == 'GO_BP'

    def test_reactome(self):
        assert _library_short_name('gmt/m2.cp.reactome.v2026.1.Mm.symbols.gmt') == 'REACTOME'

    def test_fallback_uppercases_stem(self):
        # No known substring -> uppercase stem, truncated to 20 chars.
        assert _library_short_name('path/mystery_lib.gmt') == 'MYSTERY_LIB'


# ============================================================
# Section 2: prepare_gene_symbols
# ============================================================

class TestPrepareGeneSymbols:

    def _da(self):
        return pd.DataFrame({
            'contrast':         ['KO_vs_WT', 'KO_vs_WT', 'KO_vs_WT', 'KO_vs_WT', 'X_vs_Y'],
            'gene_symbol':      ['GeneA', 'GeneA', 'GeneB', None, 'GeneC'],
            'deqms_adj_pvalue': [0.01, 0.20, 0.05, 0.04, 0.03],
            'limma_adj_pvalue': [0.02, 0.30, 0.06, 0.05, 0.04],
        })

    def test_filters_by_contrast_and_drops_unmapped(self):
        df, stats = prepare_gene_symbols(self._da(), 'KO_vs_WT')
        assert set(df['gene_symbol']) == {'GeneA', 'GeneB'}   # X_vs_Y + null excluded
        assert stats['n_unmapped'] == 1
        assert stats['n_unique'] == 2

    def test_dedup_keeps_most_significant_per_gene(self):
        df, stats = prepare_gene_symbols(self._da(), 'KO_vs_WT')
        # GeneA appears twice (0.01, 0.20); the 0.01 row is kept.
        gene_a = df[df['gene_symbol'] == 'GeneA']
        assert len(gene_a) == 1
        assert gene_a['deqms_adj_pvalue'].iloc[0] == 0.01
        assert stats['n_collapsed'] == 1

    def test_uses_deqms_pval_when_present(self):
        _df, stats = prepare_gene_symbols(self._da(), 'KO_vs_WT')
        assert stats['pval_col'] == 'deqms_adj_pvalue'

    def test_falls_back_to_limma_when_deqms_all_null(self):
        da = self._da()
        da['deqms_adj_pvalue'] = np.nan
        _df, stats = prepare_gene_symbols(da, 'KO_vs_WT')
        assert stats['pval_col'] == 'limma_adj_pvalue'


# ============================================================
# Section 3: build_ranked_series (GSEA ranking metric)
# ============================================================

class TestBuildRankedSeries:

    def _df(self):
        return pd.DataFrame({
            'gene_symbol':   ['GeneA', 'GeneB'],
            'log2_fc':       [2.0, -1.0],
            'deqms_t':       [5.0, -3.0],
            'limma_t':       [4.0, -2.0],
            'deqms_pvalue':  [1e-3, 1e-2],
            'limma_pvalue':  [1e-2, 1e-1],
        })

    def test_t_statistic_prefers_deqms_t(self):
        rnk = build_ranked_series(self._df(), 't_statistic')
        assert rnk['GeneA'] == 5.0
        assert rnk['GeneB'] == -3.0

    def test_t_statistic_falls_back_to_limma_t(self):
        df = self._df()
        df['deqms_t'] = np.nan
        rnk = build_ranked_series(df, 't_statistic')
        assert rnk['GeneA'] == 4.0

    def test_log2fc_ranking(self):
        rnk = build_ranked_series(self._df(), 'log2fc')
        assert rnk['GeneA'] == 2.0
        assert rnk['GeneB'] == -1.0

    def test_signed_log10p_sign_follows_fold_change(self):
        # GeneA: +fc, p=1e-3 -> +3;  GeneB: -fc, p=1e-2 -> -2
        rnk = build_ranked_series(self._df(), 'signed_log10p')
        assert rnk['GeneA'] == pytest.approx(3.0)
        assert rnk['GeneB'] == pytest.approx(-2.0)

    def test_unknown_ranking_exits(self):
        with pytest.raises(SystemExit):
            build_ranked_series(self._df(), 'bogus')


# ============================================================
# Section 4: build_protein_term_mapping
# ============================================================

class TestBuildProteinTermMapping:

    def _da(self):
        return pd.DataFrame({
            'protein_id':  ['P1', 'P2', 'P3'],
            'gene_symbol': ['GeneA', 'GeneB', 'GeneC'],
            'contrast':    ['KO_vs_WT', 'KO_vs_WT', 'KO_vs_WT'],
            'significant': [True, False, True],
        })

    def test_empty_enrichment_returns_empty_schema(self):
        out = build_protein_term_mapping(
            self._da(), [], [], pd.DataFrame(columns=['term_id']), {},
        )
        assert out.empty
        assert 'in_significant_set' in out.columns
        assert 'is_leading_edge' in out.columns

    def test_maps_annotated_proteins_to_tested_terms(self, tmp_path):
        gmt = tmp_path / 'lib.gmt'
        # Only TERM1 is tested; TERM2 must be ignored.
        gmt.write_text('TERM1\tdesc\tGeneA\tGeneB\nTERM2\tdesc\tGeneZ\n')
        enr = pd.DataFrame({'term_id': ['TERM1']})

        out = build_protein_term_mapping(
            self._da(), [str(gmt)], ['GO_BP'], enr, {},
        ).set_index('gene_symbol')

        assert set(out.index) == {'GeneA', 'GeneB'}          # GeneC/TERM2 excluded
        assert out.loc['GeneA', 'protein_id'] == 'P1'
        assert out.loc['GeneA', 'term_id'] == 'TERM1'
        assert out.loc['GeneA', 'library'] == 'GO_BP'
        # GeneA is significant, GeneB is not.
        assert bool(out.loc['GeneA', 'in_significant_set']) is True
        assert bool(out.loc['GeneB', 'in_significant_set']) is False
        # No GSEA supplied -> not leading edge.
        assert bool(out.loc['GeneA', 'is_leading_edge']) is False


# ============================================================
# Section 5: string helpers
# ============================================================

class TestStringHelpers:

    def test_truncate_label_short_unchanged(self):
        assert _truncate_label('short label') == 'short label'

    def test_truncate_label_long_truncated(self):
        s = 'x' * 80
        out = _truncate_label(s, maxlen=55)
        assert len(out) == 55
        assert out.endswith('...')

    def test_msigdb_gobp_to_phrase(self):
        assert _msigdb_name_to_go_lookup_phrase('GOBP_APOPTOTIC_PROCESS') == 'apoptotic process'

    def test_msigdb_gomf_to_phrase(self):
        assert _msigdb_name_to_go_lookup_phrase('GOMF_DNA_BINDING') == 'dna binding'

    def test_msigdb_no_go_prefix(self):
        # Non-GO term: no prefix stripped, underscores -> spaces, lowercased.
        assert _msigdb_name_to_go_lookup_phrase('REACTOME_SIGNALING') == 'reactome signaling'


# ============================================================
# Section 6: cluster_go_terms -- graceful degradation (no rpy2)
# ============================================================
# The real rrvgo clustering (needs R) is not run here. These cover the fallbacks:
# the three cluster columns are added and left null, and input rows preserved.
# Empty / non-GO inputs return before the rpy2 import; the GO-library case blocks
# rpy2 so the lazy `import rpy2.robjects` fails cleanly instead of starting R.

class TestClusterGoTermsDegradation:

    def test_empty_input_gets_null_columns(self):
        out = cluster_go_terms(pd.DataFrame(columns=['library', 'term_id']))
        for col in ('cluster_id', 'is_representative', 'parent_term'):
            assert col in out.columns
        assert out.empty

    def test_non_go_library_skips_with_null_columns(self):
        enr = pd.DataFrame({'library': ['REACTOME', 'REACTOME'],
                            'term_id': ['R1', 'R2'], 'adj_pvalue': [0.01, 0.02]})
        out = cluster_go_terms(enr)
        assert len(out) == 2                                  # rows preserved
        assert out['cluster_id'].isna().all()
        assert out['parent_term'].isna().all()

    def test_go_library_without_rpy2_returns_null_columns(self, monkeypatch):
        # Block rpy2 for this call so cluster_go_terms' lazy `import rpy2.robjects`
        # raises ImportError (instead of starting embedded R, which would segfault
        # where R is not linked), taking the graceful-degradation path.
        monkeypatch.setitem(sys.modules, 'rpy2', None)
        enr = pd.DataFrame({'library': ['GO_BP', 'GO_BP'],
                            'term_id': ['GOBP_A', 'GOBP_B'], 'adj_pvalue': [0.01, 0.02]})
        out = cluster_go_terms(enr)
        assert len(out) == 2                                  # input preserved
        assert out['cluster_id'].isna().all()
        assert out['is_representative'].isna().all()
        assert out['parent_term'].isna().all()

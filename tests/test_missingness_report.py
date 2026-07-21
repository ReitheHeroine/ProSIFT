#!/usr/bin/env python3
# title: test_missingness_report.py
# project: ProSIFT (PROtein Statistical Integration and Filtering Tool)
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-07-21
# last modified: 2026-07-21
#
# purpose:
#   Unit tests for Module 01 Process 4.10 missingness report
#   (bin/missingness_report.py). This module is mostly Plotly plotting, but the
#   filter-category axis and the filtered-out heatmap decide which proteins are
#   surfaced as removed, so their category bookkeeping is worth pinning.
#
#   filter_proteins.py now emits a WEAK-ANCHOR filter_status: proteins that
#   would be PARTIAL/SINGLE-GROUP but lack a fully-detected anchor group.
#   WEAK-ANCHOR is a REMOVED category (like SPARSE/ABSENT). These tests confirm
#   the module (a) keeps WEAK-ANCHOR on the filter-category bar chart with the
#   right count, (b) includes WEAK-ANCHOR proteins in the removed set that feeds
#   the filtered-out missingness heatmap, and (c) has WEAK-ANCHOR entries in the
#   module-level category order and color maps so rendering does not KeyError.
#
#   Tests are hermetic: small in-memory DataFrames, no file I/O, no network,
#   no PNG rendering (kaleido). We assert on returned Plotly figure data (bar
#   trace x/y, heatmap y labels) rather than on rendered HTML strings.
#
# inputs:
#   None (tests build inputs in-memory).
#
# outputs:
#   Test results (stdout via pytest).
#
# usage example:
#   pytest tests/test_missingness_report.py -v
#
#   copy/paste: pytest tests/test_missingness_report.py -v

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# Add bin/ to path so we can import missingness_report directly.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'bin'))

import missingness_report as mr
from missingness_report import (
    _CATEGORY_COLORS,
    _CATEGORY_ORDER,
    plot_filter_categories,
    plot_missingness_heatmap,
)

# ============================================================
# Section 0: shared in-memory inputs
# ============================================================
# A filter table exercising every category, with two WEAK-ANCHOR proteins so
# the count is distinguishable from a trivial "1 each" pass.
#   1 PASSED + 1 PARTIAL + 1 SINGLE-GROUP + 2 WEAK-ANCHOR + 1 SPARSE + 1 ABSENT
_STATUSES = [
    'PASSED', 'PARTIAL', 'SINGLE-GROUP',
    'WEAK-ANCHOR', 'WEAK-ANCHOR',
    'SPARSE', 'ABSENT',
]
_PROTEIN_IDS = [f'P{i}' for i in range(len(_STATUSES))]


def _filter_df():
    '''One row per protein: protein_id + filter_status.'''
    return pd.DataFrame({
        'protein_id': _PROTEIN_IDS,
        'filter_status': _STATUSES,
    })


def _matrix_df():
    '''Minimal validated_matrix shape: protein_id + one column per sample.
    Bare sample IDs (no abundance_ prefix), 2 WT + 2 KO. A few NaN so the
    heatmap has both detected and missing cells; exact values are irrelevant.'''
    rng = np.random.default_rng(seed=7)
    n = len(_PROTEIN_IDS)
    data = {'protein_id': _PROTEIN_IDS}
    for sid in ['WT-1', 'WT-2', 'KO-1', 'KO-2']:
        col = rng.normal(loc=20.0, scale=1.0, size=n)
        col[rng.integers(0, n)] = np.nan       # inject at least one missing
        data[sid] = col
    return pd.DataFrame(data)


def _metadata_df():
    '''Minimal validated_metadata: sample_id + group column.'''
    return pd.DataFrame({
        'sample_id': ['WT-1', 'WT-2', 'KO-1', 'KO-2'],
        'group': ['WT', 'WT', 'KO', 'KO'],
    })


# ============================================================
# Section 1: plot_filter_categories keeps WEAK-ANCHOR (check B.5)
# ============================================================

class TestFilterCategoryAxis:

    def test_weak_anchor_on_axis_with_count(self):
        '''The reindex to _CATEGORY_ORDER must retain WEAK-ANCHOR, not drop it.'''
        fig = plot_filter_categories(_filter_df(), 'run')
        bar = fig.data[0]
        x_cats = list(bar.x)
        y_vals = list(bar.y)
        assert 'WEAK-ANCHOR' in x_cats
        idx = x_cats.index('WEAK-ANCHOR')
        assert y_vals[idx] == 2                     # two WEAK-ANCHOR proteins

    def test_all_categories_present_and_sum_to_total(self):
        '''No category dropped: the six canonical categories all appear and
        their heights sum to the number of input proteins.'''
        fig = plot_filter_categories(_filter_df(), 'run')
        bar = fig.data[0]
        x_cats = list(bar.x)
        y_vals = [int(v) for v in bar.y]
        for cat in _CATEGORY_ORDER:
            assert cat in x_cats
        assert sum(y_vals) == len(_PROTEIN_IDS)


# ============================================================
# Section 2: WEAK-ANCHOR in the removed set feeding the heatmap (check B.6)
# ============================================================

class TestHeatmapRemovedSet:

    def test_weak_anchor_proteins_in_heatmap(self):
        '''WEAK-ANCHOR proteins (P3, P4) must appear as rows in the filtered-out
        heatmap; PASSED (P0) must not. This exercises the removed_categories set
        that selects rows for the heatmap.'''
        fig = plot_missingness_heatmap(
            _filter_df(), _matrix_df(), _metadata_df(), 'group', 'run'
        )
        heatmap = fig.data[0]
        y_ids = list(heatmap.y)
        assert 'P3' in y_ids and 'P4' in y_ids      # both WEAK-ANCHOR proteins
        assert 'P0' not in y_ids                     # PASSED is retained, not shown

    def test_heatmap_does_not_raise_on_weak_anchor(self):
        '''Sorting + category-boundary annotation look up each removed protein's
        category; a missing WEAK-ANCHOR rank/color would KeyError here.'''
        fig = plot_missingness_heatmap(
            _filter_df(), _matrix_df(), _metadata_df(), 'group', 'run'
        )
        assert fig is not None


# ============================================================
# Section 3: module-level category maps contain WEAK-ANCHOR (check B.7)
# ============================================================

class TestCategoryMaps:

    def test_color_map_has_weak_anchor(self):
        assert 'WEAK-ANCHOR' in _CATEGORY_COLORS

    def test_category_order_has_weak_anchor(self):
        '''_CATEGORY_ORDER is the module's category-rank map (drives reindex and
        the removed-category ordering). A WEAK-ANCHOR entry keeps rendering from
        raising KeyError on a WEAK-ANCHOR protein.'''
        assert 'WEAK-ANCHOR' in _CATEGORY_ORDER

#!/usr/bin/env python3
# title: test_qc_report_assembly.py
# project: ProSIFT (PROtein Statistical Integration and Filtering Tool)
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-07-09
# last modified: 2026-07-21
#
# purpose:
#   Unit tests for Module 03b QC Report Assembly (bin/qc_report_assembly.py).
#   Module 03b is a report-assembly module: most of it is Plotly plotting and
#   HTML rendering (not unit-tested), but a few pure-logic functions decide the
#   numbers that land in the report and are worth covering:
#     - extract_abundance      : abundance parsing + method-aware log2
#     - build_run_overview      : filter-category counts, flag bucketing, CV medians
#     - build_table_html        : DataFrame -> HTML with row truncation
#     - make_link               : path -> link or "not available"
#
#   Section 5 (added 2026-07-21) pins WEAK-ANCHOR handling. filter_proteins.py
#   now emits a WEAK-ANCHOR filter_status for proteins that would be PARTIAL or
#   SINGLE-GROUP but lack a fully-detected anchor group. WEAK-ANCHOR is a REMOVED
#   category (like SPARSE/ABSENT); only PASSED/PARTIAL/SINGLE-GROUP are RETAINED.
#   These tests confirm the overview counts and the category bar chart both know
#   about WEAK-ANCHOR (listed, counted as removed, not dropped, has a color).
#
#   Not tested here: the other plot_* Plotly functions, fig_to_div, load_params,
#   parse_args, main (I/O / rendering).
#
#   Note: extract_abundance defaults abundance_type to 'raw' via .get(), which
#   is the inconsistency flagged in the 2026-07-08 zero-handling review (now
#   that VALIDATE_INPUTS requires abundance_type). test_defaults_to_raw
#   documents the CURRENT behavior; see handoff Section 3.1 follow-up.
#
# inputs:
#   None (tests build inputs in-memory / in tmp_path).
#
# outputs:
#   Test results (stdout via pytest).
#
# usage example:
#   pytest tests/test_qc_report_assembly.py -v
#
#   copy/paste: pytest tests/test_qc_report_assembly.py -v

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# Add bin/ to path so we can import qc_report_assembly directly.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'bin'))

from qc_report_assembly import (
    build_run_overview,
    build_table_html,
    extract_abundance,
    make_link,
    plot_filter_categories,
)

# ============================================================
# Section 1: extract_abundance
# ============================================================

class TestExtractAbundance:

    def _matrix(self):
        return pd.DataFrame({
            'protein_id': ['P1', 'P2'],
            'abundance_S1': [4.0, 8.0],
            'abundance_S2': [16.0, 2.0],
            'peptide_count_S1': [3, 5],
            'peptide_count_S2': [4, 6],
        })

    def _params(self, abundance_type='log2'):
        return {'input': {'abundance_prefix': 'abundance_',
                          'peptide_count_prefix': 'peptide_count_',
                          'abundance_type': abundance_type}}

    def test_strips_prefix_and_excludes_peptides(self):
        df, sample_ids = extract_abundance(self._matrix(), self._params('log2'))
        assert sample_ids == ['S1', 'S2']
        assert list(df.columns) == ['S1', 'S2']
        assert df.index.tolist() == ['P1', 'P2']       # protein_id is the index
        assert df.loc['P1', 'S1'] == 4.0               # log2 -> passthrough

    def test_raw_converted_to_log2(self):
        df, _ = extract_abundance(self._matrix(), self._params('raw'))
        # log2([4, 8]) = [2, 3]; log2([16, 2]) = [4, 1]
        assert df.loc['P1', 'S1'] == 2.0
        assert df.loc['P2', 'S1'] == 3.0
        assert df.loc['P1', 'S2'] == 4.0

    def test_raw_zero_becomes_nan(self):
        matrix = self._matrix()
        matrix.loc[0, 'abundance_S1'] = 0.0            # a zero on raw scale
        df, _ = extract_abundance(matrix, self._params('raw'))
        assert pd.isna(df.loc['P1', 'S1'])             # 0 -> NaN before log2

    def test_raw_negative_becomes_nan(self):
        '''Non-positive handling covers negatives too (<= 0), not just zero, so
        log2 never sees a negative (which would be NaN + a RuntimeWarning).'''
        matrix = self._matrix()
        matrix.loc[0, 'abundance_S1'] = -5.0           # a negative on raw scale
        df, _ = extract_abundance(matrix, self._params('raw'))
        assert pd.isna(df.loc['P1', 'S1'])             # <= 0 -> NaN before log2

    def test_defaults_to_raw_when_abundance_type_absent(self):
        '''
        DOCUMENTS the flagged inconsistency: abundance_type defaults to 'raw'
        (so log2 is applied) when the key is absent. VALIDATE_INPUTS now requires
        the key, so this default is unreachable in a normal pipeline run, but it
        is still a silent default here. See handoff Section 3.1 follow-up.
        '''
        params = {'input': {'abundance_prefix': 'abundance_',
                            'peptide_count_prefix': 'peptide_count_'}}  # no abundance_type
        df, _ = extract_abundance(self._matrix(), params)
        assert df.loc['P1', 'S1'] == 2.0               # log2(4) -> default was 'raw'


# ============================================================
# Section 2: build_run_overview (aggregation logic)
# ============================================================
# Returns an HTML string; we assert the computed numbers appear in it.

class TestBuildRunOverview:

    def _args(self, filter_status, n_flags, cv_wt=(0.1, 0.2, 0.3)):
        metadata = pd.DataFrame({
            'sample_id': ['WT-1', 'WT-2', 'WT-3', 'KO-1', 'KO-2', 'KO-3'],
            'genotype': ['WT', 'WT', 'WT', 'KO', 'KO', 'KO'],
        })
        filter_df = pd.DataFrame({
            'protein_id': [f'P{i}' for i in range(len(filter_status))],
            'filter_status': filter_status,
        })
        sample_flags = pd.DataFrame({
            'sample_id': ['WT-1', 'WT-2', 'WT-3', 'KO-1', 'KO-2', 'KO-3'],
            'n_flags': n_flags,
        })
        cv_summary = pd.DataFrame({
            'protein_id': ['P0', 'P1', 'P2'],
            'cv_WT': list(cv_wt),
            'cv_KO': [0.4, 0.5, 0.6],
        })
        imp_summary = pd.DataFrame({
            'protein_id': ['P0', 'P1'],
            'n_mnar_imputed': [2, 0],
            'n_mar_imputed': [0, 3],
            'n_imputed_total': [2, 3],
        })
        params = {
            'design': {'group_column': 'genotype'},
            'normalization': {'method': 'median'},
            'imputation': {'mode': 'mixed', 'mnar_method': 'minprob', 'mar_method': 'knn'},
            'qc': {'min_detections_per_group': 2},
        }
        return dict(run_id='run', metadata_df=metadata, filter_df=filter_df,
                    sample_flags_df=sample_flags, cv_summary_df=cv_summary,
                    imp_summary_df=imp_summary, params=params)

    def test_retained_count_and_category_percentages(self):
        # 3 PASSED + 1 PARTIAL + 1 SINGLE-GROUP + 1 SPARSE + 1 ABSENT = 7; retained 5.
        html = build_run_overview(**self._args(
            ['PASSED', 'PASSED', 'PASSED', 'PARTIAL', 'SINGLE-GROUP', 'SPARSE', 'ABSENT'],
            [0, 0, 0, 0, 0, 0],
        ))
        assert 'PASSED: 3' in html
        assert 'SPARSE: 1' in html
        assert '5 / 7' in html                         # retained / total

    def test_flag_bucketing(self):
        # n_flags [0,0,1,2,3,4] -> no=2, mild(1-2)=2, severe(3-4)=2
        html = build_run_overview(**self._args(
            ['PASSED'], [0, 0, 1, 2, 3, 4],
        ))
        assert '0 flags: 2, 1-2 flags: 2, 3-4 flags: 2' in html

    def test_severe_samples_surfaced(self):
        # KO-1 and KO-2 have >= 3 flags -> listed in the warning line.
        html = build_run_overview(**self._args(
            ['PASSED'], [0, 0, 0, 3, 4, 0],
        ))
        assert '3+ flags' in html
        assert 'KO-1' in html and 'KO-2' in html

    def test_median_cv_per_group(self):
        # cv_WT median of [0.1, 0.2, 0.3] = 0.200
        html = build_run_overview(**self._args(['PASSED'], [0] * 6, cv_wt=(0.1, 0.2, 0.3)))
        assert 'WT: 0.200' in html

    def test_empty_filter_table_no_division_error(self):
        '''total_input = 0 must not raise ZeroDivisionError (percentage guard).'''
        html = build_run_overview(**self._args([], [0] * 6))
        assert isinstance(html, str)
        assert '0 / 0' in html                          # retained / total


# ============================================================
# Section 3: build_table_html (row truncation)
# ============================================================

class TestBuildTableHtml:

    def test_no_truncation_shows_all_rows(self):
        df = pd.DataFrame({'a': [1, 2, 3], 'b': ['x', 'y', 'z']})
        html = build_table_html(df, max_rows=50)
        assert 'Showing first' not in html
        assert html.count('<tr>') == 4                  # 1 header + 3 body

    def test_truncation_caps_rows_and_adds_note(self):
        df = pd.DataFrame({'a': list(range(60))})
        html = build_table_html(df, max_rows=50)
        assert 'Showing first 50 of 60 rows' in html
        assert html.count('<tr>') == 51                 # 1 header + 50 body

    def test_float_values_formatted_to_four_places(self):
        df = pd.DataFrame({'x': [1.234567]})
        html = build_table_html(df)
        assert '1.2346' in html

    def test_empty_dataframe_does_not_crash(self):
        html = build_table_html(pd.DataFrame({'a': []}))
        assert '<table' in html
        assert 'Showing first' not in html


# ============================================================
# Section 4: make_link (path -> link or "not available")
# ============================================================

class TestMakeLink:

    def test_none_path_is_unavailable(self):
        out = make_link(None, 'Report')
        assert 'not available' in out
        assert 'Report' in out

    def test_nonexistent_path_is_unavailable(self, tmp_path):
        out = make_link(str(tmp_path / 'missing.txt'), 'Report')
        assert 'not available' in out

    def test_existing_path_becomes_anchor(self, tmp_path):
        f = tmp_path / 'report.html'
        f.write_text('x')
        out = make_link(str(f), 'Report')
        assert out == '<a href="report.html">Report</a>'   # links by basename


# ============================================================
# Section 5: WEAK-ANCHOR handling (approved 2026-07-21)
# ============================================================
# filter_proteins.py emits WEAK-ANCHOR as a REMOVED category. These tests pin
# that both the overview counts and the filter-category bar chart know about it.

# A mixed filter table exercising every category. Deliberately > 1 of some so
# counts are distinguishable from a trivial "1 each" pass.
#   3 PASSED + 1 PARTIAL + 1 SINGLE-GROUP + 2 WEAK-ANCHOR + 1 SPARSE + 1 ABSENT
#   = 9 total; RETAINED (PASSED+PARTIAL+SINGLE-GROUP) = 5.
_WEAK_ANCHOR_STATUSES = [
    'PASSED', 'PASSED', 'PASSED',
    'PARTIAL', 'SINGLE-GROUP',
    'WEAK-ANCHOR', 'WEAK-ANCHOR',
    'SPARSE', 'ABSENT',
]


class TestWeakAnchorHandling:

    def _overview_args(self, filter_status):
        '''Same shape as TestBuildRunOverview._args, kept local so this section
        is self-contained. n_flags/cv/imputation values are irrelevant here.'''
        metadata = pd.DataFrame({
            'sample_id': ['WT-1', 'WT-2', 'WT-3', 'KO-1', 'KO-2', 'KO-3'],
            'genotype': ['WT', 'WT', 'WT', 'KO', 'KO', 'KO'],
        })
        filter_df = pd.DataFrame({
            'protein_id': [f'P{i}' for i in range(len(filter_status))],
            'filter_status': filter_status,
        })
        sample_flags = pd.DataFrame({
            'sample_id': ['WT-1', 'WT-2', 'WT-3', 'KO-1', 'KO-2', 'KO-3'],
            'n_flags': [0] * 6,
        })
        cv_summary = pd.DataFrame({
            'protein_id': ['P0', 'P1', 'P2'],
            'cv_WT': [0.1, 0.2, 0.3],
            'cv_KO': [0.4, 0.5, 0.6],
        })
        imp_summary = pd.DataFrame({
            'protein_id': ['P0', 'P1'],
            'n_mnar_imputed': [2, 0],
            'n_mar_imputed': [0, 3],
            'n_imputed_total': [2, 3],
        })
        params = {
            'design': {'group_column': 'genotype'},
            'normalization': {'method': 'median'},
            'imputation': {'mode': 'mixed', 'mnar_method': 'minprob', 'mar_method': 'knn'},
            'qc': {'min_detections_per_group': 2},
        }
        return dict(run_id='run', metadata_df=metadata, filter_df=filter_df,
                    sample_flags_df=sample_flags, cv_summary_df=cv_summary,
                    imp_summary_df=imp_summary, params=params)

    # --- Check A.1: overview lists WEAK-ANCHOR with its correct count ---
    def test_overview_lists_weak_anchor_count(self):
        html = build_run_overview(**self._overview_args(_WEAK_ANCHOR_STATUSES))
        assert 'WEAK-ANCHOR: 2' in html               # not silently omitted

    # --- Check A.2: WEAK-ANCHOR counted as REMOVED, not retained ---
    def test_weak_anchor_excluded_from_retained(self):
        html = build_run_overview(**self._overview_args(_WEAK_ANCHOR_STATUSES))
        # RETAINED = PASSED(3) + PARTIAL(1) + SINGLE-GROUP(1) = 5; total = 9.
        # If WEAK-ANCHOR were (wrongly) retained, this would read 7 / 9.
        assert '5 / 9' in html

    # --- Check A.3: per-category counts shown sum to total input ---
    def test_category_counts_sum_to_total(self):
        args = self._overview_args(_WEAK_ANCHOR_STATUSES)
        html = build_run_overview(**args)
        total_input = len(args['filter_df'])
        # Parse each "<li>CATEGORY: N (pct%)</li>" row. The Retained summary row
        # is "<li><b>Retained..." and does not match this pattern, so it is
        # excluded from the per-category sum.
        pairs = re.findall(r'<li>([A-Z-]+): ([\d,]+) \(', html)
        counts = {cat: int(n.replace(',', '')) for cat, n in pairs}
        assert sum(counts.values()) == total_input     # no category dropped
        assert counts['WEAK-ANCHOR'] == 2

    # --- Check A.4: bar chart has a WEAK-ANCHOR bar of correct height ---
    def test_plot_filter_categories_has_weak_anchor_bar(self):
        filter_df = pd.DataFrame({
            'protein_id': [f'P{i}' for i in range(len(_WEAK_ANCHOR_STATUSES))],
            'filter_status': _WEAK_ANCHOR_STATUSES,
        })
        # Must not raise on the WEAK-ANCHOR color lookup.
        fig = plot_filter_categories(filter_df, 'run')
        bar = fig.data[0]
        x_cats = list(bar.x)
        y_vals = list(bar.y)
        assert 'WEAK-ANCHOR' in x_cats                 # present on the axis
        idx = x_cats.index('WEAK-ANCHOR')
        assert y_vals[idx] == 2                         # correct height
        # A color was assigned for every bar (no missing marker color).
        assert bar.marker.color is not None
        assert len(list(bar.marker.color)) == len(x_cats)

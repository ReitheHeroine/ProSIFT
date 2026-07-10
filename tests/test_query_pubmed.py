#!/usr/bin/env python3
# title: test_query_pubmed.py
# project: ProSIFT (PROtein Statistical Integration and Filtering Tool)
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-07-09
# last modified: 2026-07-09
#
# purpose:
#   Unit tests for Module 06 QUERY_PUBMED (bin/query_pubmed.py). Scope: the PMI
#   scoring math (compute_pmi), the pure-Python core that turns co-occurrence
#   counts into the literature-association score. The ESearch API calls
#   (PubMedClient) are network/key-dependent and out of scope here (they belong
#   in a requires_network integration test); the caching layer is tested
#   separately in test_prosift_cache.py.
#
# inputs:
#   None (tests build inputs in-memory).
#
# outputs:
#   Test results (stdout via pytest).
#
# usage example:
#   pytest tests/test_query_pubmed.py -v
#
#   copy/paste: pytest tests/test_query_pubmed.py -v

import math
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'bin'))

from query_pubmed import compute_pmi


class TestComputePMI:

    def test_ground_truth(self):
        '''
        Hand-computed: co=10, symbol=100, term=100, total=10000.
        p_co=1e-3, p_symbol=1e-2, p_term=1e-2. ratio = 1e-3/(1e-2*1e-2) = 10.
        PMI = log2(10) = 3.321928...
        '''
        pmi = compute_pmi(10, 100, 100, 10000)
        assert pmi == pytest.approx(math.log2(10.0))

    def test_independence_gives_zero(self):
        '''When co-occurrence equals chance (ratio = 1), PMI = log2(1) = 0.'''
        # co=1, symbol=100, term=100, total=10000 -> ratio = 1e-4/(1e-2*1e-2) = 1
        assert compute_pmi(1, 100, 100, 10000) == pytest.approx(0.0)

    def test_enrichment_is_positive(self):
        assert compute_pmi(50, 100, 100, 10000) > 0

    def test_depletion_is_negative(self):
        # co=1, symbol=1000, term=1000, total=10000 -> ratio = 1e-4/(0.1*0.1) = 0.01
        assert compute_pmi(1, 1000, 1000, 10000) < 0

    @pytest.mark.parametrize('co,sym,term,total', [
        (0, 100, 100, 10000),     # zero co-occurrence
        (10, 0, 100, 10000),      # zero symbol total
        (10, 100, 0, 10000),      # zero term count
        (10, 100, 100, 0),        # zero total articles
    ])
    def test_any_zero_count_returns_none(self, co, sym, term, total):
        '''PMI is undefined (log of 0) if any count is zero -> None.'''
        assert compute_pmi(co, sym, term, total) is None

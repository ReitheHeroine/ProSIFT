#!/usr/bin/env python3
# title: test_query_pubmed.py
# project: ProSIFT (PROtein Statistical Integration and Filtering Tool)
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-07-09
# last modified: 2026-07-14  (added build_term_query + cache_token coverage)
#
# purpose:
#   Unit tests for Module 06 QUERY_PUBMED (bin/query_pubmed.py). Scope: the PMI
#   scoring math (compute_pmi); the term-query builder (build_term_query -- plain
#   phrase vs compound boolean); and the filesystem-safe cache-key token
#   (cache_token). The ESearch API calls
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

from query_pubmed import build_term_query, cache_token, compute_pmi


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


class TestBuildTermQuery:
    '''build_term_query: plain exact-phrase vs compound boolean pass-through.'''

    def test_plain_term_is_tiab_phrase(self):
        assert build_term_query('ketamine') == '"ketamine"[tiab]'

    def test_multiword_plain_term_is_a_single_phrase(self):
        # No boolean operator -> the whole string is one exact tiab phrase.
        assert build_term_query('down syndrome') == '"down syndrome"[tiab]'

    def test_or_compound_passed_through_parenthesised(self):
        term = '"TBI"[tiab] OR "traumatic brain injury"[tiab]'
        assert build_term_query(term) == f'({term})'

    def test_and_compound_passed_through_parenthesised(self):
        assert build_term_query('a AND b') == '(a AND b)'

    def test_bare_or_compound_is_not_quoted_as_one_phrase(self):
        # Regression: quoting the whole compound (the pre-2026-07-14 co-occurrence
        # bug) produced a nonsense literal phrase that matched nothing.
        out = build_term_query('TBI OR traumatic brain injury')
        assert out == '(TBI OR traumatic brain injury)'
        assert '"TBI OR traumatic brain injury"' not in out

    def test_operator_requires_surrounding_spaces(self):
        # ' OR '/' AND ' are operators only with spaces; a bare word is a phrase.
        assert build_term_query('corticosterone') == '"corticosterone"[tiab]'


class TestCacheToken:
    '''cache_token: filesystem-safe key that preserves existing warm caches.'''

    def test_plain_word_unchanged(self):
        # Critical: single-word terms already in the cache must be stable, or the
        # key changes and silently forces a full re-query.
        assert cache_token('ketamine') == 'ketamine'
        assert cache_token('TBI') == 'TBI'

    def test_spaces_become_underscores(self):
        assert cache_token('down syndrome') == 'down_syndrome'

    def test_apostrophe_and_operators_sanitised(self):
        assert cache_token("alzheimer's") == 'alzheimer_s'
        assert cache_token('TBI OR traumatic brain injury') == \
            'TBI_OR_traumatic_brain_injury'

    def test_no_unsafe_filename_chars_survive(self):
        tok = cache_token('"x"[tiab] OR "y":z/w')
        assert all(c not in tok for c in ('/', ':', '"', ' ', '[', ']'))

    def test_empty_or_all_unsafe_falls_back(self):
        assert cache_token('') == 'term'
        assert cache_token('  ??  ') == 'term'

    def test_variant_terms_get_distinct_tokens(self):
        # The two variants a user might add as separate terms must not collide.
        assert cache_token('TBI') != cache_token('traumatic brain injury')

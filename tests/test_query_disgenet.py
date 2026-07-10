#!/usr/bin/env python3
# title: test_query_disgenet.py
# project: ProSIFT (PROtein Statistical Integration and Filtering Tool)
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-07-09
# last modified: 2026-07-09
#
# purpose:
#   Unit tests for Module 06 QUERY_DISGENET (bin/query_disgenet.py). Scope:
#   parse_disgenet_associations, which maps the DisGeNET v1 camelCase API
#   response into ProSIFT's association schema (a mis-parse silently corrupts
#   disease annotations). The REST client (DisGeNETClient, requires an API key)
#   and caching are out of scope here.
#
# inputs:
#   None (tests build raw-response fixtures in-memory).
#
# outputs:
#   Test results (stdout via pytest).
#
# usage example:
#   pytest tests/test_query_disgenet.py -v
#
#   copy/paste: pytest tests/test_query_disgenet.py -v

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'bin'))

from query_disgenet import parse_disgenet_associations


class TestParseDisgenetAssociations:

    def _raw(self):
        return [{
            'diseaseUMLSCUI': 'C0001',
            'diseaseName': 'Test Disease',
            'diseaseType': '[disease]',
            'score': 0.8,
            'ei': 0.9,
            'numPMIDs': 5,
            'el': 'Definitive',
        }]

    def test_maps_all_fields(self):
        out = parse_disgenet_associations(self._raw())
        assert len(out) == 1
        row = out[0]
        assert row['disease_id'] == 'C0001'
        assert row['disease_name'] == 'Test Disease'
        assert row['disease_type'] == '[disease]'
        assert row['gda_score'] == 0.8
        assert row['evidence_index'] == 0.9
        assert row['n_publications'] == 5
        assert row['source'] == 'Definitive'

    def test_numeric_coercion_from_strings(self):
        '''Scores/counts arriving as strings are coerced to float/int.'''
        raw = self._raw()
        raw[0]['score'] = '0.8'
        raw[0]['numPMIDs'] = '5'
        row = parse_disgenet_associations(raw)[0]
        assert row['gda_score'] == 0.8 and isinstance(row['gda_score'], float)
        assert row['n_publications'] == 5 and isinstance(row['n_publications'], int)

    def test_missing_fields_become_none_and_zero(self):
        '''An empty association: optional fields -> None, n_publications -> 0.'''
        row = parse_disgenet_associations([{}])[0]
        assert row['disease_id'] is None
        assert row['gda_score'] is None
        assert row['evidence_index'] is None
        assert row['n_publications'] == 0        # numPMIDs defaults to 0
        assert row['source'] is None

    def test_empty_input(self):
        assert parse_disgenet_associations([]) == []

    def test_multiple_associations_preserved(self):
        raw = self._raw() + self._raw()
        raw[1]['diseaseUMLSCUI'] = 'C0002'
        out = parse_disgenet_associations(raw)
        assert [r['disease_id'] for r in out] == ['C0001', 'C0002']

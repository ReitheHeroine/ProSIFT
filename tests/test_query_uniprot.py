#!/usr/bin/env python3
# title: test_query_uniprot.py
# project: ProSIFT (PROtein Statistical Integration and Filtering Tool)
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-07-09
# last modified: 2026-07-09
#
# purpose:
#   Unit tests for Module 06 QUERY_UNIPROT (bin/query_uniprot.py). Scope:
#   _parse_uniprot_entry (extracts name / function / location / tissue /
#   keywords from a UniProt JSON entry) and build_accession_query. The paginated
#   REST retrieval (query_uniprot_batch, _request_with_retry) is network-
#   dependent and out of scope here.
#
# inputs:
#   None (tests build UniProt-entry fixtures in-memory).
#
# outputs:
#   Test results (stdout via pytest).
#
# usage example:
#   pytest tests/test_query_uniprot.py -v
#
#   copy/paste: pytest tests/test_query_uniprot.py -v

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'bin'))

from query_uniprot import _parse_uniprot_entry, build_accession_query


class TestBuildAccessionQuery:

    def test_single(self):
        assert build_accession_query(['Q9WTX5']) == '(accession:Q9WTX5)'

    def test_multiple_joined_with_or(self):
        assert build_accession_query(['Q1', 'P2']) == '(accession:Q1) OR (accession:P2)'


class TestParseUniprotEntry:

    def test_recommended_name(self):
        entry = {'proteinDescription': {'recommendedName': {'fullName': {'value': 'My Protein'}}}}
        assert _parse_uniprot_entry(entry)['protein_name'] == 'My Protein'

    def test_submission_name_fallback(self):
        '''When no recommendedName, fall back to the first submissionName.'''
        entry = {'proteinDescription': {'submissionNames': [{'fullName': {'value': 'Sub Name'}}]}}
        assert _parse_uniprot_entry(entry)['protein_name'] == 'Sub Name'

    def test_function_comment_joined(self):
        entry = {'comments': [{'commentType': 'FUNCTION',
                               'texts': [{'value': 'Does X'}, {'value': 'Does Y'}]}]}
        assert _parse_uniprot_entry(entry)['function_description'] == 'Does X; Does Y'

    def test_subcellular_location(self):
        entry = {'comments': [{
            'commentType': 'SUBCELLULAR LOCATION',
            'subcellularLocations': [
                {'location': {'value': 'Cytoplasm'}},
                {'location': {'value': 'Nucleus'}},
            ],
        }]}
        assert _parse_uniprot_entry(entry)['subcellular_location'] == 'Cytoplasm; Nucleus'

    def test_tissue_specificity(self):
        entry = {'comments': [{'commentType': 'TISSUE SPECIFICITY',
                               'texts': [{'value': 'Expressed in brain'}]}]}
        assert _parse_uniprot_entry(entry)['tissue_expression'] == 'Expressed in brain'

    def test_keywords_joined(self):
        entry = {'keywords': [{'name': 'Kinase'}, {'name': 'ATP-binding'}]}
        assert _parse_uniprot_entry(entry)['keywords'] == 'Kinase; ATP-binding'

    def test_empty_entry_all_none(self):
        out = _parse_uniprot_entry({})
        assert all(out[f] is None for f in (
            'protein_name', 'function_description', 'subcellular_location',
            'tissue_expression', 'keywords',
        ))

    def test_unrelated_comment_types_ignored(self):
        '''A comment type we do not extract leaves the fields null.'''
        entry = {'comments': [{'commentType': 'SIMILARITY',
                               'texts': [{'value': 'Belongs to family X'}]}]}
        out = _parse_uniprot_entry(entry)
        assert out['function_description'] is None
        assert out['subcellular_location'] is None

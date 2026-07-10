#!/usr/bin/env python3
# title: test_query_dgidb.py
# project: ProSIFT (PROtein Statistical Integration and Filtering Tool)
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-07-09
# last modified: 2026-07-09
#
# purpose:
#   Unit tests for Module 06 QUERY_DGIDB (bin/query_dgidb.py). Scope:
#   parse_dgidb_response, which flattens the DGIdb GraphQL response into
#   drug-interaction rows (approval-status mapping, first interaction type,
#   deduplicated sources/PMIDs). The GraphQL client (query_dgidb_gene) is
#   network-dependent and out of scope here.
#
# inputs:
#   None (tests build GraphQL-response fixtures in-memory).
#
# outputs:
#   Test results (stdout via pytest).
#
# usage example:
#   pytest tests/test_query_dgidb.py -v
#
#   copy/paste: pytest tests/test_query_dgidb.py -v

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'bin'))

from query_dgidb import parse_dgidb_response


def _response(interactions):
    '''Wrap a list of interaction dicts in the GraphQL response envelope.'''
    return {'data': {'genes': {'nodes': [{'interactions': interactions}]}}}


class TestParseDgidbResponse:

    def test_no_genes_returns_empty(self):
        assert parse_dgidb_response({'data': {'genes': {'nodes': []}}}, 'GENE') == []

    def test_gene_with_no_interactions_returns_empty(self):
        assert parse_dgidb_response(_response([]), 'GENE') == []

    def test_parses_interaction_fields(self):
        raw = _response([{
            'drug': {'name': 'Aspirin', 'conceptId': 'chembl:CHEMBL25', 'approved': True},
            'interactionTypes': [{'type': 'inhibitor'}],
            'interactionScore': 1.5,
            'interactionClaims': [
                {'source': {'fullName': 'DrugBank'},
                 'publications': [{'pmid': 111}, {'pmid': 222}]},
            ],
        }])
        row = parse_dgidb_response(raw, 'GENE')[0]
        assert row['drug_name'] == 'Aspirin'
        assert row['drug_concept_id'] == 'chembl:CHEMBL25'
        assert row['interaction_type'] == 'inhibitor'
        assert row['interaction_score'] == 1.5
        assert row['approval_status'] == 'approved'
        assert row['n_sources'] == 1
        assert row['sources'] == 'DrugBank'
        assert row['pmids'] == '111; 222'

    def test_approval_status_mapping(self):
        def status(approved):
            raw = _response([{'drug': {'name': 'D', 'approved': approved},
                              'interactionClaims': []}])
            return parse_dgidb_response(raw, 'G')[0]['approval_status']
        assert status(True) == 'approved'
        assert status(False) == 'not_approved'
        assert status(None) is None            # unknown -> None

    def test_first_interaction_type_used(self):
        raw = _response([{
            'drug': {'name': 'D', 'approved': True},
            'interactionTypes': [{'type': 'inhibitor'}, {'type': 'antagonist'}],
            'interactionClaims': [],
        }])
        assert parse_dgidb_response(raw, 'G')[0]['interaction_type'] == 'inhibitor'

    def test_sources_and_pmids_deduped_and_sorted(self):
        raw = _response([{
            'drug': {'name': 'D', 'approved': True},
            'interactionClaims': [
                {'source': {'fullName': 'TTD'},
                 'publications': [{'pmid': 222}, {'pmid': 111}]},
                {'source': {'fullName': 'DrugBank'},
                 'publications': [{'pmid': 222}]},          # duplicate pmid
            ],
        }])
        row = parse_dgidb_response(raw, 'G')[0]
        assert row['sources'] == 'DrugBank; TTD'            # sorted, deduped
        assert row['pmids'] == '111; 222'                    # sorted, deduped
        assert row['n_sources'] == 2

    def test_no_claims_gives_null_sources_pmids(self):
        raw = _response([{'drug': {'name': 'D', 'approved': True},
                          'interactionClaims': []}])
        row = parse_dgidb_response(raw, 'G')[0]
        assert row['sources'] is None
        assert row['pmids'] is None
        assert row['n_sources'] == 0

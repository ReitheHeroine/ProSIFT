#!/usr/bin/env python3
# title: test_query_ctd.py
# project: ProSIFT (PROtein Statistical Integration and Filtering Tool)
# author: Reina Hastings
# contact: reinahastings13@gmail.com
# date created: 2026-07-09
# last modified: 2026-07-09
#
# purpose:
#   Unit tests for Module 06 QUERY_CTD (bin/query_ctd.py). Scope:
#   build_entrez_lookups (mapping table -> Entrez lookups) and filter_ctd_file
#   (organism-aware filtering of the bulk CTD TSV). CTD is a local bulk-file
#   query (no live API), so filter_ctd_file is exercised end to end against a
#   small gzipped fixture -- the real behavior, not a mock.
#
# inputs:
#   None (tests build a mapping DataFrame and a gzipped CTD fixture in tmp_path).
#
# outputs:
#   Test results (stdout via pytest).
#
# usage example:
#   pytest tests/test_query_ctd.py -v
#
#   copy/paste: pytest tests/test_query_ctd.py -v

import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'bin'))

from query_ctd import build_entrez_lookups, filter_ctd_file

# CTD column order (bin/query_ctd.py CTD_COLUMNS): ChemicalName, ChemicalID,
# CasRN, GeneSymbol, GeneID, GeneForms, Organism, OrganismID, Interaction,
# InteractionActions, PubMedIDs.


def _write_ctd(path, rows):
    '''Write a gzipped, headerless, tab-separated CTD file (with a # comment).'''
    with gzip.open(path, 'wt') as fh:
        fh.write('# CTD_chem_gene_ixns fixture\n')
        for r in rows:
            fh.write('\t'.join(r) + '\n')


def _row(chem, mesh, cas, gene_sym, gene_id, organism_id, pmids):
    return [chem, mesh, cas, gene_sym, gene_id, 'protein',
            'organism', organism_id, f'{chem} interacts {gene_sym}',
            'increases^expression', pmids]


# ============================================================
# build_entrez_lookups
# ============================================================

class TestBuildEntrezLookups:

    def test_builds_mouse_and_human_lookups(self):
        mapping = pd.DataFrame({
            'protein_id':            ['P1', 'P2'],
            'entrez_id_mouse':       [12345, np.nan],
            'human_ortholog_entrez': [np.nan, 67890],
        })
        mouse, human, all_entrez = build_entrez_lookups(mapping)
        assert mouse == {'12345': ['P1']}
        assert human == {'67890': ['P2']}
        assert all_entrez == {'12345', '67890'}

    def test_missing_entrez_ids_skipped(self):
        mapping = pd.DataFrame({
            'protein_id':            ['P1'],
            'entrez_id_mouse':       [np.nan],
            'human_ortholog_entrez': [np.nan],
        })
        mouse, human, all_entrez = build_entrez_lookups(mapping)
        assert mouse == {} and human == {} and all_entrez == set()

    def test_multiple_proteins_same_entrez(self):
        mapping = pd.DataFrame({
            'protein_id':            ['P1', 'P2'],
            'entrez_id_mouse':       [12345, 12345],
            'human_ortholog_entrez': [np.nan, np.nan],
        })
        mouse, _human, _all = build_entrez_lookups(mapping)
        assert mouse == {'12345': ['P1', 'P2']}


# ============================================================
# filter_ctd_file
# ============================================================

class TestFilterCtdFile:

    def test_matches_mouse_by_organism_and_parses_fields(self, tmp_path):
        path = tmp_path / 'ctd.tsv.gz'
        _write_ctd(path, [
            _row('Aspirin', 'MESH:D001241', '50-78-2', 'GeneA', '12345', '10090', '111|222'),
            _row('Other', 'MESH:D000000', '', 'GeneZ', '99999', '9606', ''),  # not in lookup
        ])
        rows = filter_ctd_file(path, {'12345': ['P1']}, {}, {'12345'})

        assert len(rows) == 1                      # the 99999 row is filtered out
        r = rows[0]
        assert r['protein_id'] == 'P1'
        assert r['query_organism'] == 'mouse'
        assert r['gene_symbol_queried'] == 'GeneA'
        assert r['chemical_name'] == 'Aspirin'
        assert r['chemical_mesh_id'] == 'D001241'  # 'MESH:' prefix stripped
        assert r['chemical_cas_rn'] == '50-78-2'
        assert r['n_publications'] == 2            # '111|222' -> 2
        assert r['pmids'] == '111; 222'
        assert r['ctd_query_status'] == 'success'

    def test_empty_cas_and_pmids_become_none(self, tmp_path):
        path = tmp_path / 'ctd.tsv.gz'
        _write_ctd(path, [
            _row('Chem', 'MESH:D999', '', 'GeneA', '12345', '10090', ''),
        ])
        r = filter_ctd_file(path, {'12345': ['P1']}, {}, {'12345'})[0]
        assert r['chemical_cas_rn'] is None
        assert r['n_publications'] == 0
        assert r['pmids'] is None

    def test_no_matching_rows(self, tmp_path):
        path = tmp_path / 'ctd.tsv.gz'
        _write_ctd(path, [
            _row('Chem', 'MESH:D1', '', 'GeneZ', '99999', '10090', '1'),
        ])
        assert filter_ctd_file(path, {'12345': ['P1']}, {}, {'12345'}) == []

    def test_gene_matches_via_fallback_when_organism_mismatches(self, tmp_path):
        '''A mouse Entrez ID appearing under a non-mouse OrganismID still matches
        via the gene-id fallback (organism-agnostic second pass).'''
        path = tmp_path / 'ctd.tsv.gz'
        _write_ctd(path, [
            _row('Chem', 'MESH:D1', '', 'GeneA', '12345', '9606', '5'),  # mouse eid, human org
        ])
        rows = filter_ctd_file(path, {'12345': ['P1']}, {}, {'12345'})
        assert len(rows) == 1
        assert rows[0]['protein_id'] == 'P1'
        assert rows[0]['query_organism'] == 'mouse'

# =============================================================================
# Title:         helper-fixture.R
# Project:       ProSIFT (PROtein Statistical Integration and Filtering Tool)
# Author:        Reina Hastings (reinahastings13@gmail.com)
# Created:       2026-07-15
# Last Modified: 2026-07-15
# Purpose:       Builds a small synthetic Module 07 results database in-code for
#                the Module 08 frontend gate. Only the columns db.R actually
#                reads are populated. Values reuse the Python suite's synthetic
#                run (P1..P5, KO_vs_WT, GOBP_X/Y) for cross-language consistency,
#                and add the adversarial rows the invariants need (placeholder
#                query-status rows, a DEqMS-NULL protein, a 0.05 boundary term,
#                an imputed-but-flagged sample, a null-gene protein, and single-
#                analysis terms). No binary DB is committed; only this builder.
# Inputs:        None.
# Outputs:       build_fixture_db(path) writes a SQLite file at `path`;
#                fixture_path() returns a cached temp DB path; fixture_con()
#                opens a read-only connection to it (via db.R open_results_db).
# Usage:         con <- fixture_con(); on.exit(close_results_db(con))
# =============================================================================

# --- Table builders (mirror the Module 07 schema, db.R-read columns only) ----

.fx_proteins <- function() {
  # P6 has a NULL gene symbol (A9: null gene sorts last) and no rows in most
  # child tables (A11: LEFT JOIN keeps it with NA statistics).
  data.frame(
    protein_id             = c('P1', 'P2', 'P3', 'P4', 'P5', 'P6'),
    gene_symbol            = c('Gene1', 'Gene2', 'Gene3', 'Gene4', 'Gene5', NA),
    human_ortholog_symbol  = c('GENE1', 'GENE2', NA, 'GENE4', 'GENE5', 'GENE6'),
    human_ortholog_entrez  = c('111', '222', NA, '444', '555', '666'),
    ortholog_mapping_status = c('one_to_one', 'one_to_one', 'no_ortholog',
                                'one_to_one', 'one_to_one', 'one_to_one'),
    detection_category     = c('PASSED', 'PASSED', 'PASSED', 'SINGLE-GROUP',
                               'PARTIAL', 'PASSED'),
    imputation_fraction    = c(0.0, 0.5, 0.5, 0.5, 1.0, 0.0),
    stringsAsFactors = FALSE
  )
}

.fx_differential_abundance <- function() {
  # P2 has deqms_adj_pvalue NULL -> COALESCE falls back to limma (A3).
  # A second contrast (AKO_vs_WT, sorts before KO_vs_WT) exercises A15 ordering
  # and db_protein_contrasts. P6 has no row at all (A11).
  data.frame(
    protein_id       = c('P1', 'P2', 'P3', 'P4', 'P5', 'P1'),
    contrast         = c('KO_vs_WT', 'KO_vs_WT', 'KO_vs_WT', 'KO_vs_WT',
                         'KO_vs_WT', 'AKO_vs_WT'),
    log2_fc          = c(2.0, -0.1, 3.0, -2.5, 0.0, 1.0),
    deqms_adj_pvalue = c(0.008, NA, 0.004, 0.015, 0.99, 0.02),
    limma_adj_pvalue = c(0.01, 0.9, 0.005, 0.02, 0.99, 0.03),
    direction        = c('up', 'ns', 'up', 'down', 'ns', 'up'),
    significant      = c(1L, 0L, 1L, 1L, 0L, 1L),
    stringsAsFactors = FALSE
  )
}

.fx_uniprot_annotations <- function() {
  # P6 is intentionally absent (LEFT JOIN -> NA protein_name).
  data.frame(
    protein_id           = c('P1', 'P2', 'P3', 'P4', 'P5'),
    protein_name         = paste0('name', 1:5),
    function_description  = paste0('func', 1:5),
    subcellular_location = rep('Cytoplasm', 5),
    tissue_expression    = rep('brain', 5),
    keywords             = rep('kw', 5),
    uniprot_query_status = rep('success', 5),
    stringsAsFactors = FALSE
  )
}

.fx_disease_associations <- function() {
  # P1: two genuine records (A2 count, A12 order by gda_score desc).
  # P2: a no_associations placeholder -> must NOT be counted (A1).
  # P3: a no_ortholog placeholder -> must NOT be counted (A1).
  data.frame(
    protein_id            = c('P1', 'P1', 'P2', 'P3'),
    disgenet_query_status = c('success', 'success', 'no_associations', 'no_ortholog'),
    disease_name          = c('DiseaseA', 'DiseaseB', NA, NA),
    gda_score             = c(0.4, 0.6, NA, NA),
    disease_type          = c('[disease]', '[disease]', NA, NA),
    evidence_index        = c(1.0, 2.0, NA, NA),
    n_publications        = c(3L, 5L, NA, NA),
    stringsAsFactors = FALSE
  )
}

.fx_drug_interactions <- function() {
  # P1: two genuine records (A4 count, order by drug_name -> DrugA, DrugX).
  # P4: a no_interactions placeholder -> must NOT be counted (A1).
  data.frame(
    protein_id         = c('P1', 'P1', 'P4'),
    dgidb_query_status = c('success', 'success', 'no_interactions'),
    drug_name          = c('DrugX', 'DrugA', NA),
    interaction_type   = c('inhibitor', 'agonist', NA),
    approval_status    = c('approved', NA, NA),
    sources            = c('SourceA', 'SourceB', NA),
    stringsAsFactors = FALSE
  )
}

.fx_pubmed_cooccurrence <- function() {
  # P1: two genuine records; top_pmi = MAX(normalized_score) = 2.5 (A4).
  # P2: a no_hits placeholder carrying a HIGH score (9.9) that must be excluded
  #     from top_pmi (A1) -> a strong test that placeholders are not aggregated.
  data.frame(
    protein_id          = c('P1', 'P1', 'P2'),
    pubmed_query_status = c('success', 'success', 'no_hits'),
    search_term         = c('ketamine', 'TBI', 'ketamine'),
    mouse_hit_count     = c(1L, 0L, NA),
    human_hit_count     = c(2L, 3L, NA),
    mouse_symbol_used   = c('Gene1', 'Gene1', 'Gene2'),
    human_symbol_used   = c('GENE1', 'GENE1', 'GENE2'),
    normalized_score    = c(1.3, 2.5, 9.9),
    stringsAsFactors = FALSE
  )
}

.fx_chemical_interactions <- function() {
  # P1: two genuine records (order by n_publications desc -> ChemY, ChemX).
  # ChemY uses the multi-action CTD grammar for prettify_ctd_actions (B1).
  # P4: a no_interactions placeholder -> excluded by db_protein_chemicals.
  data.frame(
    protein_id          = c('P1', 'P1', 'P4'),
    ctd_query_status    = c('success', 'success', 'no_interactions'),
    chemical_name       = c('ChemX', 'ChemY', NA),
    interaction_actions = c('increases^expression',
                            'decreases^expression|increases^abundance', NA),
    interaction_text    = c('increases expression', 'decreases expression', NA),
    n_publications      = c(2L, 5L, NA),
    stringsAsFactors = FALSE
  )
}

.fx_sample_qc_flags <- function() {
  # S_WT-1 fires extreme_median; S_KO-1 fires pca_outlier (on different samples,
  # matching the Python fixture). Used with the imputation status below for A8.
  data.frame(
    sample_id            = c('S_WT-1', 'S_WT-2', 'S_KO-1', 'S_KO-2'),
    flag_low_detection   = c(0L, 0L, 0L, 0L),
    flag_extreme_median  = c(1L, 0L, 0L, 0L),
    flag_pca_outlier     = c(0L, 0L, 1L, 0L),
    flag_low_correlation = c(0L, 0L, 0L, 0L),
    n_flags              = c(1L, 0L, 1L, 0L),
    stringsAsFactors = FALSE
  )
}

.fx_sample_abundances <- function() {
  # P1 is OBSERVED in S_WT-1, S_WT-2, S_KO-2 but IMPUTED (mar) in S_KO-1.
  # S_KO-1 carries the pca_outlier flag; because P1 is imputed there, that flag
  # must NOT appear in db_protein_qc_flags(P1) (A8). "group" is the reserved
  # column name db.R reads as grp.
  data.frame(
    sample_id          = c('S_WT-1', 'S_WT-2', 'S_KO-1', 'S_KO-2'),
    group              = c('WT', 'WT', 'KO', 'KO'),
    protein_id         = rep('P1', 4),
    imputation_status  = c('observed', 'observed', 'mar', 'observed'),
    stringsAsFactors = FALSE
  )
}

.fx_enrichment_results <- function() {
  # GOBP_X: ORA sig + GSEA sig (A5 both flags 1).
  # GOBP_Y: ORA adj_pvalue EXACTLY 0.05 -> sig 0 (A6 strict boundary); GSEA ns.
  # GOBP_Z: ORA only -> gsea_sig COALESCE 0 (A7).
  # REACTOME_W: GSEA only -> ora_sig COALESCE 0 (A7); different library (C6).
  data.frame(
    term_id          = c('GOBP_X', 'GOBP_X', 'GOBP_Y', 'GOBP_Y', 'GOBP_Z', 'REACTOME_W'),
    term_name        = c('GOBP_X', 'GOBP_X', 'GOBP_Y', 'GOBP_Y', 'GOBP_Z', 'REACTOME_W'),
    library          = c('GO_BP', 'GO_BP', 'GO_BP', 'GO_BP', 'GO_BP', 'REACTOME'),
    gene_set_size    = c(20L, 20L, 30L, 30L, 15L, 40L),
    contrast         = rep('KO_vs_WT', 6),
    analysis_type    = c('ORA', 'GSEA', 'ORA', 'GSEA', 'ORA', 'GSEA'),
    adj_pvalue       = c(0.01, 0.02, 0.05, 0.20, 0.03, 0.04),
    odds_ratio       = c(3.0, NA, 2.0, NA, 4.0, NA),
    overlap_size     = c(3L, NA, 5L, NA, 4L, NA),
    enrichment_score = c(NA, 1.5, NA, -0.5, NA, -1.2),
    stringsAsFactors = FALSE
  )
}

.fx_protein_term_mapping <- function() {
  # P1 in GOBP_X across TWO contrasts (collapse to one row, is_leading_edge MAX
  # = 1 for A10) and in GOBP_Y but not-significant (excluded from A10).
  # GOBP_X / KO_vs_WT members: P1, P2, P3, and P6 (empty-string gene, ordered
  # last for A13). P2 carries the DEqMS-NULL DA row -> COALESCE fallback (A13).
  data.frame(
    protein_id         = c('P1', 'P1', 'P1', 'P2', 'P3', 'P6'),
    gene_symbol        = c('Gene1', 'Gene1', 'Gene1', 'Gene2', 'Gene3', ''),
    term_id            = c('GOBP_X', 'GOBP_X', 'GOBP_Y', 'GOBP_X', 'GOBP_X', 'GOBP_X'),
    term_name          = c('GOBP_X', 'GOBP_X', 'GOBP_Y', 'GOBP_X', 'GOBP_X', 'GOBP_X'),
    library            = rep('GO_BP', 6),
    in_significant_set = c(1L, 1L, 0L, 1L, 1L, 1L),
    is_leading_edge    = c(0L, 1L, 0L, 0L, 1L, 0L),
    contrast           = c('KO_vs_WT', 'AKO_vs_WT', 'KO_vs_WT', 'KO_vs_WT',
                           'KO_vs_WT', 'KO_vs_WT'),
    stringsAsFactors = FALSE
  )
}

.fx_run_metadata <- function() {
  data.frame(run_id = 'SYN_RUN', stringsAsFactors = FALSE)
}


# --- Assembly ---------------------------------------------------------------

#' Write the full synthetic fixture database to `path`.
build_fixture_db <- function(path) {
  con <- DBI::dbConnect(RSQLite::SQLite(), path)
  on.exit(DBI::dbDisconnect(con), add = TRUE)
  tables <- list(
    proteins                = .fx_proteins(),
    differential_abundance  = .fx_differential_abundance(),
    uniprot_annotations     = .fx_uniprot_annotations(),
    disease_associations    = .fx_disease_associations(),
    drug_interactions       = .fx_drug_interactions(),
    pubmed_cooccurrence     = .fx_pubmed_cooccurrence(),
    chemical_interactions   = .fx_chemical_interactions(),
    sample_qc_flags         = .fx_sample_qc_flags(),
    sample_abundances       = .fx_sample_abundances(),
    enrichment_results      = .fx_enrichment_results(),
    protein_term_mapping    = .fx_protein_term_mapping(),
    run_metadata            = .fx_run_metadata()
  )
  for (nm in names(tables)) {
    DBI::dbWriteTable(con, nm, tables[[nm]], overwrite = TRUE)
  }
  invisible(path)
}


# --- Cached read-only accessors for the unit / testServer suites ------------

.fx_cache <- new.env(parent = emptyenv())

#' Path to a lazily-built, process-cached fixture database (tests read only).
fixture_path <- function() {
  if (is.null(.fx_cache$path) || !file.exists(.fx_cache$path %||% '')) {
    p <- tempfile('prosift_fixture_', fileext = '.db')
    build_fixture_db(p)
    .fx_cache$path <- p
  }
  .fx_cache$path
}

#' Open a read-only connection to the fixture (exercises db.R open_results_db).
fixture_con <- function() open_results_db(fixture_path())

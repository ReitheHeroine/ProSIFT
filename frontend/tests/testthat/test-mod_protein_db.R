# =============================================================================
# Title:         test-mod_protein_db.R
# Project:       ProSIFT (PROtein Statistical Integration and Filtering Tool)
# Author:        Reina Hastings (reinahastings13@gmail.com)
# Created:       2026-07-15
# Last Modified: 2026-07-15
# Purpose:       Group C reactive-logic checks for the Protein Database View
#                module (frontend/R/mod_protein_db.R) using shiny::testServer.
#                Covers the significance-only filter (C3), click-through row ->
#                protein_id mapping plus the gene fallback (C4), and the
#                isolate() regression that stops a filter change from firing a
#                spurious click-through (C1).
# Inputs:        The synthetic fixture (helper-fixture.R); mod server (sourced).
# Outputs:       testthat assertions.
# Usage:         Rscript frontend/run_tests.R    (or testthat::test_file(...))
# =============================================================================

library(shiny)

# --- Shared fixture data (the db_protein_table for KO_vs_WT) -----------------

.pdb_table <- local({
  con <- fixture_con()
  on.exit(close_results_db(con))
  db_protein_table(con, 'KO_vs_WT')
})


# --- C1: the isolate regression (filter change must not fire a click) --------

test_that('C1: toggling the filter with a row selected does not invalidate the click reactive', {
  testServer(mod_protein_db_server,
             args = list(protein_data = reactive(.pdb_table)), {
    session$setInputs(sig_only = FALSE, table_rows_selected = 1L)

    fires <- 0L
    observe({ session$returned(); fires <<- fires + 1L })
    session$flushReact()
    base <- fires

    # A data change (the significance filter) with the SAME selection must not
    # invalidate the clicked-protein reactive -> no spurious navigation.
    session$setInputs(sig_only = TRUE)
    session$flushReact()
    expect_equal(fires, base)

    # A genuine new selection event DOES invalidate it.
    session$setInputs(table_rows_selected = 2L)
    session$flushReact()
    expect_gt(fires, base)
  })
})


# --- C3: significance-only filter --------------------------------------------

test_that('C3: sig_only filters to significant proteins and updates the count', {
  testServer(mod_protein_db_server,
             args = list(protein_data = reactive(.pdb_table)), {
    session$setInputs(sig_only = FALSE)
    expect_equal(nrow(filtered_data()), 6L)
    expect_equal(output$row_count, '6 proteins')

    # Significant in KO_vs_WT: P1, P3, P4.
    session$setInputs(sig_only = TRUE)
    expect_equal(nrow(filtered_data()), 3L)
    expect_true(all(filtered_data()$significant == 1L))
    expect_equal(output$row_count, '3 proteins')
  })
})


# --- C4: click-through maps to the right protein; gene fallback --------------

test_that('C4: clicked reactive returns the selected protein_id, NULL when none', {
  testServer(mod_protein_db_server,
             args = list(protein_data = reactive(.pdb_table)), {
    # No selection yet -> NULL.
    expect_null(session$returned())

    session$setInputs(sig_only = FALSE, table_rows_selected = 1L)
    # Rows are ordered gene-ascending, so row 1 is P1 (Gene1).
    expect_equal(session$returned(), 'P1')

    # Gene fallback: the null-gene protein P6 displays its accession, not blank.
    dd <- display_data()
    expect_equal(dd$Gene[dd$protein_id == 'P6'], 'P6')
  })
})

# =============================================================================
# Title:         test-mod_protein_profile.R
# Project:       ProSIFT (PROtein Statistical Integration and Filtering Tool)
# Author:        Reina Hastings (reinahastings13@gmail.com)
# Created:       2026-07-15
# Last Modified: 2026-07-15
# Purpose:       Group C reactive-logic checks for the Protein Profile View
#                module (frontend/R/mod_protein_profile.R) using testServer
#                (C7). Verifies the select-prompt / not-found empty states, that
#                a known protein renders all nine profile sections, that the
#                chem_data reactive matches the query, and that the returned
#                back / term_selected navigation events fire.
# Inputs:        The synthetic fixture (helper-fixture.R); mod server (sourced).
# Outputs:       testthat assertions.
# Usage:         Rscript frontend/run_tests.R    (or testthat::test_file(...))
# =============================================================================

library(shiny)

# Render a renderUI output value to a single HTML string for grepl checks.
.as_html <- function(x) paste(as.character(x), collapse = '\n')


# --- C7: profile module states, sections, and navigation events --------------

test_that('C7: profile handles empty states, renders sections, and returns events', {
  con <- fixture_con()
  on.exit(close_results_db(con))
  sel <- reactiveVal(NULL)

  testServer(mod_protein_profile_server,
             args = list(con = reactive(con), selected_protein = sel), {
    # NULL selection -> the select-a-protein prompt.
    expect_match(.as_html(output$content), 'Select a protein', fixed = TRUE)

    # A known protein renders every profile section.
    sel('P1')
    session$flushReact()
    html_p1 <- .as_html(output$content)
    for (s in c('Differential abundance', 'Data quality',
                'Functional annotation', 'Enriched terms',
                'Disease associations', 'Druggability',
                'Chemical interactions', 'Literature co-occurrence')) {
      expect_true(grepl(s, html_p1, fixed = TRUE), info = s)
    }

    # chem_data reactive matches db_protein_chemicals (P1 has 2 success rows).
    expect_equal(nrow(chem_data()), 2L)

    # Unknown protein -> the not-found placeholder.
    sel('NOPE')
    session$flushReact()
    expect_match(.as_html(output$content), 'not found', fixed = TRUE)

    # Navigation events surface through the returned list.
    sel('P1')
    session$flushReact()
    session$setInputs(back_to_db = 1L)
    expect_equal(session$returned$back(), 1L)

    session$setInputs(term_link = 'GOBP_X')
    expect_equal(session$returned$term_selected(), 'GOBP_X')
  })
})


# --- Fix 2a: disabled-database state vs. found-nothing; fail-open -------------

test_that('Fix 2a: a disabled database shows the not-enabled state, not a zero', {
  # 'disabled' variant: dgidb absent from databases.enabled, drug table emptied.
  con <- fixture_con('disabled')
  on.exit(close_results_db(con))
  testServer(mod_protein_profile_server,
             args = list(con = reactive(con), selected_protein = reactive('P1')), {
    h <- .as_html(output$content)
    expect_match(h, 'DGIdb query was not enabled for this run', fixed = TRUE)
    expect_false(grepl('No drug interactions found', h, fixed = TRUE))
  })
})

test_that('Fix 2a: fail-open (no run_parameters) shows the normal drug state', {
  con <- fixture_con('no_params')
  on.exit(close_results_db(con))
  testServer(mod_protein_profile_server,
             args = list(con = reactive(con), selected_protein = reactive('P1')), {
    h <- .as_html(output$content)
    expect_false(grepl('was not enabled for this run', h, fixed = TRUE))
  })
})

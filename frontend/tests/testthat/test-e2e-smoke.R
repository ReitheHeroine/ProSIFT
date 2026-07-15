# =============================================================================
# Title:         test-e2e-smoke.R
# Project:       ProSIFT (PROtein Statistical Integration and Filtering Tool)
# Author:        Reina Hastings (reinahastings13@gmail.com)
# Created:       2026-07-15
# Last Modified: 2026-07-15
# Purpose:       Group D end-to-end smoke tests for the Module 08 Shiny app,
#                driven under headless Chrome via shinytest2. Formalises two of
#                the scratchpad smoke flows against the in-code synthetic
#                fixture (not the cluster databases): D2 (protein click-through
#                -> profile card -> back) and D4 (the protein <-> term
#                navigation triangle). Both skip on CRAN and skip gracefully
#                when Chrome / shinytest2 is unavailable, so the headless
#                Group A/B/C gate stays green regardless.
# Inputs:        The frontend app (PROSIFT_FRONTEND_ROOT) and a fixture DB
#                staged under a temp PROSIFT_RESULTS_DIR (build_fixture_db()).
# Outputs:       testthat assertions (browser-driven).
# Usage:         NOT_CRAN=true Rscript frontend/run_tests.R
# =============================================================================

# --- Guards + shared launch helper -------------------------------------------

.e2e_skip_checks <- function() {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed('shinytest2')
  testthat::skip_if_not_installed('chromote')
  chrome <- tryCatch(chromote::find_chrome(), error = function(e) NA_character_)
  testthat::skip_if(is.null(chrome) || is.na(chrome) || !nzchar(chrome),
                    'Chrome not available for headless E2E')
}

# Stage a fresh fixture DB where discover_runs() will find it, launch the app,
# and register cleanup. Returns the AppDriver. env is scoped so the temp dir and
# the PROSIFT_RESULTS_DIR override are torn down when the test frame exits.
.e2e_launch <- function(env, name) {
  Sys.setenv(NOT_CRAN = 'true')
  chromote::set_chrome_args(c(chromote::default_chrome_args(), '--no-sandbox'))

  tmp <- tempfile('prosift_e2e_')
  dir.create(file.path(tmp, 'SYN_RUN'), recursive = TRUE)
  build_fixture_db(file.path(tmp, 'SYN_RUN', 'prosift_results.db'))

  old_dir <- Sys.getenv('PROSIFT_RESULTS_DIR', unset = NA_character_)
  Sys.setenv(PROSIFT_RESULTS_DIR = tmp)
  withr::defer({
    if (is.na(old_dir)) Sys.unsetenv('PROSIFT_RESULTS_DIR')
    else Sys.setenv(PROSIFT_RESULTS_DIR = old_dir)
    unlink(tmp, recursive = TRUE)
  }, envir = env)

  app <- shinytest2::AppDriver$new(
    app_dir = PROSIFT_FRONTEND_ROOT, name = name,
    load_timeout = 60000, timeout = 25000)
  withr::defer(app$stop(), envir = env)
  app$wait_for_idle(2000, 30000)
  app
}


# --- D2: protein click-through -> profile card -> back -----------------------

test_that('D2: click-through opens the full profile card and back returns', {
  .e2e_skip_checks()
  app <- .e2e_launch(environment(), 'prosift_d2')

  # Startup chain: a run and a contrast are populated from the fixture.
  expect_true(nzchar(app$get_value(input = 'run') %||% ''))
  app$set_inputs(contrast = 'KO_vs_WT')
  app$wait_for_idle(1000, 20000)

  # Click the first table row -> the app switches to the profile tab.
  app$set_inputs(`db-table_rows_selected` = 1, allow_no_input_binding_ = TRUE)
  app$wait_for_idle(1500, 25000)
  expect_identical(app$get_value(input = 'main_tabs'), 'profile')

  # All eight profile sections render.
  prof <- app$get_html('#profile-content') %||% ''
  for (sec in c('Differential abundance', 'Data quality',
                'Functional annotation', 'Enriched terms',
                'Disease associations', 'Druggability',
                'Chemical interactions', 'Literature co-occurrence')) {
    expect_true(grepl(sec, prof, fixed = TRUE), info = sec)
  }

  # Back link returns to the database tab.
  app$click('profile-back_to_db')
  app$wait_for_idle(1000, 15000)
  expect_identical(app$get_value(input = 'main_tabs'), 'db')
})


# --- D4: the protein <-> term navigation triangle ----------------------------

test_that('D4: member -> protein and profile-term-link -> term navigation works', {
  .e2e_skip_checks()
  app <- .e2e_launch(environment(), 'prosift_d4')

  app$set_inputs(contrast = 'KO_vs_WT')
  app$wait_for_idle(1000, 20000)

  # Open the biological-process tab -> the term list renders.
  app$set_inputs(main_tabs = 'process')
  app$wait_for_idle(1500, 25000)
  expect_true(grepl('terms', app$get_value(output = 'process-term_count') %||% ''))

  # Click a term row -> the term profile (stats + member table) renders.
  app$set_inputs(`process-term_table_rows_selected` = 1, allow_no_input_binding_ = TRUE)
  app$wait_for_idle(1500, 25000)
  view <- app$get_html('#process-view') %||% ''
  expect_true(grepl('Enrichment statistics', view, fixed = TRUE))
  expect_true(grepl('Member proteins', view, fixed = TRUE))

  # Term -> protein leg: click a member row -> the protein profile tab opens.
  app$set_inputs(`process-member_table_rows_selected` = 1, allow_no_input_binding_ = TRUE)
  app$wait_for_idle(1500, 25000)
  expect_identical(app$get_value(input = 'main_tabs'), 'profile')

  # Protein -> term leg: a profile enriched-term link opens that term's profile,
  # driven the way the onclick does (setInputValue with event priority).
  app$run_js("Shiny.setInputValue('profile-term_link', 'GOBP_X', {priority: 'event'});")
  app$wait_for_idle(1500, 25000)
  expect_identical(app$get_value(input = 'main_tabs'), 'process')
  view2 <- app$get_html('#process-view') %||% ''
  expect_true(grepl('Member proteins', view2, fixed = TRUE))

  # Back to the term list.
  app$click('process-back_to_list')
  app$wait_for_idle(1200, 20000)
  expect_true(grepl('terms', app$get_value(output = 'process-term_count') %||% ''))
})

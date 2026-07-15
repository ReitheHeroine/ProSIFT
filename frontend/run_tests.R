#!/usr/bin/env Rscript
# =============================================================================
# Title:         run_tests.R
# Project:       ProSIFT (PROtein Statistical Integration and Filtering Tool)
# Author:        Reina Hastings (reinahastings13@gmail.com)
# Created:       2026-07-15
# Last Modified: 2026-07-15
# Purpose:       Runner for the Module 08 frontend test suite. Sources the R/
#                function libraries (via the testthat helpers) and runs every
#                testthat file under frontend/tests/testthat. Groups A/B/C run
#                headless (testthat + DBI/RSQLite + shiny only); Group D
#                (shinytest2 E2E) runs only when NOT_CRAN=true and Chrome is
#                available, and skips gracefully otherwise.
# Inputs:        None. Optional env vars: NOT_CRAN=true to enable the E2E tests;
#                PROSIFT_RESULTS_DIR is set internally by the E2E fixture.
# Outputs:       A testthat summary to stdout; non-zero exit on any failure.
# Usage:         Rscript frontend/run_tests.R
#                copy/paste:  NOT_CRAN=true Rscript frontend/run_tests.R
# =============================================================================

# --- Locate this script's directory (the frontend root) ---------------------

# Resolve the frontend dir from the --file= argument so the runner works from
# any launch CWD; fall back to the current directory.
.args <- commandArgs(trailingOnly = FALSE)
.file_arg <- sub('^--file=', '', .args[grepl('^--file=', .args)])
frontend_dir <- if (length(.file_arg) == 1 && nzchar(.file_arg)) {
  dirname(normalizePath(.file_arg))
} else {
  normalizePath(getwd())
}

# --- Run the suite ----------------------------------------------------------

library(testthat)

test_path <- file.path(frontend_dir, 'tests', 'testthat')
message('Running ProSIFT Module 08 frontend tests from: ', test_path)

results <- test_dir(test_path, reporter = 'summary', stop_on_failure = FALSE)

# --- Exit code reflects failures --------------------------------------------

df <- as.data.frame(results)
n_fail <- sum(df$failed) + sum(df$error)
if (n_fail > 0) quit(status = 1, save = 'no') else quit(status = 0, save = 'no')

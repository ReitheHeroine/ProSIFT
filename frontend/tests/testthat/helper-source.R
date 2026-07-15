# =============================================================================
# Title:         helper-source.R
# Project:       ProSIFT (PROtein Statistical Integration and Filtering Tool)
# Author:        Reina Hastings (reinahastings13@gmail.com)
# Created:       2026-07-15
# Last Modified: 2026-07-15
# Purpose:       testthat helper that makes the Module 08 frontend functions
#                available to the test suite. The app in frontend/ is a plain
#                Shiny app (no DESCRIPTION), so the R/ function libraries are
#                sourced here into the global environment before any test runs.
# Inputs:        None (locates frontend/R relative to the test working dir).
# Outputs:       Sourced side effect: db.R / ui_helpers.R / config.R / mod_*.R
#                functions defined globally; PROSIFT_FRONTEND_ROOT set globally.
# Usage:         Auto-sourced by testthat::test_dir() (files matching ^helper).
# =============================================================================

# --- Locate the frontend root and source the R/ function libraries ----------

local({
  # Walk up from the test working directory until we find R/db.R. testthat runs
  # tests with the working directory set to tests/testthat, so the frontend root
  # is normally two levels up, but the search keeps this robust to the launch
  # context (run_tests.R sets a different CWD than an ad-hoc test_file() call).
  find_root <- function(start) {
    d <- normalizePath(start, winslash = '/', mustWork = FALSE)
    for (i in seq_len(8)) {
      if (file.exists(file.path(d, 'R', 'db.R'))) return(d)
      parent <- dirname(d)
      if (identical(parent, d)) break
      d <- parent
    }
    stop('helper-source: could not locate frontend/R from ', start)
  }

  root <- find_root(getwd())
  assign('PROSIFT_FRONTEND_ROOT', root, envir = globalenv())

  # Every file in R/ is a pure function/constant library (no top-level side
  # effects, no shinyApp() call), so sourcing them all is safe and gives the
  # tests the db.R queries, the ui_helpers formatters, the config constants,
  # and the module server functions.
  r_files <- list.files(file.path(root, 'R'), pattern = '\\.R$', full.names = TRUE)
  for (f in r_files) sys.source(f, envir = globalenv())
})

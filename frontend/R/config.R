# =============================================================================
# Title:         config.R
# Project:       ProSIFT (PROtein Statistical Integration and Filtering Tool)
# Author:        Reina Hastings (reinahastings13@gmail.com)
# Created:       2026-07-13
# Last Modified: 2026-07-13
# Purpose:       Module 08 frontend configuration. Holds shared constants (brand
#                accent, display defaults) and the results-directory / run
#                discovery logic that scans for Module 07 SQLite databases.
# Inputs:        Optional env var PROSIFT_RESULTS_DIR (a directory to scan for
#                prosift_results.db files). Falls back to the local dev database
#                and the project results/ directory.
# Outputs:       Functions and constants sourced by app.R (auto-loaded from R/).
# Usage:         Sourced automatically by Shiny when the app in frontend/ starts.
# =============================================================================

# --- Brand / display constants ---------------------------------------------

# Teal accent established from the 2026-04-14 mockup (logo, active tabs,
# leading-edge annotations). See Module 08 spec Section 4.4.3.
PROSIFT_TEAL <- '#1D9E75'

# Default rows per page in the protein table.
PROSIFT_PAGE_LENGTH <- 25L

# The database filename every Module 07 run publishes.
PROSIFT_DB_FILENAME <- 'prosift_results.db'


# --- Run discovery ----------------------------------------------------------

#' Candidate root directories to scan for result databases.
#'
#' Order of precedence: an explicit PROSIFT_RESULTS_DIR env var, then the local
#' dev database staged under the app, then the project-level results/ tree. Only
#' directories that exist are returned. Paths are resolved relative to the app
#' directory (`app_dir`) so the app works regardless of the launch CWD.
prosift_result_roots <- function(app_dir = getwd()) {
  env_dir <- Sys.getenv('PROSIFT_RESULTS_DIR', unset = '')
  candidates <- c(
    if (nzchar(env_dir)) env_dir,
    file.path(app_dir, '.devdata'),
    normalizePath(file.path(app_dir, '..', 'results'), mustWork = FALSE),
    # Cluster runs rsynced in (real Module 07 outputs).
    normalizePath(file.path(app_dir, '..', 'results_cluster'), mustWork = FALSE)
  )
  candidates[dir.exists(candidates)]
}

#' Discover available runs by scanning for prosift_results.db files.
#'
#' Returns a data.frame with columns `run_id` (label, read from each database's
#' run_metadata table when available, else the parent directory name) and `path`
#' (absolute path to the database file). Empty data.frame if none found.
discover_runs <- function(roots = prosift_result_roots()) {
  empty <- data.frame(run_id = character(0), path = character(0),
                      stringsAsFactors = FALSE)
  if (length(roots) == 0) return(empty)

  # Step 1: collect every prosift_results.db under the candidate roots.
  db_paths <- unlist(lapply(roots, function(r) {
    list.files(r, pattern = paste0('^', PROSIFT_DB_FILENAME, '$'),
               recursive = TRUE, full.names = TRUE)
  }), use.names = FALSE)
  db_paths <- unique(normalizePath(db_paths, mustWork = FALSE))
  if (length(db_paths) == 0) return(empty)

  # Step 2: label each database by its run_metadata.run_id, falling back to the
  # parent directory name if that table or column is unavailable.
  run_ids <- vapply(db_paths, function(p) {
    label <- tryCatch({
      con <- open_results_db(p)
      on.exit(close_results_db(con), add = TRUE)
      val <- DBI::dbGetQuery(con, 'SELECT run_id FROM run_metadata LIMIT 1')$run_id
      if (length(val) == 1 && nzchar(val)) val else NA_character_
    }, error = function(e) NA_character_)
    if (is.na(label)) basename(dirname(p)) else label
  }, character(1), USE.NAMES = FALSE)

  # Step 3: disambiguate duplicate labels so the selector never shows two
  # entries with the same name pointing at different databases (e.g. the same
  # run present in both results/ and .devdata/). Tag duplicates with their
  # parent directory; if those still collide, fall back to the full directory.
  if (anyDuplicated(run_ids)) {
    is_dup <- run_ids %in% run_ids[duplicated(run_ids)]
    labels <- ifelse(is_dup,
                     paste0(run_ids, ' [', basename(dirname(db_paths)), ']'),
                     run_ids)
    if (anyDuplicated(labels)) {
      still <- labels %in% labels[duplicated(labels)]
      labels[still] <- paste0(run_ids[still], ' [', dirname(db_paths)[still], ']')
    }
    run_ids <- labels
  }

  data.frame(run_id = run_ids, path = db_paths, stringsAsFactors = FALSE)
}

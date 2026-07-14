# =============================================================================
# Title:         ui_helpers.R
# Project:       ProSIFT (PROtein Statistical Integration and Filtering Tool)
# Author:        Reina Hastings (reinahastings13@gmail.com)
# Created:       2026-07-14
# Last Modified: 2026-07-14
# Purpose:       Small shared UI formatters and building blocks used by more
#                than one Module 08 view module (profile card, biological
#                process view): value formatters, coloured badges, stat cards,
#                and section wrappers. Styling lives in www/prosift.css.
# Inputs:        Scalars from query results (formatters expect length-1 values).
# Outputs:       Character strings or shiny tag objects; sourced by Shiny.
# =============================================================================

# --- Value formatters (expect scalars) --------------------------------------

fmt_fc <- function(x) if (is.na(x)) '-' else sprintf('%+.2f', x)
fmt_p  <- function(x) if (is.na(x)) '-' else formatC(x, format = 'g', digits = 2)

na_dash <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(as.character(x))) {
    '-'
  } else {
    as.character(x)
  }
}

# MSigDB-style term names -> readable ("GOBP_ADAPTIVE_THERMOGENESIS" ->
# "adaptive thermogenesis").
clean_term_name <- function(x) {
  if (is.na(x)) return('-')
  x <- sub('^(GOBP|GOCC|GOMF|REACTOME|KEGG|WP|HP)_', '', x)
  tolower(gsub('_', ' ', x))
}


# --- Coloured badges (styling in prosift.css) -------------------------------

dir_badge <- function(d) {
  d <- if (is.na(d)) 'ns' else d
  cls <- switch(d, up = 'badge-up', down = 'badge-down', 'badge-ns')
  shiny::span(class = paste('badge', cls), d)
}

det_badge <- function(cat) {
  shiny::span(class = 'badge det-badge', na_dash(cat))
}


# --- Layout building blocks -------------------------------------------------

stat_card <- function(label, value) {
  shiny::div(class = 'stat-card',
    shiny::div(class = 'stat-label', label),
    shiny::div(class = 'stat-value', value))
}

section <- function(title, ...) {
  shiny::div(class = 'section',
    shiny::div(class = 'section-title', title),
    shiny::tagList(...))
}

empty_note <- function(msg) shiny::div(class = 'empty-note', msg)

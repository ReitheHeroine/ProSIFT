# =============================================================================
# Title:         mod_protein_db.R
# Project:       ProSIFT (PROtein Statistical Integration and Filtering Tool)
# Author:        Reina Hastings (reinahastings13@gmail.com)
# Created:       2026-07-13
# Last Modified: 2026-07-13
# Purpose:       Module 08 Protein Database View (spec Section 4.1). A Shiny
#                module rendering the sortable/filterable overview table of all
#                proteins for the selected contrast, with click-through to a
#                protein profile. Presentation only; all SQL lives in db.R.
# Inputs:        protein_data - a reactive returning the db_protein_table()
#                data.frame for the current run + contrast.
# Outputs:       UI + server. The server returns a reactive giving the
#                protein_id of the clicked row (NULL when none selected).
# Usage:         mod_protein_db_ui('db'); sel <- mod_protein_db_server('db', dat)
# =============================================================================

# --- UI ---------------------------------------------------------------------

mod_protein_db_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::div(
      class = 'view-controls',
      shiny::checkboxInput(ns('sig_only'), 'Significant only', value = FALSE),
      shiny::span(class = 'row-count', shiny::textOutput(ns('row_count'), inline = TRUE))
    ),
    DT::DTOutput(ns('table'))
  )
}


# --- Server -----------------------------------------------------------------

mod_protein_db_server <- function(id, protein_data) {
  shiny::moduleServer(id, function(input, output, session) {

    # Step 1: shape the raw query result into display columns. Kept separate
    # from the significance filter so the badge/format logic runs once.
    display_data <- shiny::reactive({
      df <- protein_data()
      shiny::req(df)
      # Column order here is the display order (identity, statistics, QC, then
      # the annotation-richness indicators). protein_id and significant are
      # hidden helpers (click-through + the significance filter).
      data.frame(
        protein_id   = df$protein_id,
        # Fall back to the accession when a protein has no gene symbol (matches
        # the profile-header fallback), so the cell is never blank.
        Gene         = ifelse(is.na(df$gene_symbol) | !nzchar(df$gene_symbol),
                              df$protein_id, df$gene_symbol),
        `Protein name` = df$protein_name,
        `Adj p-value` = df$adj_pvalue,
        `log2 FC`    = df$log2_fc,
        Direction    = ifelse(is.na(df$direction), 'ns', df$direction),
        Detection    = df$detection_category,
        `% imputed`  = round(df$imputation_fraction * 100, 1),
        Diseases     = df$disease_count,
        Drugs        = ifelse(!is.na(df$drug_count) & df$drug_count > 0, 'Yes', '--'),
        PubMed       = df$top_pmi,
        significant  = ifelse(is.na(df$significant), 0L, df$significant),
        check.names  = FALSE,
        stringsAsFactors = FALSE
      )
    })

    # Step 2: apply the "significant only" toggle.
    filtered_data <- shiny::reactive({
      d <- display_data()
      if (isTRUE(input$sig_only)) d <- d[d$significant == 1L, , drop = FALSE]
      d
    })

    output$row_count <- shiny::renderText({
      sprintf('%d proteins', nrow(filtered_data()))
    })

    # Step 3: render the DT table. Column 0 (protein_id) and the trailing
    # `significant` helper column are hidden; they back click-through and the
    # filter but are not shown. Server-side processing (the renderDT default)
    # keeps ~8,000 rows responsive (spec Section 4.1.4).
    output$table <- DT::renderDT({
      d <- filtered_data()
      hide_targets <- which(names(d) %in% c('protein_id', 'significant')) - 1L
      # Render NA in the richness columns as '--' at display time only, so the
      # underlying value stays numeric and sorts correctly (spec Section 4.1.1).
      dash_na <- DT::JS(
        "function(data, type){return type === 'display' && data === null ? '--' : data;}")
      pmi_na <- DT::JS(paste0(
        "function(data, type){return type === 'display' ? ",
        "(data === null ? '--' : Number(data).toFixed(2)) : data;}"))
      diseases_col <- which(names(d) == 'Diseases') - 1L
      pubmed_col <- which(names(d) == 'PubMed') - 1L
      adjp_col <- which(names(d) == 'Adj p-value') - 1L
      # Centre every column except the two text identifiers (and the hidden
      # helpers), which read better left-aligned.
      center_cols <- which(!names(d) %in%
        c('protein_id', 'Gene', 'Protein name', 'significant')) - 1L

      dt <- DT::datatable(
        d,
        rownames = FALSE,
        filter = 'top',
        selection = 'single',
        class = 'stripe hover row-border',
        options = list(
          pageLength = PROSIFT_PAGE_LENGTH,
          lengthMenu = list(c(25, 50, 100, -1), c('25', '50', '100', 'All')),
          scrollX = TRUE,
          # Persist sort/search/page/column state so it survives tab switches
          # and the view re-rendering (spec Section 3.2).
          stateSave = TRUE,
          columnDefs = list(
            list(visible = FALSE, targets = hide_targets),
            list(className = 'dt-center', targets = center_cols),
            list(targets = diseases_col, render = dash_na),
            list(targets = pubmed_col, render = pmi_na)),
          # Default to most-significant-first (adjusted p-value ascending).
          order = list(list(adjp_col, 'asc'))
        )
      )

      # log2 FC: green when up, red when down (spec Section 4.1.1).
      dt <- DT::formatStyle(
        dt, 'log2 FC',
        color = DT::styleInterval(0, c('#A32D2D', '#0F6E56'))
      )
      # Direction badge colouring via text style (avoids raw HTML in the table).
      dt <- DT::formatStyle(
        dt, 'Direction', fontWeight = 'bold',
        color = DT::styleEqual(c('up', 'down', 'ns'),
                               c('#0F6E56', '#A32D2D', '#8A8F98'))
      )
      # Adjusted p-value: two significant figures, numeric sort preserved.
      dt <- DT::formatSignif(dt, 'Adj p-value', digits = 2)
      dt
    }, server = TRUE)

    # Step 4: expose the clicked protein_id to the parent for navigation. Depend
    # only on the selection event; read filtered_data() with isolate() so that
    # toggling a filter with a row selected does not fire a spurious click.
    shiny::reactive({
      sel <- input$table_rows_selected
      if (length(sel) == 0) return(NULL)
      shiny::isolate(filtered_data())$protein_id[sel]
    })
  })
}

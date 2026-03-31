/*
 * ProSIFT -- modules/local/enrichment/main.nf
 * Process: ENRICHMENT
 * Module 05: Enrichment Analysis.
 * Runs ORA (gseapy.enrich) and preranked GSEA (gseapy.prerank) against local
 * MSigDB GMT files. Produces a unified enrichment results table, a protein-term
 * mapping table, lollipop plots (PNG + HTML) per contrast per library, GSEA
 * running score plots (PNG) for top significant terms, and a summary text file.
 */

nextflow.enable.dsl = 2

process ENRICHMENT {

    tag "$meta.run_id"

    publishDir "${params.outdir}/${meta.run_id}/enrichment", mode: 'copy'

    input:
    tuple val(meta),
          path(diff_abundance_results),
          path(params_yml)

    output:
    tuple val(meta), path("*.enrichment_results.parquet"),      emit: enrichment_results
    tuple val(meta), path("*.enrichment_results.csv"),          emit: enrichment_results_csv
    tuple val(meta), path("*.protein_term_mapping.parquet"),    emit: protein_term_mapping
    tuple val(meta), path("*.enrichment_summary.txt"),          emit: summary_txt
    // Plots: one pair per contrast per library (glob matches all combinations)
    tuple val(meta), path("*.ora_lollipop.png"),                emit: ora_lollipop_png
    tuple val(meta), path("*.ora_lollipop.html"),               emit: ora_lollipop_html
    tuple val(meta), path("*.gsea_lollipop.png"),               emit: gsea_lollipop_png
    tuple val(meta), path("*.gsea_lollipop.html"),              emit: gsea_lollipop_html
    tuple val(meta), path("*.gsea_running_score.*.png"),        emit: gsea_running_score_png

    script:
    """
    enrichment.py \\
        --results ${diff_abundance_results} \\
        --params  ${params_yml} \\
        --run-id  ${meta.run_id} \\
        --outdir  .
    """

    stub:
    """
    # Extract the first contrast name from params.yml (same approach as Module 04)
    contrast=\$(grep -A 20 'contrasts:' '${params_yml}' \\
                | grep -m1 '^[[:space:]]*-' \\
                | sed 's/^[[:space:]]*-[[:space:]]*//' \\
                | tr -d '"')

    touch ${meta.run_id}.enrichment_results.parquet
    touch ${meta.run_id}.enrichment_results.csv
    touch ${meta.run_id}.protein_term_mapping.parquet
    touch ${meta.run_id}.enrichment_summary.txt
    touch ${meta.run_id}.\${contrast}.GO_BP.ora_lollipop.png
    touch ${meta.run_id}.\${contrast}.GO_BP.ora_lollipop.html
    touch ${meta.run_id}.\${contrast}.GO_BP.gsea_lollipop.png
    touch ${meta.run_id}.\${contrast}.GO_BP.gsea_lollipop.html
    touch ${meta.run_id}.\${contrast}.GO_BP.gsea_running_score.STUB_TERM.png
    """

}

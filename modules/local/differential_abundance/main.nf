/*
 * ProSIFT -- modules/local/differential_abundance/main.nf
 * Process: DIFFERENTIAL_ABUNDANCE
 * Module 04: Differential Abundance.
 * Parses and validates contrasts from params.yml, fits a linear model per
 * protein using limma + empirical Bayes via rpy2, applies DEqMS peptide-
 * count-aware correction, and produces a per-protein results table, text
 * summary, and volcano / MA diagnostic plots (PNG + HTML) per contrast.
 */

nextflow.enable.dsl = 2

process DIFFERENTIAL_ABUNDANCE {

    tag "$meta.run_id"

    publishDir "${params.outdir}/${meta.run_id}/differential_abundance", mode: 'copy'

    input:
    tuple val(meta),
          path(imputed_matrix),
          path(validated_metadata),
          path(id_mapping),
          path(params_yml)

    output:
    // Primary output -- PRIMARY results table consumed by Modules 05/06/07.
    // NB: the primary globs must NOT match the sensitivity files, which use a
    // `.sensitivity.*` infix (module spec Section 2.2 filename constraint).
    tuple val(meta), path("*.diff_abundance_results.parquet"), emit: results_table
    tuple val(meta), path("*.diff_abundance_results.csv"),     emit: results_table_csv
    tuple val(meta), path("*.diff_abundance_summary.txt"),     emit: summary_txt
    tuple val(meta), path("*.analysis_provenance.txt"),        emit: provenance
    // Sensitivity table (dual reporting): written only when quarantine is active.
    tuple val(meta), path("*.diff_abundance_results.sensitivity.parquet"), emit: results_table_sensitivity, optional: true
    tuple val(meta), path("*.diff_abundance_results.sensitivity.csv"),     emit: results_table_sensitivity_csv, optional: true
    // Plots: one pair per contrast (glob matches all contrasts). The `.sensitivity`
    // plot globs are separate and optional so the primary globs stay single-file.
    tuple val(meta), path("*.volcano_plot.png"),               emit: volcano_png
    tuple val(meta), path("*.volcano_plot.html"),              emit: volcano_html
    tuple val(meta), path("*.ma_plot.png"),                    emit: ma_png
    tuple val(meta), path("*.ma_plot.html"),                   emit: ma_html
    tuple val(meta), path("*.volcano_plot.sensitivity.png"),   emit: volcano_png_sensitivity,  optional: true
    tuple val(meta), path("*.volcano_plot.sensitivity.html"),  emit: volcano_html_sensitivity, optional: true
    tuple val(meta), path("*.ma_plot.sensitivity.png"),        emit: ma_png_sensitivity,       optional: true
    tuple val(meta), path("*.ma_plot.sensitivity.html"),       emit: ma_html_sensitivity,      optional: true

    script:
    """
    differential_abundance.py \\
        --matrix     ${imputed_matrix} \\
        --metadata   ${validated_metadata} \\
        --id-mapping ${id_mapping} \\
        --params     ${params_yml} \\
        --run-id     ${meta.run_id} \\
        --outdir     .
    """

    stub:
    """
    # Extract the first contrast name from params.yml so stub filenames match
    # the real output filenames. Strips leading whitespace, '- ', and quotes.
    contrast=\$(grep -A 20 'contrasts:' '${params_yml}' \\
                | grep -m1 '^[[:space:]]*-' \\
                | sed 's/^[[:space:]]*-[[:space:]]*//' \\
                | tr -d '"')

    # Always-on outputs (non-optional emits). Sensitivity + sensitivity-plot
    # outputs are optional and intentionally omitted from the stub.
    touch ${meta.run_id}.diff_abundance_results.parquet
    touch ${meta.run_id}.diff_abundance_results.csv
    touch ${meta.run_id}.diff_abundance_summary.txt
    touch ${meta.run_id}.analysis_provenance.txt
    touch ${meta.run_id}.\${contrast}.volcano_plot.png
    touch ${meta.run_id}.\${contrast}.volcano_plot.html
    touch ${meta.run_id}.\${contrast}.ma_plot.png
    touch ${meta.run_id}.\${contrast}.ma_plot.html
    """

}

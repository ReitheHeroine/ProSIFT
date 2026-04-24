/*
 * ProSIFT -- modules/local/missingness_report/main.nf
 * Process: MISSINGNESS_REPORT
 * Module 01, Process 4.10: Missingness visualization report.
 * Reads FILTER_PROTEINS and VALIDATE_INPUTS outputs to produce four
 * interactive plots characterizing pre-filter missingness patterns.
 * Advisory only -- errorStrategy 'ignore' so pipeline continues on failure.
 */

nextflow.enable.dsl = 2

process MISSINGNESS_REPORT {

    tag "$meta.run_id"

    // Advisory process: failure does not stop the pipeline.
    errorStrategy 'ignore'

    publishDir "${params.outdir}/${meta.run_id}/validation",
        mode: 'copy'

    input:
    tuple val(meta), path(filter_table), path(validated_matrix), path(validated_metadata), path(params_yml)

    output:
    tuple val(meta), path("*.missingness_report.html"),   emit: report,   optional: true
    tuple val(meta), path("*.missingness_plots/"),        emit: plot_dir, optional: true

    script:
    """
    missingness_report.py \\
        --filter-table ${filter_table} \\
        --matrix       ${validated_matrix} \\
        --metadata     ${validated_metadata} \\
        --params       ${params_yml} \\
        --run-id       ${meta.run_id} \\
        --outdir       .
    """

    stub:
    """
    touch ${meta.run_id}.missingness_report.html
    mkdir -p ${meta.run_id}.missingness_plots
    touch ${meta.run_id}.missingness_plots/${meta.run_id}.filter_categories.png
    touch ${meta.run_id}.missingness_plots/${meta.run_id}.per_sample_detection.png
    touch ${meta.run_id}.missingness_plots/${meta.run_id}.missingness_heatmap.png
    touch ${meta.run_id}.missingness_plots/${meta.run_id}.missingness_histogram.png
    """

}

/*
 * ProSIFT -- modules/local/impute/main.nf
 * Process: IMPUTE
 * Module 03: Imputation.
 * Classifies missing values as MNAR or MAR using Module 01's detection
 * filter table (DEP-style: SINGLE-GROUP proteins = MNAR, PASSED proteins
 * with sporadic NaN = MAR). Applies MinProb for MNAR and protein-wise KNN
 * for MAR. Produces a fully imputed matrix (no NaN), per-protein imputation
 * summary, and two diagnostic plots.
 */

nextflow.enable.dsl = 2

process IMPUTE {

    tag "$meta.run_id"

    // Outputs published to the same normalization/ subdirectory as NORMALIZE.
    // imputed_matrix.parquet is the primary artifact consumed by Module 04.
    publishDir "${params.outdir}/${meta.run_id}/normalization", mode: 'copy'

    input:
    tuple val(meta),
          path(normalized_matrix),
          path(validated_metadata),
          path(filter_table),
          path(params_yml)

    output:
    // Primary output -- fed into Module 04 Differential Abundance
    tuple val(meta), path("*.imputed_matrix.parquet"),          emit: imputed_matrix
    // Per-cell imputation mask (observed/mnar/mar) -- used by Module 03b QC Report Assembly
    tuple val(meta), path("*.imputation_mask.parquet"),         emit: imputation_mask
    // Diagnostic data files
    tuple val(meta), path("*.imputation_summary.parquet"),      emit: imputation_summary
    tuple val(meta), path("*.imputation_summary.csv"),          emit: imputation_summary_csv
    tuple val(meta), path("*.imputation_summary.txt"),          emit: imputation_summary_txt
    // Plots (PNG + interactive HTML, 2 plot types)
    tuple val(meta), path("*.imputation_distributions.png"),    emit: plot_distributions_png
    tuple val(meta), path("*.imputation_distributions.html"),   emit: plot_distributions_html
    tuple val(meta), path("*.imputation_fractions.png"),        emit: plot_fractions_png
    tuple val(meta), path("*.imputation_fractions.html"),       emit: plot_fractions_html

    script:
    """
    impute.py \\
        --matrix       ${normalized_matrix} \\
        --metadata     ${validated_metadata} \\
        --filter-table ${filter_table} \\
        --params       ${params_yml} \\
        --run-id       ${meta.run_id} \\
        --outdir       .
    """

    stub:
    """
    touch ${meta.run_id}.imputed_matrix.parquet
    touch ${meta.run_id}.imputation_mask.parquet
    touch ${meta.run_id}.imputation_summary.parquet
    touch ${meta.run_id}.imputation_summary.csv
    touch ${meta.run_id}.imputation_summary.txt
    touch ${meta.run_id}.imputation_distributions.png
    touch ${meta.run_id}.imputation_distributions.html
    touch ${meta.run_id}.imputation_fractions.png
    touch ${meta.run_id}.imputation_fractions.html
    """

}

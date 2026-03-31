/*
 * ProSIFT -- modules/local/normalize/main.nf
 * Process: NORMALIZE
 * Module 03: Normalization.
 * Applies the normalization method specified in params.yml (median, quantile,
 * VSN, or none) to the filtered protein abundance matrix. Produces a
 * normalized matrix in Parquet format, per-protein CV summaries, a plain-text
 * normalization summary, and three post-normalization diagnostic plots
 * (intensity distributions, PCA, sample correlation heatmap).
 */

nextflow.enable.dsl = 2

process NORMALIZE {

    tag "$meta.run_id"

    // Outputs published to a dedicated normalization/ subdirectory.
    // normalized_matrix.parquet feeds directly into IMPUTE (Module 03b).
    // cv_summary and plot files are user-facing diagnostics.
    publishDir "${params.outdir}/${meta.run_id}/normalization", mode: 'copy'

    input:
    tuple val(meta),
          path(filtered_matrix),
          path(validated_metadata),
          path(params_yml)

    output:
    // Primary data output -- fed into IMPUTE
    tuple val(meta), path("*.normalized_matrix.parquet"),        emit: normalized_matrix
    // Diagnostic data files
    tuple val(meta), path("*.normalization_summary.txt"),        emit: normalization_summary
    tuple val(meta), path("*.cv_summary.parquet"),               emit: cv_summary
    tuple val(meta), path("*.cv_summary.csv"),                   emit: cv_summary_csv
    tuple val(meta), path("*.postnorm_pca_results.parquet"),     emit: pca_results
    tuple val(meta), path("*.postnorm_correlation_matrix.parquet"), emit: correlation_matrix
    // Plots (PNG + interactive HTML, 3 plot types)
    tuple val(meta), path("*.postnorm_boxplots.png"),            emit: plot_boxplots_png
    tuple val(meta), path("*.postnorm_boxplots.html"),           emit: plot_boxplots_html
    tuple val(meta), path("*.postnorm_pca_scatter.png"),         emit: plot_pca_png
    tuple val(meta), path("*.postnorm_pca_scatter.html"),        emit: plot_pca_html
    tuple val(meta), path("*.postnorm_correlation_heatmap.png"), emit: plot_corr_png
    tuple val(meta), path("*.postnorm_correlation_heatmap.html"), emit: plot_corr_html

    script:
    """
    normalize.py \\
        --matrix   ${filtered_matrix} \\
        --metadata ${validated_metadata} \\
        --params   ${params_yml} \\
        --run-id   ${meta.run_id} \\
        --outdir   .
    """

    stub:
    """
    touch ${meta.run_id}.normalized_matrix.parquet
    touch ${meta.run_id}.normalization_summary.txt
    touch ${meta.run_id}.cv_summary.parquet
    touch ${meta.run_id}.cv_summary.csv
    touch ${meta.run_id}.postnorm_pca_results.parquet
    touch ${meta.run_id}.postnorm_correlation_matrix.parquet
    touch ${meta.run_id}.postnorm_boxplots.png
    touch ${meta.run_id}.postnorm_boxplots.html
    touch ${meta.run_id}.postnorm_pca_scatter.png
    touch ${meta.run_id}.postnorm_pca_scatter.html
    touch ${meta.run_id}.postnorm_correlation_heatmap.png
    touch ${meta.run_id}.postnorm_correlation_heatmap.html
    """

}

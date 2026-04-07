/*
 * ProSIFT -- modules/local/prenorm_qc/main.nf
 * Process: PRENORM_QC
 * Module 02: Pre-normalization QC/EDA.
 * Computes per-sample summaries, box plots, density plots, Q-Q plot, PCA,
 * and correlation clustermap on the post-filter pre-normalization data.
 * Produces per-sample outlier flags and a consolidated HTML QC report.
 * Does not modify any data.
 */

nextflow.enable.dsl = 2

process PRENORM_QC {

    tag "$meta.run_id"

    // All outputs published to a dedicated qc/ subdirectory.
    // Data files (parquet/csv) are inputs to QC Report Assembly (Phase 2).
    // PNG and HTML plots are user-facing diagnostics.
    publishDir "${params.outdir}/${meta.run_id}/qc", mode: 'copy'

    input:
    tuple val(meta),
          path(filtered_matrix),
          path(validated_metadata),
          path(id_mapping),
          path(params_yml)

    output:
    // Data files (used by downstream QC Report Assembly, Phase 2)
    tuple val(meta), path("*.sample_summary.parquet"),      emit: sample_summary
    tuple val(meta), path("*.sample_summary.csv"),          emit: sample_summary_csv
    tuple val(meta), path("*.sample_flags.parquet"),        emit: sample_flags
    tuple val(meta), path("*.sample_flags.csv"),            emit: sample_flags_csv
    tuple val(meta), path("*.correlation_matrix.parquet"),  emit: correlation_matrix
    tuple val(meta), path("*.pca_results.parquet"),         emit: pca_results
    // Report
    tuple val(meta), path("*.prenorm_qc_report.html"),      emit: qc_report
    // Plots (PNG + interactive HTML, 6 plot types)
    tuple val(meta), path("*.sample_intensity_summary.png"),  emit: plot_sample_summary_png
    tuple val(meta), path("*.sample_intensity_summary.html"), emit: plot_sample_summary_html
    tuple val(meta), path("*.intensity_boxplots.png"),        emit: plot_boxplots_png
    tuple val(meta), path("*.intensity_boxplots.html"),       emit: plot_boxplots_html
    tuple val(meta), path("*.density_plot.png"),              emit: plot_density_png
    tuple val(meta), path("*.density_plot.html"),             emit: plot_density_html
    tuple val(meta), path("*.qq_plot.png"),                   emit: plot_qq_png
    tuple val(meta), path("*.qq_plot.html"),                  emit: plot_qq_html
    tuple val(meta), path("*.pca_scatter.png"),               emit: plot_pca_png
    tuple val(meta), path("*.pca_scatter.html"),              emit: plot_pca_html
    tuple val(meta), path("*.correlation_heatmap.png"),       emit: plot_clustermap_png
    tuple val(meta), path("*.correlation_heatmap.html"),      emit: plot_clustermap_html

    script:
    """
    prenorm_qc.py \\
        --matrix   ${filtered_matrix} \\
        --metadata ${validated_metadata} \\
        --mapping  ${id_mapping} \\
        --params   ${params_yml} \\
        --run-id   ${meta.run_id} \\
        --outdir   .
    """

    stub:
    """
    touch ${meta.run_id}.sample_summary.parquet
    touch ${meta.run_id}.sample_summary.csv
    touch ${meta.run_id}.sample_flags.parquet
    touch ${meta.run_id}.sample_flags.csv
    touch ${meta.run_id}.correlation_matrix.parquet
    touch ${meta.run_id}.pca_results.parquet
    touch ${meta.run_id}.prenorm_qc_report.html
    touch ${meta.run_id}.sample_intensity_summary.png
    touch ${meta.run_id}.sample_intensity_summary.html
    touch ${meta.run_id}.intensity_boxplots.png
    touch ${meta.run_id}.intensity_boxplots.html
    touch ${meta.run_id}.density_plot.png
    touch ${meta.run_id}.density_plot.html
    touch ${meta.run_id}.qq_plot.png
    touch ${meta.run_id}.qq_plot.html
    touch ${meta.run_id}.pca_scatter.png
    touch ${meta.run_id}.pca_scatter.html
    touch ${meta.run_id}.correlation_heatmap.png
    touch ${meta.run_id}.correlation_heatmap.html
    """

}

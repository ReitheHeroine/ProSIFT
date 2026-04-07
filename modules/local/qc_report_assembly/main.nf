/*
 * ProSIFT -- modules/local/qc_report_assembly/main.nf
 * Process: QC_REPORT_ASSEMBLY
 * Module 03b: QC Report Assembly.
 * Compiles diagnostic outputs from Modules 01-03 into a single consolidated
 * HTML QC report. Re-renders all plots from upstream data files. Generates
 * one new visualization (CV density plot). Advisory checkpoint between data
 * preparation and statistical analysis.
 */

nextflow.enable.dsl = 2

process QC_REPORT_ASSEMBLY {

    tag "$meta.run_id"

    // Advisory process: failure does not block Module 04.
    errorStrategy 'ignore'

    publishDir "${params.outdir}/${meta.run_id}/qc", mode: 'copy'

    input:
    tuple val(meta),
          path(validated_metadata),
          path(detection_filter_table),
          path(filtered_matrix),
          path(sample_summary),
          path(sample_flags),
          path(normalized_matrix),
          path(cv_summary),
          path(normalization_summary),
          path(imputed_matrix),
          path(imputation_mask),
          path(imputation_summary),
          path(imputation_summary_txt),
          path(params_yml),
          path(missingness_report),
          path(prenorm_qc_report)

    output:
    tuple val(meta), path("*.qc_report.html"), emit: qc_report, optional: true

    script:
    """
    qc_report_assembly.py \\
        --metadata          ${validated_metadata} \\
        --filter-table      ${detection_filter_table} \\
        --filtered-matrix   ${filtered_matrix} \\
        --sample-summary    ${sample_summary} \\
        --sample-flags      ${sample_flags} \\
        --norm-matrix       ${normalized_matrix} \\
        --cv-summary        ${cv_summary} \\
        --norm-summary      ${normalization_summary} \\
        --imputed-matrix    ${imputed_matrix} \\
        --imp-mask          ${imputation_mask} \\
        --imp-summary       ${imputation_summary} \\
        --imp-summary-txt   ${imputation_summary_txt} \\
        --params            ${params_yml} \\
        --run-id            ${meta.run_id} \\
        --outdir            . \\
        --missingness-report ${missingness_report} \\
        --prenorm-report    ${prenorm_qc_report}
    """

    stub:
    """
    touch ${meta.run_id}.qc_report.html
    """

}

/*
 * ProSIFT -- modules/local/filter_proteins/main.nf
 * Process: FILTER_PROTEINS
 * Module 01, Process 4.4: Per-run detection filter and missingness overview.
 * Produces filtered matrix (downstream input), detection filter table
 * (diagnostic artifact), and partial validation report (Part 2 of 2).
 */

nextflow.enable.dsl = 2

process FILTER_PROTEINS {

    tag "$meta.run_id"

    // All three outputs are published: filtered matrix is the main input
    // for Module 02; table and report are diagnostic artifacts.
    publishDir "${params.outdir}/${meta.run_id}/validation",
        mode: 'copy'

    input:
    tuple val(meta), path(validated_matrix), path(validated_metadata), path(params_yml)

    output:
    tuple val(meta), path('*.filtered_matrix.parquet'),    emit: matrix
    tuple val(meta), path('*.detection_filter_table.csv'), emit: filter_table
    tuple val(meta), path('*.validation_report_part2.txt'), emit: report

    script:
    """
    filter_proteins.py \\
        --matrix   ${validated_matrix} \\
        --metadata ${validated_metadata} \\
        --params   ${params_yml} \\
        --run-id   ${meta.run_id} \\
        --outdir   .
    """

    stub:
    """
    touch ${meta.run_id}.filtered_matrix.parquet
    touch ${meta.run_id}.detection_filter_table.csv
    touch ${meta.run_id}.validation_report_part2.txt
    """

}

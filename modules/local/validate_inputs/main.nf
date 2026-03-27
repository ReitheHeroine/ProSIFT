/*
 * ProSIFT -- modules/local/validate_inputs/main.nf
 * Process: VALIDATE_INPUTS
 * Module 01, Processes 4.1-4.3: Validate matrix, validate metadata,
 * cross-validate samples. Produces validated Parquet intermediates and
 * a partial validation report (Part 1 of 2).
 */

nextflow.enable.dsl = 2

process VALIDATE_INPUTS {

    tag "$meta.run_id"

    // Only the report is user-facing; matrix/metadata are intermediates
    // passed directly to FILTER_PROTEINS via channel.
    publishDir "${params.outdir}/${meta.run_id}/validation",
        mode: 'copy',
        pattern: '*.validation_report_part1.txt'

    input:
    tuple val(meta), path(abundance_matrix), path(metadata), path(params_yml)

    output:
    tuple val(meta), path('*.validated_matrix.parquet'),    emit: matrix
    tuple val(meta), path('*.validated_metadata.parquet'),  emit: metadata
    tuple val(meta), path('*.validation_report_part1.txt'), emit: report

    script:
    """
    validate_inputs.py \\
        --abundance ${abundance_matrix} \\
        --metadata  ${metadata} \\
        --params    ${params_yml} \\
        --run-id    ${meta.run_id} \\
        --outdir    .
    """

    stub:
    """
    touch ${meta.run_id}.validated_matrix.parquet
    touch ${meta.run_id}.validated_metadata.parquet
    touch ${meta.run_id}.validation_report_part1.txt
    """

}

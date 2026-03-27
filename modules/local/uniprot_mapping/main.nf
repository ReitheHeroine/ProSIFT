/*
 * ProSIFT -- modules/local/uniprot_mapping/main.nf
 * Process: UNIPROT_MAPPING
 * Module 01, Processes 4.5-4.9: UniProt bulk ID mapping, edge case handling,
 * ortholog column stubs (Phase 3), and cache management.
 */

nextflow.enable.dsl = 2

process UNIPROT_MAPPING {

    tag "$meta.run_id"

    // Published outputs: the mapping Parquet (used by all downstream modules)
    // and the report are both user-facing.
    publishDir "${params.outdir}/${meta.run_id}/validation", mode: 'copy'

    input:
    tuple val(meta), path(filtered_matrix), path(params_yml)

    output:
    tuple val(meta), path("*.id_mapping.parquet"),             emit: mapping_table
    tuple val(meta), path("*.validation_report_part3.txt"),    emit: report

    script:
    // The cache directory lives in params.outdir (outside the Nextflow work dir)
    // so it persists across pipeline re-runs. This is intentional: the cache is
    // designed to survive between runs and avoid repeated UniProt API calls.
    def cachedir = "${params.outdir}/cache/id_mapping"
    """
    uniprot_mapping.py \\
        --matrix   ${filtered_matrix} \\
        --params   ${params_yml} \\
        --run_id   ${meta.run_id} \\
        --cachedir ${cachedir} \\
        --outdir   .
    """

    stub:
    """
    touch ${meta.run_id}.id_mapping.parquet
    touch ${meta.run_id}.validation_report_part3.txt
    """

}

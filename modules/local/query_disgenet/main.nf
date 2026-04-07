/*
 * ProSIFT -- modules/local/query_disgenet/main.nf
 * Process: QUERY_DISGENET
 * Module 06: Database Queries.
 * Queries DisGeNET REST API for gene-disease associations via human orthologs.
 * Proteins without orthologs are flagged and skipped. Rate limited (50/min
 * on Free Academic plan). Requires API key via environment variable.
 */

nextflow.enable.dsl = 2

process QUERY_DISGENET {

    tag "$meta.run_id"

    publishDir "${params.outdir}/${meta.run_id}/databases", mode: 'copy'

    input:
    tuple val(meta), path(mapping_table), path(params_yml)

    output:
    tuple val(meta), path("*.disgenet_associations.parquet"), emit: disgenet_associations

    script:
    def cachedir = "${params.outdir}/cache/databases/disgenet"
    """
    query_disgenet.py \\
        --mapping  ${mapping_table} \\
        --params   ${params_yml} \\
        --run-id   ${meta.run_id} \\
        --cachedir ${cachedir} \\
        --outdir   .
    """

    stub:
    """
    touch ${meta.run_id}.disgenet_associations.parquet
    """

}

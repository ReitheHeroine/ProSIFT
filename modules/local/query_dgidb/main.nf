/*
 * ProSIFT -- modules/local/query_dgidb/main.nf
 * Process: QUERY_DGIDB
 * Module 06: Database Queries.
 * Queries DGIdb GraphQL API for drug-gene interactions via human orthologs.
 * Proteins without orthologs are flagged and skipped. No API key required.
 * No documented rate limits.
 */

nextflow.enable.dsl = 2

process QUERY_DGIDB {

    tag "$meta.run_id"

    publishDir "${params.outdir}/${meta.run_id}/databases", mode: 'copy'

    input:
    tuple val(meta), path(mapping_table), path(params_yml)

    output:
    tuple val(meta), path("*.dgidb_interactions.parquet"), emit: dgidb_interactions

    script:
    def cachedir = "${params.outdir}/cache/databases/dgidb"
    """
    query_dgidb.py \\
        --mapping  ${mapping_table} \\
        --params   ${params_yml} \\
        --run-id   ${meta.run_id} \\
        --cachedir ${cachedir} \\
        --outdir   .
    """

    stub:
    """
    touch ${meta.run_id}.dgidb_interactions.parquet
    """

}

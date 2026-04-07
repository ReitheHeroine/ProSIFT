/*
 * ProSIFT -- modules/local/query_pubmed/main.nf
 * Process: QUERY_PUBMED
 * Module 06: Database Queries.
 * Queries NCBI ESearch API for literature co-occurrence of each protein with
 * user-defined search terms. Computes PMI-based normalized scores. One row per
 * protein per search term. Rate limited (10/sec with API key).
 */

nextflow.enable.dsl = 2

process QUERY_PUBMED {

    tag "$meta.run_id"

    publishDir "${params.outdir}/${meta.run_id}/databases", mode: 'copy'

    input:
    tuple val(meta), path(mapping_table), path(params_yml)

    output:
    tuple val(meta), path("*.pubmed_cooccurrence.parquet"), emit: pubmed_cooccurrence

    script:
    def cachedir = "${params.outdir}/cache/databases/pubmed"
    """
    query_pubmed.py \\
        --mapping  ${mapping_table} \\
        --params   ${params_yml} \\
        --run-id   ${meta.run_id} \\
        --cachedir ${cachedir} \\
        --outdir   .
    """

    stub:
    """
    touch ${meta.run_id}.pubmed_cooccurrence.parquet
    """

}

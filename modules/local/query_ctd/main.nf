/*
 * ProSIFT -- modules/local/query_ctd/main.nf
 * Process: QUERY_CTD
 * Module 06: Database Queries.
 * Filters pre-downloaded CTD bulk file (CTD_chem_gene_ixns.tsv.gz) for
 * chemical-gene interactions matching run proteins. Queries by both mouse
 * and human ortholog Entrez IDs. Fully offline after initial download.
 */

nextflow.enable.dsl = 2

process QUERY_CTD {

    tag "$meta.run_id"

    publishDir "${params.outdir}/${meta.run_id}/databases", mode: 'copy'

    input:
    tuple val(meta), path(mapping_table), path(params_yml)

    output:
    tuple val(meta), path("*.ctd_interactions.parquet"), emit: ctd_interactions

    script:
    def cachedir = "${params.outdir}/cache/databases/ctd"
    """
    query_ctd.py \\
        --mapping  ${mapping_table} \\
        --params   ${params_yml} \\
        --run-id   ${meta.run_id} \\
        --cachedir ${cachedir} \\
        --outdir   .
    """

    stub:
    """
    touch ${meta.run_id}.ctd_interactions.parquet
    """

}

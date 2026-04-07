/*
 * ProSIFT -- modules/local/query_uniprot/main.nf
 * Process: QUERY_UNIPROT
 * Module 06: Database Queries.
 * Queries UniProt REST API for curated protein annotations (name, function,
 * subcellular location, tissue expression, keywords). One row per protein.
 * Uses application-level caching for cross-run reuse.
 */

nextflow.enable.dsl = 2

process QUERY_UNIPROT {

    tag "$meta.run_id"

    publishDir "${params.outdir}/${meta.run_id}/databases", mode: 'copy'

    input:
    tuple val(meta), path(mapping_table), path(params_yml)

    output:
    tuple val(meta), path("*.uniprot_annotations.parquet"), emit: uniprot_annotations

    script:
    def cachedir = "${params.outdir}/cache/databases/uniprot"
    """
    query_uniprot.py \\
        --mapping  ${mapping_table} \\
        --params   ${params_yml} \\
        --run-id   ${meta.run_id} \\
        --cachedir ${cachedir} \\
        --outdir   .
    """

    stub:
    """
    touch ${meta.run_id}.uniprot_annotations.parquet
    """

}

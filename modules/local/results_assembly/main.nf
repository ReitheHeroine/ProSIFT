/*
 * ProSIFT -- modules/local/results_assembly/main.nf
 * Process: RESULTS_ASSEMBLY
 * Module 07: Results Assembly (SQLite).
 * The convergence point of the pipeline. Collects the Parquet/CSV outputs of
 * the analytical spine (Modules 01-05) and the database query layer (Module 06)
 * and assembles them into a single SQLite database (12 tables) plus CSV exports
 * and a human-readable assembly summary. Pure ETL: no statistics, no API calls.
 * The SQLite database is the sole data source for the Module 08 frontend.
 */

nextflow.enable.dsl = 2

process RESULTS_ASSEMBLY {

    tag "$meta.run_id"

    publishDir "${params.outdir}/${meta.run_id}", mode: 'copy'

    input:
    tuple val(meta),
          path(mapping_table),
          path(detection_table),
          path(sample_flags),
          path(imputed_matrix),
          path(imputation_mask),
          path(diff_abundance),
          path(enrichment_results),
          path(protein_term_mapping),
          path(uniprot_annotations),
          path(pubmed_cooccurrence),
          path(disgenet_associations),
          path(dgidb_interactions),
          path(ctd_interactions),
          path(params_yml)

    output:
    tuple val(meta), path('prosift_results.db'),               emit: database
    tuple val(meta), path('*.diff_abundance_results.csv'),     emit: diff_abundance_csv
    tuple val(meta), path('*.enrichment_results.csv'),         emit: enrichment_csv
    tuple val(meta), path('*.significant_proteins.csv'),       emit: significant_csv
    tuple val(meta), path('*.assembly_summary.txt'),           emit: summary_txt

    script:
    """
    results_assembly.py \\
        --mapping-table        ${mapping_table} \\
        --detection-table      ${detection_table} \\
        --sample-flags         ${sample_flags} \\
        --imputed-matrix       ${imputed_matrix} \\
        --imputation-mask      ${imputation_mask} \\
        --diff-abundance       ${diff_abundance} \\
        --enrichment-results   ${enrichment_results} \\
        --protein-term-mapping ${protein_term_mapping} \\
        --uniprot              ${uniprot_annotations} \\
        --pubmed               ${pubmed_cooccurrence} \\
        --disgenet             ${disgenet_associations} \\
        --dgidb                ${dgidb_interactions} \\
        --ctd                  ${ctd_interactions} \\
        --params               ${params_yml} \\
        --run-id               ${meta.run_id} \\
        --outdir               . \\
        --pipeline-version     ${workflow.manifest.version ?: 'dev'}
    """

    stub:
    """
    touch prosift_results.db
    touch ${meta.run_id}.diff_abundance_results.csv
    touch ${meta.run_id}.enrichment_results.csv
    touch ${meta.run_id}.significant_proteins.csv
    touch ${meta.run_id}.assembly_summary.txt
    """

}

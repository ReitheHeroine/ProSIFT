/*
 * ProSIFT -- workflows/prosift.nf
 * Top-level workflow. Imports and calls module processes in order.
 */

nextflow.enable.dsl = 2

// Module 01: Input Validation & ID Mapping
include { VALIDATE_INPUTS  } from '../modules/local/validate_inputs/main'
include { FILTER_PROTEINS  } from '../modules/local/filter_proteins/main'
include { UNIPROT_MAPPING  } from '../modules/local/uniprot_mapping/main'

// Module 01: Missingness Report (Process 4.10, advisory)
include { MISSINGNESS_REPORT } from '../modules/local/missingness_report/main'

// Module 02: Pre-Normalization QC/EDA
include { PRENORM_QC         } from '../modules/local/prenorm_qc/main'

// Module 03: Normalization & Imputation
include { NORMALIZE          } from '../modules/local/normalize/main'
include { IMPUTE             } from '../modules/local/impute/main'

// Module 04: Differential Abundance
include { DIFFERENTIAL_ABUNDANCE } from '../modules/local/differential_abundance/main'

// Module 05: Enrichment Analysis
include { ENRICHMENT             } from '../modules/local/enrichment/main'

workflow PROSIFT {

    // --- Build per-run input channel from samplesheet ---
    // Each emission: [meta, abundance_matrix, metadata, params_yml]
    // meta is a Map; run_id is the only key used by Module 01.
    //
    // ch_input is a queue channel (from Channel.fromPath) and can only be
    // consumed once per branch. multiMap is required here. Process output
    // channels (VALIDATE_INPUTS.out.*, FILTER_PROTEINS.out.*) are broadcast
    // in DSL2 and do not need multiMap even with multiple downstream subscribers.
    //
    // Branches:
    //   validate          -> VALIDATE_INPUTS
    //   filter_params     -> FILTER_PROTEINS (joined after validation)
    //   map_params        -> UNIPROT_MAPPING (joined after detection filter)
    //   prenorm_params    -> PRENORM_QC (joined after ID mapping)
    //   normalize_params  -> NORMALIZE (joined after QC)
    //   impute_params     -> IMPUTE (joined after NORMALIZE)
    Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            def meta = [run_id: row.run_id]
            [ meta,
              file(row.abundance, checkIfExists: true),
              file(row.metadata,  checkIfExists: true),
              file(row.params,    checkIfExists: true) ]
        }
        .multiMap { meta, abund, meta_csv, params_yml ->
            // validate branch: all four inputs for VALIDATE_INPUTS
            validate:         [ meta, abund, meta_csv, params_yml ]
            // filter_params: params_yml for the FILTER_PROTEINS join
            filter_params:    [ meta, params_yml ]
            // map_params: params_yml for the UNIPROT_MAPPING join
            map_params:       [ meta, params_yml ]
            // prenorm_params: params_yml for the PRENORM_QC join
            prenorm_params:   [ meta, params_yml ]
            // normalize_params: params_yml for the NORMALIZE join
            normalize_params: [ meta, params_yml ]
            // impute_params: params_yml for the IMPUTE join
            impute_params:    [ meta, params_yml ]
            // da_params: params_yml for the DIFFERENTIAL_ABUNDANCE join
            da_params:        [ meta, params_yml ]
            // enrich_params: params_yml for the ENRICHMENT join
            enrich_params:    [ meta, params_yml ]
        }
        .set { ch_input }

    // --- Module 01, Processes 4.1-4.3: Input validation ---
    VALIDATE_INPUTS(ch_input.validate)

    // --- Module 01, Process 4.4: Detection filter ---
    // Process output channels (VALIDATE_INPUTS.out.*) are broadcast in DSL2 --
    // multiple downstream joins can subscribe without multiMap.
    VALIDATE_INPUTS.out.matrix
        .join(VALIDATE_INPUTS.out.metadata)
        .join(ch_input.filter_params)
        .set { ch_filter_input }

    FILTER_PROTEINS(ch_filter_input)

    // --- Module 01, Process 4.10: Missingness report (advisory) ---
    // Joins: filter_table (FILTER_PROTEINS) + validated_matrix (pre-filter, VALIDATE_INPUTS)
    //        + validated_metadata (VALIDATE_INPUTS).
    FILTER_PROTEINS.out.filter_table
        .join(VALIDATE_INPUTS.out.matrix)
        .join(VALIDATE_INPUTS.out.metadata)
        .set { ch_missingness_input }

    MISSINGNESS_REPORT(ch_missingness_input)

    // --- Module 01, Processes 4.5-4.9: UniProt ID mapping ---
    FILTER_PROTEINS.out.matrix
        .join(ch_input.map_params)
        .set { ch_mapping_input }

    UNIPROT_MAPPING(ch_mapping_input)

    // --- Module 02: Pre-normalization QC/EDA ---
    // Joins: filtered_matrix (post-filter) + validated_metadata + id_mapping + params_yml
    FILTER_PROTEINS.out.matrix
        .join(VALIDATE_INPUTS.out.metadata)
        .join(UNIPROT_MAPPING.out.mapping_table)
        .join(ch_input.prenorm_params)
        .set { ch_prenorm_input }

    PRENORM_QC(ch_prenorm_input)

    // --- Module 03: Normalization ---
    // Joins: filtered_matrix (post-filter) + validated_metadata + params_yml
    // PRENORM_QC runs in parallel -- NORMALIZE does not depend on its outputs.
    FILTER_PROTEINS.out.matrix
        .join(VALIDATE_INPUTS.out.metadata)
        .join(ch_input.normalize_params)
        .set { ch_normalize_input }

    NORMALIZE(ch_normalize_input)

    // --- Module 03: Imputation ---
    // Joins: normalized_matrix (NORMALIZE) + validated_metadata (VALIDATE_INPUTS)
    //        + filter_table (FILTER_PROTEINS, for MNAR/MAR classification) + params_yml
    NORMALIZE.out.normalized_matrix
        .join(VALIDATE_INPUTS.out.metadata)
        .join(FILTER_PROTEINS.out.filter_table)
        .join(ch_input.impute_params)
        .set { ch_impute_input }

    IMPUTE(ch_impute_input)

    // --- Module 04: Differential Abundance ---
    // Joins: imputed_matrix (IMPUTE) + validated_metadata (VALIDATE_INPUTS)
    //        + mapping_table (UNIPROT_MAPPING) + params_yml
    IMPUTE.out.imputed_matrix
        .join(VALIDATE_INPUTS.out.metadata)
        .join(UNIPROT_MAPPING.out.mapping_table)
        .join(ch_input.da_params)
        .set { ch_da_input }

    DIFFERENTIAL_ABUNDANCE(ch_da_input)

    // --- Module 05: Enrichment Analysis ---
    // Joins: diff_abundance_results (DIFFERENTIAL_ABUNDANCE) + params_yml
    // GMT library paths are resolved from params.yml (not Nextflow file() staging)
    // so no additional path inputs are needed in the channel.
    DIFFERENTIAL_ABUNDANCE.out.results_table
        .join(ch_input.enrich_params)
        .set { ch_enrich_input }

    ENRICHMENT(ch_enrich_input)

}

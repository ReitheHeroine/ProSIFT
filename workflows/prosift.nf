/*
 * ProSIFT -- workflows/prosift.nf
 * Top-level workflow. Imports and calls module processes in order.
 */

nextflow.enable.dsl = 2

// Module 01: Input Validation & ID Mapping
include { VALIDATE_INPUTS  } from '../modules/local/validate_inputs/main'
include { FILTER_PROTEINS  } from '../modules/local/filter_proteins/main'
include { UNIPROT_MAPPING  } from '../modules/local/uniprot_mapping/main'

// Module 02: Pre-Normalization QC/EDA
include { PRENORM_QC       } from '../modules/local/prenorm_qc/main'

workflow PROSIFT {

    // --- Build per-run input channel from samplesheet ---
    // Each emission: [meta, abundance_matrix, metadata, params_yml]
    // meta is a Map; run_id is the only key used by Module 01.
    //
    // Four multiMap branches are needed because each downstream join
    // consumes a channel exactly once:
    //   validate       -> VALIDATE_INPUTS
    //   filter_params  -> FILTER_PROTEINS (joined after validation)
    //   map_params     -> UNIPROT_MAPPING (joined after detection filter)
    //   prenorm_params -> PRENORM_QC (joined after ID mapping)
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
            validate:       [ meta, abund, meta_csv, params_yml ]
            // filter_params: params_yml for the FILTER_PROTEINS join
            filter_params:  [ meta, params_yml ]
            // map_params: params_yml for the UNIPROT_MAPPING join
            map_params:     [ meta, params_yml ]
            // prenorm_params: params_yml for the PRENORM_QC join
            prenorm_params: [ meta, params_yml ]
        }
        .set { ch_input }

    // --- Module 01, Processes 4.1-4.3: Input validation ---
    VALIDATE_INPUTS(ch_input.validate)

    // VALIDATE_INPUTS.out.metadata is consumed by two downstream joins:
    // FILTER_PROTEINS and PRENORM_QC. multiMap creates two independent
    // branch channels, each receiving all items.
    VALIDATE_INPUTS.out.metadata
        .multiMap { meta, path ->
            filter:  [ meta, path ]
            prenorm: [ meta, path ]
        }
        .set { ch_validated_metadata }

    // --- Module 01, Process 4.4: Detection filter ---
    // Join validated matrix + validated metadata + params_yml into one channel.
    // .join() matches on the first element (meta) across all three sources.
    VALIDATE_INPUTS.out.matrix
        .join(ch_validated_metadata.filter)
        .join(ch_input.filter_params)
        .set { ch_filter_input }

    FILTER_PROTEINS(ch_filter_input)

    // --- Module 01, Processes 4.5-4.9: UniProt ID mapping ---
    FILTER_PROTEINS.out.matrix
        .join(ch_input.map_params)
        .set { ch_mapping_input }

    UNIPROT_MAPPING(ch_mapping_input)

    // --- Module 02: Pre-normalization QC/EDA ---
    // Joins: filtered_matrix (post-filter) + validated_metadata + id_mapping + params_yml
    FILTER_PROTEINS.out.matrix
        .join(ch_validated_metadata.prenorm)
        .join(UNIPROT_MAPPING.out.mapping_table)
        .join(ch_input.prenorm_params)
        .set { ch_prenorm_input }

    PRENORM_QC(ch_prenorm_input)

}

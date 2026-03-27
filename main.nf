#!/usr/bin/env nextflow

/*
 * ProSIFT -- PROtein Statistical Integration and Filtering Tool
 * Author: Reina Hastings, Blanco-Suárez Lab, SDSU
 *
 * Entry point. Validates params and calls the main workflow.
 */

nextflow.enable.dsl = 2

include { PROSIFT } from './workflows/prosift'

workflow {
    PROSIFT()
}

#!/usr/bin/env nextflow

/*
================================================================================
Coriell Institute for Medical Research

Contributors:
Anthony Pompetti <apompetti@coriell.org>
================================================================================
*/

/*
================================================================================
Enable Nextflow DSL2
================================================================================
*/
nextflow.enable.dsl=2

/*
================================================================================
Configurable variables for pipeline
================================================================================
*/

params.sample_table = false
params.input_type = false
params.genome = false
params.db = params.genome ? params.genomes[ params.genome ].db ?: false : false

params.star = false // specify whether to align using STAR aligner (default: false)

/*
================================================================================
Include modules to main pipeline
================================================================================
*/

/*
================================================================================
Include functions to main pipeline
================================================================================
*/

include { createInputChannel } from './functions/main.nf'

/*
================================================================================
Include subworkflows to main pipeline
================================================================================
*/

include { PREPROCESS } from './subworkflows/preprocess/main.nf'
include { ALIGN_QUANTIFY } from './subworkflows/align_quantify/main.nf'

/*
================================================================================
Workflow declaration
================================================================================
*/

workflow {

    // Ingest sample table to create input channel
    input_ch = createInputChannel(params.sample_table, params.input_type)

    // Run PREPROCESS subworkflow on input_ch (sample_table)
    PREPROCESS(input_ch)

    // Run ALIGN_QUANTIFY subworkflow on PREPROCESS.out.fastp
    ALIGN_QUANTIFY(PREPROCESS.out.fastp)

}
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

/*
================================================================================
Include modules to main pipeline
================================================================================
*/

include { PREPROCESS_READS } from '../../modules/preprocess/main.nf'
include { FASTP } from '../../modules/fastp/main.nf'
include { MULTIQC } from '../../modules/multiqc/main.nf'

/*
================================================================================
Include functions to main pipeline
================================================================================
*/


/*
================================================================================
Workflow declaration
================================================================================
*/

workflow PREPROCESS {

    take:
        input_ch

    main:
        // Preprocess the sample table to change the files listed inside. Concatenates any multilane files.
        PREPROCESS_READS(input_ch)

        // Trim and filter reads 
        FASTP(PREPROCESS_READS.out)

        // Run multiqc on the fastp output
        MULTIQC(
            FASTP.out.json.collect(),
            "multiqc/preprocess",
            "PREPROCESS Report",
            "multiqc_report"
        )

        emit:
            fastp = FASTP.out.reads
         
}
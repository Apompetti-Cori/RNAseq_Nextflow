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
params.multiqc_config = "${projectDir}/modules/multiqc/multiqc_config.yaml"
params.screenConfig = "${projectDir}/modules/fastq_screen/fastq_screen.conf"

/*
================================================================================
Include modules to main pipeline
================================================================================
*/

include { PREPROCESS_READS } from '../../modules/preprocess/main.nf'
include { FASTP } from '../../modules/fastp/main.nf'
include { FASTQ_SCREEN } from '../../modules/fastq_screen/main.nf'
include { MULTIQC } from '../../modules/multiqc/main.nf'

/*
================================================================================
Include functions to main pipeline
================================================================================
*/

/*
================================================================================
Include subworkflows to main pipeline
================================================================================
*/

include { SETUP } from '../../subworkflows/setup/main.nf'


/*
================================================================================
Workflow declaration
================================================================================
*/

workflow PREPROCESS {

    take:

    main:

        // Create an empty channel for multiqc input
        multiqc_ch = Channel.empty()

        // Run SETUP: Ingest sample table to create input channel
        SETUP()

        // Preprocess the sample table to change the files listed inside. Concatenates any multilane files.
        PREPROCESS_READS(SETUP.out.input_ch)

        // Trim and filter reads 
        FASTP(PREPROCESS_READS.out)
        multiqc_ch = multiqc_ch.mix(FASTP.out.json)

        // Screen reads for contamination if specified
        if(params.screen){

        FASTQ_SCREEN(
            PREPROCESS_READS.out,
            Channel.fromPath( "${params.screenConfig}" )
        )
        multiqc_ch = multiqc_ch.mix(FASTQ_SCREEN.out.meta_files) // add fastq screen meta files to multiqc channel

        }

        // Run multiqc on the multiqc channel
        MULTIQC(
            multiqc_ch.collect(),
            Channel.fromPath( "${params.multiqc_config}" ),
            "multiqc/preprocess",
            "PREPROCESS_Report",
            "multiqc_report"
        )

    emit:
        fastp = FASTP.out.reads
         
}
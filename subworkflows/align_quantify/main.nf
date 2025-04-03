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

include { SALMON_QUANT } from '../../modules/salmon/quant/main.nf'
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

workflow ALIGN_QUANTIFY {

    take:
        reads_ch

    main:

        // Quantify the reads using Salmon
        SALMON_QUANT(
            reads_ch,
            channel.fromPath( "${params.db}" )
        )

        MULTIQC(
            SALMON_QUANT.out.meta_files.collect(),
            Channel.fromPath( "${params.multiqc_config}" ),
            "multiqc/align_quantify",
            "ALIGN_QUANTIFY_Report",
            "multiqc_report"
        )

    emit:
        salmon = SALMON_QUANT.out.quants

}
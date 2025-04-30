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

/*
================================================================================
Include modules to main pipeline
================================================================================
*/

include { SALMON_QUANT } from '../../modules/salmon/quant/main.nf'
include { STAR } from '../../modules/star/main.nf'
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

        // Create an empty channel for multiqc input
        multiqc_ch = Channel.empty()

        // Quantify the reads using Salmon pseudo-aligner
        SALMON_QUANT(
            reads_ch,
            channel.fromPath( "${params.db}" )
        )
        multiqc_ch = multiqc_ch.mix(SALMON_QUANT.out.meta_files) // add salmon meta files to multiqc channel

        // Align and quantify using STAR aligner (only if use specifies for STAR alignments)
        if(params.star){

        STAR(
            reads_ch,
            channel.fromPath( "${params.db}" )
        )
        multiqc_ch = multiqc_ch.mix(STAR.out.meta_files) // add star meta files to multiqc channel

        }

        // Perform multiqc on meta files
        MULTIQC(
            multiqc_ch.collect(),
            Channel.fromPath( "${params.multiqc_config}" ),
            "multiqc/align_quantify",
            "ALIGN_QUANTIFY_Report",
            "multiqc_report"
        )

    emit:
        salmon = SALMON_QUANT.out.quants

}
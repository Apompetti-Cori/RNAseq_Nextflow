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
Configurable variables for module
================================================================================
*/

params.outdir = "./nfoutput"
params.pubdir = "multiqc"
params.fileprefix = "multiqc_report"
params.multiqc_config = "./modules/multiqc/multiqc_config.yaml"
params.multiqc_report_title = "MultiQC Report"

/*
================================================================================
Module declaration
================================================================================
*/


process MULTIQC {

    memory '8 GB'
    cpus 1

    conda 'bioconda::multiqc'

    publishDir "${params.outdir}/${params.pubdir}"

    input:
    path('multiqc_input/*')

    output:
    path("*.html")

    script:
    """
    multiqc multiqc_input/ \
        --config ${params.multiqc_config} \
        --title ${params.multiqc_report_title} \
        --filename ${params.fileprefix}.html
    """
}
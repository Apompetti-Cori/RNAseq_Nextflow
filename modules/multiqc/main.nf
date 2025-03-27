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
params.multiqc_config = "./modules/multiqc/multiqc_config.yaml"

/*
================================================================================
Module declaration
================================================================================
*/


process MULTIQC {

    memory '8 GB'
    cpus 1

    conda 'bioconda::multiqc'

    publishDir "${params.outdir}/${pubdir}"

    input:
    path('multiqc_input/*')
    val(pubdir)
    val(report_title)
    val(fileprefix)

    output:
    path("*.html")

    script:
    """
    multiqc multiqc_input/ \
        --config ${params.multiqc_config} \
        --title ${report_title} \
        --filename ${fileprefix}.html
    """
}
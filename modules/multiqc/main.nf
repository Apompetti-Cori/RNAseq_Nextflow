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

/*
================================================================================
Module declaration
================================================================================
*/


process MULTIQC {

    memory '8 GB'
    cpus 1

    tag { "${report_title}" }

    conda 'bioconda::multiqc'
    
    publishDir "${params.outdir}", mode: 'link',  saveAs: { "${pubdir}/${it}" }

    input:
    path('multiqc_input/*')
    path(multiqc_config)
    val(pubdir)
    val(report_title)
    val(fileprefix)

    output:
    path("*.html")

    script:
    """
    multiqc multiqc_input/ \
        --config ${multiqc_config} \
        --title ${report_title} \
        --filename ${fileprefix}.html
    """
}
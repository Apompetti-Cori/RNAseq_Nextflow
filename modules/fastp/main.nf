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
params.pubdir = "fastp"

/*
================================================================================
Module declaration
================================================================================
*/

process FASTP {

    maxForks 4
    memory '8 GB'
    cpus 4

    // Set batch name and sample id to tag
    tag { meta.batch == '' ? "${meta.id}" : "${meta.batch}_${meta.id}" }

    // Check batch and save output accordingly
    publishDir "${params.outdir}/${meta.id}",  saveAs: { meta.batch == '' ? "${it}/${params.pubdir}" : "${meta.batch}/${it}/${params.pubdir}" }, mode: 'link'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fq.gz"), emit: reads
    path("*.html"), emit: html
    path("*.json"), emit: json

    script:

    if ( meta.single_end ){

        """
        fastp \
        --in1 ${reads} \
        --out1 ${meta.id}.trimmed.fq.gz \
        -h ${meta.id}.fastp.html \
        -j ${meta.id}.fastp.json \
        -w ${task.cpus}
        """

    }
    else {

        """
        fastp \
        --in1 ${reads[1]} \
        --in2 ${reads[2]} \
        --out1 ${meta.id}.trimmed.1.gz \
        --out2 ${meta.id}.trimmed.2.gz \
        -h ${meta.id}.fastp.html \
        -j ${meta.id}.fastp.json \
        -w ${task.cpus}
        """

    }

}
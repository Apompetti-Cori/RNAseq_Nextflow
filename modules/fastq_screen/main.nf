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
params.pubdir = "fastq_screen"

/*
================================================================================
Module declaration
================================================================================
*/

process FASTQ_SCREEN {

    maxForks 4
    memory '8 GB'
    cpus 4

    // Set batch name and sample id to tag
    tag { meta.batch == '' ? "${meta.id}" : "${meta.batch}_${meta.id}" }

    conda 'bioconda::bbmap'

    // Check batch and save output accordingly (Don't save trimmed fastq files)
    publishDir "${params.outdir}", enabled: false, mode: 'link', saveAs: {
      filename ->
        if (filename.endsWith(".gz")) {
            return null
        } else {
            return meta.batch == '' ? "${meta.id}/${params.pubdir}/${filename}" : "${meta.batch}/${meta.id}/${params.pubdir}/${filename}"
        }
    }

    input:
    tuple val(meta), path(reads)
    each path(config)

    output:
    path("*.txt"), emit: meta_files

    script:
    """
    fastq_screen \
    --aligner 'bwa' \
    --threads ${task.cpus} \
    --conf ${config} \
    ${reads}
    """

}
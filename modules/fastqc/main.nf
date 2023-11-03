#!/usr/bin/env nextflow

/*
Coriell Institute for Medical Research
RNAseq Pipeline. Started November 2023.

Contributors:
Anthony Pompetti <apompetti@coriell.org>

Methodology adapted from:
Gennaro Calendo
*/

/*
Enable Nextflow DSL2
*/
nextflow.enable.dsl=2

/*
Define local params 
*/
params.outdir = "./results"
params.pubdir = "fastqc"

/*
Run fastqc on fastq files
*/
process FASTQC {
    maxForks 4
    memory '8 GB'
    cpus 2
    
    publishDir "${params.outdir}/${params.pubdir}", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    path("*.{html,zip}")

    script:
    """
    fastqc -t $task.cpus $reads
    """
}
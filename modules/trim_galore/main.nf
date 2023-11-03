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
params.pubdir = "trim_galore"

/*
Run trim_galore on each read stored within the reads_ch channel
*/
process TRIM_GALORE {
    maxForks 4
    memory '8 GB'
    cpus 4

    publishDir "${params.outdir}/trim_galore", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fq.gz"), emit: reads
    path("*trimming_report.txt"), emit: report

    script:
    def singleEnd = meta.single_end ? '' : '--paired'

    """
    trim_galore \
    ${singleEnd} \
    --length 35 \
    --quality 28 \
    --phred33 \
    --cores ${task.cpus} \
    ${reads}
    """
}
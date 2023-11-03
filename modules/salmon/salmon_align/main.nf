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
params.pubdir = "salmon_align"
params.db = false

/*
Run salmon on each read stored within the reads_ch channel
*/
process SALMON_ALIGN {
    maxForks 4
    memory '8 GB'
    cpus 4

    publishDir "${params.outdir}/salmon_align", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("_quants"), emit: reads
    path("meta_info.json"), emit: report

    script:
    def singleEnd = meta.single_end ? '' : '--paired'

    """
    salmon quant -i ${params.db} \
                 -l 'A' \
                 -1 ${reads[0]} \
                 -2 ${reads[1]} \
                 -o ${meta.id}_quants \
                 --threads ${task.cpus} \
                 --validateMappings \
                 --gcBias \
                 --seqBias \
                 --numBootstraps ${params.salmon_boots};
    """
}
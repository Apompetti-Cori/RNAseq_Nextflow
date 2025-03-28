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

params.outdir = "./results"
params.pubdir = "salmon_align"
params.db = false
params.salmon_boots = 30

/*
================================================================================
Module declaration
================================================================================
*/

process SALMON_ALIGN {
    maxForks 4
    memory '8 GB'
    cpus 4

    // Check batch and save output accordingly (only save quant files)
    publishDir "${params.outdir}", mode: 'link',  saveAs: {
      filename -> filename.endsWith(".gz") ?
      { meta.batch == '' ? "${it}/${meta.id}/${params.pubdir}" : "${it}/${meta.batch}/${meta.id}/${params.pubdir}" } :
      null
    }

    input:
    tuple val(meta), path(reads)

    output:
    path("*_quants.tar.gz"), emit: quants
    path("*_quants/aux_info/meta_info.json"), emit: meta_info
    path("*_quants/libParams/flenDist.txt"), emit: flenDist
    
    script:
    if(meta.single_end){
      """
      salmon quant -i ${params.db} \
                   -l 'A' \
                   -r ${reads} \
                   -o ${meta.id}_quants \
                   --threads ${task.cpus} \
                   --validateMappings \
                   --gcBias \
                   --seqBias \
                   --numBootstraps ${params.salmon_boots};

      tar -zcvf ${meta.id}_quants.tar.gz ${meta.id}_quants
      """
    }
    else{
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

      tar -zcvf ${meta.id}_quants.tar.gz ${meta.id}_quants
      """
    }
}
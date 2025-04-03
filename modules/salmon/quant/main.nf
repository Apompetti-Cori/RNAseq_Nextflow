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
params.pubdir = "salmon"
params.salmon_boots = 30

/*
================================================================================
Module declaration
================================================================================
*/

process SALMON_QUANT {
    maxForks 4
    memory '8 GB'
    cpus 8

    // Set batch name and sample id to tag
    tag { meta.batch == '' ? "${meta.id}" : "${meta.batch}_${meta.id}" }

    // Check batch and save output accordingly (only save quant files)
    publishDir "${params.outdir}", mode: 'link', saveAs: { 
      filename ->
        if (filename.endsWith(".gz")) {
            return meta.batch == '' ? "${meta.id}/${params.pubdir}/${filename}" : "${meta.batch}/${meta.id}/${params.pubdir}/${filename}"
        } else {
            return null
        }
    }

    input:
    tuple val(meta), path(reads)
    each path(db)

    output:
    tuple val(meta), path("*_quants.tar.gz"), emit: quants
    path("*_quants"), emit: meta_files
    
    script:
    if(meta.single_end){
      """
      salmon quant -i ${db}/salmon_idx \
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
      salmon quant -i ${db}/salmon_idx \
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
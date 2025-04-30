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
params.pubdir = "star"

/*
================================================================================
Module declaration
================================================================================
*/

process STAR {
    maxForks 4
    memory '8 GB'
    cpus 16

    // Set batch name and sample id to tag
    tag { meta.batch == '' ? "${meta.id}" : "${meta.batch}_${meta.id}" }

    // Check batch and save output accordingly (only save bam files)
    publishDir "${params.outdir}", mode: 'link', pattern: "*.bam{,.bai}", saveAs: { 
      filename -> {
        return meta.batch == '' ? "${meta.id}/${params.pubdir}/${filename}" : "${meta.batch}/${meta.id}/${params.pubdir}/${filename}"
      }
    }

    input:
    tuple val(meta), path(reads)
    each path(db)

    output:
    tuple val(meta), path("*.bam*"), emit: bam
    tuple path("*Log.final.out"), path("*ReadsPerGene.out.tab"), emit: meta_files

    script:
    """
    STAR \
        --readFilesIn ${reads} \
        --genomeDir ${db}/gencode.v38.STAR_idx \
        --readFilesCommand zcat \
        --runThreadN ${task.cpus} \
        --twopassMode Basic \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverLmax 0.1 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --outFilterType BySJout \
        --outFilterScoreMinOverLread 0.33 \
        --outFilterMatchNminOverLread 0.33 \
        --limitSjdbInsertNsj 1200000 \
        --outFileNamePrefix ${meta.id} \
        --outSAMstrandField intronMotif \
        --outFilterIntronMotifs None \
        --alignSoftClipAtReferenceEnds Yes \
        --quantMode GeneCounts \
        --outSAMtype BAM SortedByCoordinate \
        --genomeLoad NoSharedMemory \
        --chimSegmentMin 15 \
        --chimJunctionOverhangMin 15 \
        --chimOutType Junctions \
        --chimOutJunctionFormat 1 \
        --chimMainSegmentMultNmax 1 \
        --outSAMattributes NH HI AS nM NM ch

    samtools index ${meta.id}Aligned.sortedByCoord.out.bam
    """
}
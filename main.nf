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

//Configurable variables for pipeline
params.fastq_folder = "${workflow.projectDir}/fastq"
params.reads = "${params.fastq_folder}/*{_L00,}{1,2,3,4}{_R,_}{1,2}*.{fastq,fq}.gz"
params.singleEnd = false
params.multiLane = false
params.amplifytargets = false
params.genome = false
params.db = params.genomes ? params.genomes[ params.genome ].db ?:false : false

//Include modules to main pipeline
include { FASTQC as PRETRIM_FASTQC } from './modules/fastqc/main' addParams(pubdir: 'pretrim_fastqc')
include { TRIM_GALORE } from './modules/trim_galore/main'
include { FASTQC as POSTTRIM_FASTQC } from './modules/fastqc/main' addParams(pubdir: 'posttrim_fastqc')
include { SALMON_ALIGN } from './modules/salmon/salmon_align/main'
include { MULTIQC } from './modules/multiqc/main'

//Provide sample table in csv format to have pipeline process samples via sample table
params.sample_table = "${workflow.projectDir}/sample_table.csv"
if ( params.sample_table ){
    // Channel for the samplesheet
    ch_samplesheet = Channel.fromPath(params.sample_table)
    
    //ch_samplesheet.view()
    
    // Parse it line by line
    reads_ch = ch_samplesheet.splitCsv(header:true).map {

        // This is the read1 and read2 entry
        r1_L1 = it['r1_L1']
        r1_L2 = it['r1_L2']
        r1_L3 = it['r1_L3']
        r1_L4 = it['r1_L4']
        r2_L1 = it['r2_L1']
        r2_L2 = it['r2_L2']
        r2_L3 = it['r2_L3']
        r2_L4 = it['r2_L4']

        // Detect wiether single-end or paired-end
        is_singleEnd = r2_L1.toString()=='' ? true : false
        is_multiLane = r1_L2.toString()=='' ? false : true
        
        // The "meta" map, which is a Nextflow/Groovy map with id (the sample name) and a single_end logical entry
        meta = [id: it['sample'], single_end: is_singleEnd, multi_lane: is_multiLane]
        
        // We return a nested map, the first entry is the meta map, the second one is the read(s)
        if ( is_singleEnd ){
            if ( r1_L4.toString()!='' ){
                [meta, [r1_L1, r1_L2, r1_L3, r1_L4], []]
            }
            else if ( r1_L3.toString()!='' ){
                [meta, [r1_L1, r1_L2, r1_L3], []]
            }
            else if ( r1_L2.toString()!='' ){
                [meta, [r1_L1, r1_L2], []]
            }
            else {
                [meta, [r1_L1], []]
            }
        }
        else{
            if (r2_L4.toString()!=''){
                [meta, [r1_L1, r1_L2, r1_L3, r1_L4], [r2_L1, r2_L2, r2_L3, r2_L4]]
            }
            else if ( r2_L3.toString()!='' ){
                [meta, [r1_L1, r1_L2, r1_L3], [r2_L1, r2_L2, r2_L3]]
            }
            else if ( r2_L2.toString()!='' ){
                [meta, [r1_L1, r1_L2], [r2_L1, r2_L2]]
            }
            else {
                [meta, [r1_L1], [r2_L1]]
            }
        }
    }
}
else {
    Channel
    .fromFilePairs(params.reads, size: params.singleEnd ? (params.multiLane ?: 1) : 2)
    .ifEmpty {exit 1, "Cannot find any reads matching ${params.reads}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --singleEnd on the command line."}
    .set{ reads_ch }
}


process PREPROCESS {
    input:
    tuple val(meta), path(reads_1), path(reads_2)

    output:
    tuple val(meta), path("*{_mL,_sL,_sL_1,_sL_2,_mL_1,_mL_2}.fq.gz")

    script:
    if ( meta.multi_lane ){
        if (meta.single_end){
            """
            cat ${reads_1} > ${meta.id}_mL.fq.gz
            """
        }
        else {
            """
            cat ${reads_1} > ${meta.id}_mL_1.fq.gz
            cat ${reads_2} > ${meta.id}_mL_2.fq.gz
            """
        }
    }
    else {
        if (meta.single_end){
            """
            mv ${reads_1} ${meta.id}_sL.fq.gz
            """
        }
        else {
            """
            mv ${reads_1} ${meta.id}_sL_1.fq.gz
            mv ${reads_2} ${meta.id}_sL_2.fq.gz
            """
        }
    }
}


workflow {
    //reads_ch.view()
    PREPROCESS(reads_ch)
    reads_ch = PREPROCESS.out
    //reads_ch.view()
    //Run fastqc on raw reads
    PRETRIM_FASTQC(reads_ch)
    //Run trim_galore on raw reads
    TRIM_GALORE(reads_ch)

    //Run fastqc on trimmed reads, specifies trim_galore[0] because second input channel is not need for this process
    POSTTRIM_FASTQC(TRIM_GALORE.out.reads)
    
    //Run salmon_align on trimmed reads
    //State Dependency: Wait until TRIM_GALORE is done to run bismark align
    //State Dependency: Wait until POSTTRIM_FASTQC is done to run bismark align
    state = POSTTRIM_FASTQC.out.collect()
    SALMON_ALIGN(TRIM_GALORE.out.reads.collect(flat: false).flatMap(), state)

    //Run multiqc on pretrim fastqc output, trim_galore trimming report, posttrim fastqc output, bismark conversion output
    MULTIQC(PRETRIM_FASTQC.out.collect().combine(POSTTRIM_FASTQC.out.collect()).combine(TRIM_GALORE.out.report.collect()).combine(BISMARK_ALIGN.out.report.collect()).combine(CONV_STATS_CREATE.out.report.collect()))
}


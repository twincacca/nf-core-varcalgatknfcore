/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTP } from '../modules/nf-core/fastp/main'
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_varcalgatknfcore_pipeline'
include { BWA_INDEX } from '../modules/nf-core/bwa/index/main'
include { BWA_MEM } from '../modules/nf-core/bwa/mem/main'
include { PICARD_ADDORREPLACEREADGROUPS } from '../modules/nf-core/picard/addorreplacereadgroups/main'
include { PICARD_SORTSAM } from '../modules/nf-core/picard/sortsam/main'
include { PICARD_MARKDUPLICATES } from '../modules/nf-core/picard/markduplicates/main'
include { GATK4_BASERECALIBRATOR } from '../modules/nf-core/gatk4/baserecalibrator/main'
include { SAMTOOLS_INDEX } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_DICT } from '../modules/nf-core/samtools/dict/main'
include { SAMTOOLS_FAIDX } from '../modules/nf-core/samtools/faidx/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VARCALGATKNFCORE {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()



    //
    // MODULE: FASTP
    //
    // Create a new channel (ch_adapters) using a ternary expression to use either a supplied fasta file or [].
    ch_adapters = params.adapters ? params.adapters : []

    FASTP (
        ch_samplesheet,
        ch_adapters,
        params.discard_trimmed_pass,  // <--------- need to add this one also! outdated tutorial
        params.save_trimmed_fail,
        params.save_merged
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]}.ifEmpty([]))
    ch_versions      = ch_versions.mix(FASTP.out.versions.first())


    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )


    ch_genome_fasta = Channel.fromPath(params.fasta).map { it -> [[id:it[0].simpleName], it] }.collect()
    //
    // MODULE: BWA_INDEX
    //
    BWA_INDEX (
        ch_genome_fasta // tuple val(meta), path(fasta)
    )

    //
    // MODULE: BWA_MEM
    //
    BWA_MEM (
    ch_samplesheet, //tuple val(meta) , path(reads)
    BWA_INDEX.out.index, //tuple val(meta2), path(index)
    ch_genome_fasta, //tuple val(meta3), path(fasta)
    params.sort_bam
    )

    // //
    // // MODULE: PICARD_ADDORREPLACEREADGROUPS
    // //
    PICARD_ADDORREPLACEREADGROUPS (
    BWA_MEM.out.bam, // tuple val(meta), path(reads)
    ch_genome_fasta, // tuple val(meta2), path(fasta)
    BWA_INDEX.out.index // tuple val(meta3), path(fasta_index)
    )

    // //
    // // MODULE: PICARD_SORTSAM
    // //
    PICARD_SORTSAM (
    PICARD_ADDORREPLACEREADGROUPS.out.bam, // tuple val(meta), path(bam)
    params.sort_order // val sort_order
    )

    // //
    // // MODULE: PICARD_MARKDUPLICATES
    // //
    PICARD_MARKDUPLICATES (
    PICARD_SORTSAM.out.bam, // tuple val(meta), path(reads)
    ch_genome_fasta, // tuple val(meta2), path(fasta)
    BWA_INDEX.out.index // tuple val(meta3), path(fai)
    )

    // // //
    // // // MODULE: SAMTOOLS_INDEX
    // // //
    // SAMTOOLS_INDEX (
    // ch_genome_fasta // tuple val(meta), path(input)
    // )
    // // SAMTOOLS_INDEX.out.fai[1].view() // vorrei prendere il 2nd element ma si deve fare cosi:
    // fasta_fai = SAMTOOLS_INDEX.out.fai.subscribe { result ->
    // // Assuming result is a list or tuple, extract the second element
    // def second_argument = result[1]
    // println(second_argument)
    // }



    // //
    // // MODULE: SAMTOOLS_DICT
    // //
    SAMTOOLS_FAIDX (
    ch_genome_fasta // tuple val(meta), path(input)
    )


    // //
    // // MODULE: SAMTOOLS_DICT
    // //
    SAMTOOLS_DICT (
    ch_genome_fasta // tuple val(meta), path(input)
    )


    // transform info from mark dup into something that fits GATK4_BASERECALIBRATOR 
    picard_mark_dup_bam = PICARD_MARKDUPLICATES.out.bam
    picard_mark_dup_bai = PICARD_MARKDUPLICATES.out.bai
    picard_mark_dup_metrics = PICARD_MARKDUPLICATES.out.metrics
    
    tup_meta_bam_bai_int = (picard_mark_dup_bam.join(picard_mark_dup_bai)).join(picard_mark_dup_metrics)
    .map{ meta, bam, bai, met -> [ meta, bam, bai, [] ] } // [] = no intervals provided
    //tup_meta_bam_bai_int.view()
    //> [[id:WT_REP1, single_end:false], /home/antoinebuetti/Desktop/work/varCalling_NextFlow_nfcore/nf-core-varcalgatknfcore/work/da/08b2b9033a88d112e10107ebb488fa/WT_REP1.dup.bam, /home/antoinebuetti/Desktop/work/varCalling_NextFlow_nfcore/nf-core-varcalgatknfcore/work/da/08b2b9033a88d112e10107ebb488fa/WT_REP1.dup.bai, []]

    // //
    // // MODULE: GATK4_BASERECALIBRATOR
    // //
    GATK4_BASERECALIBRATOR (
    tup_meta_bam_bai_int, // tuple val(meta), path(input), path(input_index), path(intervals)
    params.fasta, // path  fasta
    SAMTOOLS_FAIDX.out.fai.map{ meta, f -> [f] }, // path  fai // remove meta with map
    SAMTOOLS_DICT.out.dict.map{ meta, f -> [f] }, // path  dict  // remove meta with map
    params.vcf, // path  known_sites // adesso li presi local!!!
    params.vcf_tbi // path  known_sites_tbi
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

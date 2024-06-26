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
// include { GATK4SPARK_BASERECALIBRATOR } from '../modules/nf-core/gatk4spark/baserecalibrator/main' 
include { SAMTOOLS_INDEX } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_DICT } from '../modules/nf-core/samtools/dict/main'
include { SAMTOOLS_FAIDX } from '../modules/nf-core/samtools/faidx/main'
include { GATK4_APPLYBQSR } from '../modules/nf-core/gatk4/applybqsr/main'
// include { GATK4SPARK_APPLYBQSR } from '../modules/nf-core/gatk4spark/applybqsr/main'                                                                                                         
include { GATK4_MUTECT2 } from '../modules/nf-core/gatk4/mutect2/main'
include { GATK4_FILTERMUTECTCALLS } from '../modules/nf-core/gatk4/filtermutectcalls'
include { GATK4_HAPLOTYPECALLER } from '../modules/nf-core/gatk4/haplotypecaller/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_forhaplotypecaller } from '../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_forbasecalibrator } from '../modules/nf-core/bcftools/index/main'

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

    //
    // MODULE: BWA_INDEX
    //
    ch_genome_fasta = Channel.fromPath(params.fasta).map { it -> [[id:it[0].simpleName], it] }.collect()
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



    // //
    // // MODULE: BCFTOOLS_INDEX
    // //
    ch_vcf_forbasecalibrator = Channel.fromPath(params.vcf_forbasecalibrator).map { it -> [[id:it[0].simpleName], it] }.collect()
    // ch_vcf_forbasecalibrator.view() 
    //> [['id':'nf-core'], /nf-core/test-datasets/raredisease/reference/dbsnp_-138-.vcf.gz]
    BCFTOOLS_INDEX_forbasecalibrator (
    ch_vcf_forbasecalibrator // tuple val(meta), path(vcf)
    )

    // //
    // // MODULE: GATK4_BASERECALIBRATOR
    // //
    // transform info from mark dup into something that fits GATK4_BASERECALIBRATOR 
    picard_mark_dup_bam = PICARD_MARKDUPLICATES.out.bam
    picard_mark_dup_bai = PICARD_MARKDUPLICATES.out.bai
    picard_mark_dup_metrics = PICARD_MARKDUPLICATES.out.metrics
    
    tup_meta_bam_bai_int = (picard_mark_dup_bam.join(picard_mark_dup_bai)).join(picard_mark_dup_metrics)
    .map{ meta, bam, bai, met -> [ meta, bam, bai, [] ] } // [] = no intervals provided
    //tup_meta_bam_bai_int.view()
    //> [[id:WT_REP1, single_end:false], /home/antoinebuetti/Desktop/work/varCalling_NextFlow_nfcore/nf-core-varcalgatknfcore/work/da/08b2b9033a88d112e10107ebb488fa/WT_REP1.dup.bam, /home/antoinebuetti/Desktop/work/varCalling_NextFlow_nfcore/nf-core-varcalgatknfcore/work/da/08b2b9033a88d112e10107ebb488fa/WT_REP1.dup.bai, []]

    GATK4_BASERECALIBRATOR (
    tup_meta_bam_bai_int, // tuple val(meta), path(input), path(input_index), path(intervals)
    params.fasta, // path  fasta
    SAMTOOLS_FAIDX.out.fai.map{ meta, f -> [f] }, // path  fai
    SAMTOOLS_DICT.out.dict.map{ meta, f -> [f] }, // path  dict 
    params.vcf_forbasecalibrator, // path  known_sites 
    BCFTOOLS_INDEX_forbasecalibrator.out.tbi.map{ meta, f -> [f] } // path  known_sites_tbi
    // params.vcf_forbasecalibrator_tbi
    )

    // //
    // // MODULE: GATK4_APPLYBQSR
    // //
    tup_meta_bam_bai_tab_int = tup_meta_bam_bai_int.join(GATK4_BASERECALIBRATOR.out.table)
    .map{ meta, bam, bai, interv, tab -> [ meta, bam, bai, tab, interv ] } 
    //.view()
    // //> [[id:WT_REP1, single_end:false], /home/antoinebuetti/Desktop/work/varCalling_NextFlow_nfcore/nf-core-varcalgatknfcore/work/c6/1184c8d59299ea09ec468fa6fe252f/WT_REP1.dup.bam, /home/antoinebuetti/Desktop/work/varCalling_NextFlow_nfcore/nf-core-varcalgatknfcore/work/c6/1184c8d59299ea09ec468fa6fe252f/WT_REP1.dup.bai, /home/antoinebuetti/Desktop/work/varCalling_NextFlow_nfcore/nf-core-varcalgatknfcore/work/8a/ef7e750ec4ab147c0f6ae129b3e0e9/WT_REP1.table, []]

    GATK4_APPLYBQSR (
    tup_meta_bam_bai_tab_int, // tuple val(meta), path(input), path(input_index), path(bqsr_table), path(intervals)
    params.fasta, // path  fasta
    SAMTOOLS_FAIDX.out.fai.map{ meta, f -> [f] }, // path  fai 
    SAMTOOLS_DICT.out.dict.map{ meta, f -> [f] }, // path  dict 
    )

    // //
    // // MODULE: GATK4_MUTECT2
    // //
    ch_vcf_germline_resource = params.vcf_germline_resource ? params.vcf_germline_resource : []
    ch_vcf_germline_resource_tbi = params.vcf_germline_resource_tbi ? params.vcf_germline_resource_tbi : []
    ch_vcf_panel_of_normals = params.vcf_panel_of_normals ? params.vcf_panel_of_normals : []
    ch_vcf_panel_of_normals_tbi = params.vcf_panel_of_normals_tbi ? params.vcf_panel_of_normals_tbi : []

    GATK4_MUTECT2 (
    tup_meta_bam_bai_int, // tuple val(meta), path(input), path(input_index), path(intervals)
    ch_genome_fasta, // tuple val(meta2), path(fasta)
    SAMTOOLS_FAIDX.out.fai, // tuple val(meta3), path(fai)
    SAMTOOLS_DICT.out.dict, // tuple val(meta4), path(dict)
    ch_vcf_germline_resource, // path(germline_resource) .vcf.gz
    ch_vcf_germline_resource_tbi, // path(germline_resource_tbi) .vcf.gz.tbi // VCF file of sites observed in normal.
    ch_vcf_panel_of_normals, // path(panel_of_normals) .vcf.gz // Population vcf of germline sequencing containing allele fractions.
    ch_vcf_panel_of_normals_tbi // path(panel_of_normals_tbi) .vcf.gz.tbi
    )


    // //
    // // MODULE: GATK4_FILTERMUTECTCALLS
    // //
    tup_meta_vcf_vcftbi = GATK4_MUTECT2.out.vcf.join(GATK4_MUTECT2.out.tbi, failOnMismatch:true, failOnDuplicate:true)
    //tup_meta_vcf_vcftbi.view()
    //> [[id:WT_REP1, single_end:false], .../WT_REP1.vcf.gz, .../WT_REP1.vcf.gz.tbi]
    tup_meta_vcf_vcftbi_stats = tup_meta_vcf_vcftbi.join(GATK4_MUTECT2.out.stats, failOnMismatch:true, failOnDuplicate:true)
    //tup_meta_vcf_vcftbi_stats.view()
    //> [[id:WT_REP1, single_end:false], .../WT_REP1.vcf.gz, .../WT_REP1.vcf.gz.tbi, .../WT_REP1.vcf.gz.stats]
    tup_meta_vcf_vcftbi_stats_oribias_segm_table_est    = tup_meta_vcf_vcftbi_stats.map {
                        meta, vcf, tbi, stats ->
                        return [meta, vcf, tbi, stats, [], [], [], []]
                    }
    //tup_meta_vcf_vcftbi_stats_oribias_segm_table_est.view()
    //> [[id:WT_REP1, single_end:false], .../WT_REP1.vcf.gz, .../WT_REP1.vcf.gz.tbi, .../WT_REP1.vcf.gz.stats, [], [], [], []]

    GATK4_FILTERMUTECTCALLS (
    tup_meta_vcf_vcftbi_stats_oribias_segm_table_est,// tuple val(meta), path(vcf), path(vcf_tbi), path(stats), path(orientationbias), path(segmentation), path(table), val(estimate)
    ch_genome_fasta, // tuple val(meta2), path(fasta)
    SAMTOOLS_FAIDX.out.fai, // tuple val(meta3), path(fai)
    SAMTOOLS_DICT.out.dict, // tuple val(meta4), path(dict)
    )

    // //
    // // MODULE: BCFTOOLS_INDEX
    // //
    ch_vcf_forhaplotypecaller = Channel.fromPath(params.vcf_forhaplotypecaller).map { it -> [[id:it[0].simpleName], it] }.collect()
    // ch_vcf_forhaplotypecaller.view() 
    //> [['id':'nf-core'], /nf-core/test-datasets/raredisease/reference/dbsnp_-138-.vcf.gz]
    // ch_vcf_forhaplotypecaller_tbi = Channel.fromPath(params.vcf_forhaplotypecaller_tbi).map { it -> [[id:it[0].simpleName], it] }.collect()

    BCFTOOLS_INDEX_forhaplotypecaller (
    ch_vcf_forhaplotypecaller // tuple val(meta), path(vcf)
    )

    // //
    // // MODULE: GATK4_HAPLOTYPECALLER
    // //
    tup_meta_bam_bai_int_dragstr = tup_meta_bam_bai_int.map{ meta, bam, bai, interv -> [ meta, bam, bai, interv, [] ] } 
    // tup_meta_bam_bai_int_dragstr.view()
    //> [[id:hugelymodelbat, single_end:false], .../hugelymodelbat.dup.bam,.../hugelymodelbat.dup.bai, [], []]

    GATK4_HAPLOTYPECALLER (
    tup_meta_bam_bai_int_dragstr, // tuple val(meta),  path(input), path(input_index), path(intervals), path(dragstr_model)
    ch_genome_fasta, // tuple val(meta2), path(fasta)
    SAMTOOLS_FAIDX.out.fai, // tuple val(meta3), path(fai)
    SAMTOOLS_DICT.out.dict, // tuple val(meta4), path(dict)
    ch_vcf_forhaplotypecaller, // tuple val(meta5), path(dbsnp)
    BCFTOOLS_INDEX_forhaplotypecaller.out.tbi // tuple val(meta6), path(dbsnp_tbi)
    // ch_vcf_forhaplotypecaller_tbi
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

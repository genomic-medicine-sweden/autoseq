/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap                                    } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                                } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                              } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                              } from '../subworkflows/local/utils_nfcore_autoseq_pipeline'


include { FASTQC                                              } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                             } from '../modules/nf-core/multiqc/main'
include { FASTP                                               } from '../modules/nf-core/fastp/main'
include { CAT_FASTQ                                           } from '../modules/nf-core/cat/fastq/main'
include { SAMTOOLS_INDEX                                      } from '../modules/nf-core/samtools/index/main'

include { ALIGNMENT                                           } from '../subworkflows/local/fastq_align_bwamem2/main.nf'
include { FASTQ_CREATE_UMI_CONSENSUS_FGBIO as UMI_PROCESSING  } from '../subworkflows/nf-core/fastq_create_umi_consensus_fgbio/main'
include { BAM_QC_PICARD_SAMTOOLS  as BAM_QC                   } from '../subworkflows/local/bam_qc_picard_samtools/main.nf'
include { CALL_SOMATIC_SNVS                                   } from '../subworkflows/local/bam_call_somatic_snvs/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow AUTOSEQ {

    take:
    ch_samplesheet     // channel: samplesheet read in from --input
    ch_genome_fasta
    ch_genome_fai
    ch_dict
    ch_bwamem2_index
    ch_targets_bed
    ch_targets_bed_gz
    ch_interval_list
    ch_interval_list_slopped20
    ch_jumble_ref

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )

    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: Run FastP
    //
    FASTP (
        ch_samplesheet,
        [], // adapter_fasta: not used in this pipeline
        params.discard_trimmed_pass,
        params.save_trimmed_fail,
        params.save_merged
    )

    ch_versions = ch_versions.mix(FASTP.out.versions.first())
    ch_input_reads = FASTP.out.reads

    //
    // MODULE: ALIGNMENT
    //

    if (params.umi_structure) {

        ch_input_reads
            .map { meta, reads ->
                def id = "${meta.case_id}.${meta.sample_name}".toString()
                meta   = meta + [id: id]
                return [meta.sample_name, [meta , reads]]
            }
            .groupTuple()
            .map { sample_name, grouped_reads ->
                def metas = grouped_reads.collect{it[0]}
                def files = grouped_reads.collect{it[1]}.flatten()
                return [metas[0], files]
            }
            .set { ch_input_reads }

        CAT_FASTQ (
            ch_input_reads
        )


        UMI_PROCESSING(
            CAT_FASTQ.out.reads,
            ch_genome_fasta,
            ch_bwamem2_index,
            ch_dict,
            "paired",
            "bwa-mem2",
            params.duplex,
            params.min_reads,
            params.min_baseq,
            params.max_base_error_rate
        )

        SAMTOOLS_INDEX(
            UMI_PROCESSING.out.mappedconsensusbam
        )

        ch_aligned_bam = UMI_PROCESSING.out.mappedconsensusbam
            .join(SAMTOOLS_INDEX.out.bai)

        ch_versions = ch_versions.mix(UMI_PROCESSING.out.versions.first())

    } else {

        ALIGNMENT(
            ch_input_reads,
            ch_genome_fasta,
            ch_genome_fai,
            ch_bwamem2_index
        )

        ch_multiqc_files = ch_multiqc_files.mix(ALIGNMENT.out.dedup_metrics.collect{it[1]}.ifEmpty([]))
        ch_versions = ch_versions.mix(ALIGNMENT.out.versions.first())
        ch_aligned_bam = ALIGNMENT.out.dedup_bam
            .join(ALIGNMENT.out.dedup_bai)
    }

    // Module: QC

    BAM_QC(
        ch_aligned_bam,
        ch_genome_fasta,
        ch_genome_fai,
        ch_dict,
        ch_interval_list_slopped20
    )

    // Module: SOMATIC SNV CALLING

    // Branch samples by tumor/normal
    ch_aligned_bam
        .branch { meta, bam, bai ->
            tumor: meta.sample_type == "tumor"
            normal: meta.sample_type == "normal"
        }
        .set { samples }

    // Prepare tumor channel with case_id as key
    tumor_ch = samples.tumor
        .map { meta, bam, bai ->
            [meta.case_id, meta, bam, bai]
        }

    // Prepare normal channel with case_id as key
    normal_ch = samples.normal
        .map { meta, bam, bai ->
            [meta.case_id, meta, bam, bai]
        }

    // Join tumor and normal by case_id and create somatic calling format
    ch_input_paired = tumor_ch
        .join(normal_ch, by: 0)
        .combine(ch_interval_list_slopped20)
        .map { case_id, tumor_meta, tumor_bam, tumor_bai, normal_meta, normal_bam, normal_bai, meta_intervals, intervals_file ->
            // Create comprehensive meta map
            def meta = [
                id: case_id,
                case_id: case_id,
                tumor_id: tumor_meta.id,
                normal_id: normal_meta.id,
                tumor_sample: tumor_meta.sample_name,
                normal_sample: normal_meta.sample_name
            ]

            // Return in desired format: [meta, [tumor.bam, normal.bam], [tumor.bai, normal.bai], intervals]
            [meta, [tumor_bam, normal_bam], [tumor_bai, normal_bai], intervals_file]
        }

    CALL_SOMATIC_SNVS (
        ch_input_paired,
        ch_genome_fasta,
        ch_genome_fai,
        ch_dict,
        [],
        [],
        [],
        [],
        ch_interval_list_slopped20.collect{it[1]}
    )

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'autoseq_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BAM_QC.out.multiple_metrics.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BAM_QC.out.hs_metrics.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BAM_QC.out.flagstat.collect{it[1]}.ifEmpty([]))



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
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

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
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

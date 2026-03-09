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

include { READ_ALIGNMENT                                      } from '../subworkflows/local/fastq_align_bwamem2/main.nf'
include { FASTQ_CREATE_UMI_CONSENSUS_FGBIO as UMI_PROCESSING  } from '../subworkflows/nf-core/fastq_create_umi_consensus_fgbio/main'
include { BAM_QC_PICARD_SAMTOOLS  as ALIGNMENT_QC             } from '../subworkflows/local/bam_qc_picard_samtools/main.nf'
include { SOMATIC_SNV_CALLING                                 } from '../subworkflows/local/call_somatic_snvs/main.nf'
include { CNV_CALLING                                         } from '../subworkflows/local/call_cnvs/main.nf'
include { SVS_CALLING                                         } from '../subworkflows/local/call_svs/main.nf'
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
    ch_interval_list_slopped20
    ch_jumble_ref
    ch_sage_known_hotspots_somatic
    ch_sage_highconf_regions
    ch_sage_pon
    ch_ensembl_vep_cache
    ch_ensembl_data_resources
    ch_curation_ann
    ch_germline_resource
    ch_germline_resource_tbi
    ch_genome_gridss_index
    ch_pon_breakends    // channel: panel of normals breakends for SV calling
    ch_pon_breakpoints  // channel: panel of normals breakpoints for SV calling
    ch_known_fusions    // channel: known fusions for SV filtering
    ch_repeatmasker_annotations // channel: repeatmasker annotations for SV filtering
    ch_gridss_config   // path: optional GRIDSS config file

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )

    ch_versions = ch_versions.mix(FASTQC.out.versions)

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

    ch_versions = ch_versions.mix(FASTP.out.versions)
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

        ch_versions = ch_versions.mix(UMI_PROCESSING.out.versions)

    } else {

        READ_ALIGNMENT(
            ch_input_reads,
            ch_genome_fasta,
            ch_genome_fai,
            ch_bwamem2_index
        )

        ch_multiqc_files = ch_multiqc_files.mix(READ_ALIGNMENT.out.dedup_metrics.collect{it[1]}.ifEmpty([]))
        ch_versions = ch_versions.mix(READ_ALIGNMENT.out.versions)
        ch_aligned_bam = READ_ALIGNMENT.out.dedup_bam
            .join(READ_ALIGNMENT.out.dedup_bai)
    }

    //
    // MODULE: QC of aligned BAM files
    //

    ALIGNMENT_QC(
        ch_aligned_bam,
        ch_genome_fasta,
        ch_genome_fai,
        ch_dict,
        ch_interval_list_slopped20
    )

    ch_versions = ch_versions.mix(ALIGNMENT_QC.out.versions)

    //
    // MODULE: Somatic SNV and INDELs Calling
    //

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

    // genome version for sage and gripss
    def genome_version = params.genome.equals("GRCh37") ? '37' : '38'

    SOMATIC_SNV_CALLING (
        ch_input_paired,
        ch_genome_fasta,
        ch_genome_fai,
        ch_dict,
        ch_germline_resource.collect{it -> it[1]},  // germline_resource
        ch_germline_resource_tbi.collect{it -> it[1]},  // germline_resource_tbi
        [],
        [],
        ch_interval_list_slopped20.collect{ it -> it[1] },
        ch_sage_known_hotspots_somatic,
        ch_sage_highconf_regions,
        ch_sage_pon,
        ch_ensembl_data_resources,
        ch_ensembl_vep_cache,
        genome_version
    )

    ch_versions = ch_versions.mix(SOMATIC_SNV_CALLING.out.versions)

    //
    // MODULE: CNV Calling
    //

    CNV_CALLING(
        ch_aligned_bam,
        ch_jumble_ref,
        ch_curation_ann
    )

    //
    // MODULE: SV Calling
    //

    SVS_CALLING(
        ch_aligned_bam,
        ch_genome_fasta,
        ch_genome_fai,
        ch_genome_gridss_index,
        ch_dict,
        ch_pon_breakends,
        ch_pon_breakpoints,
        ch_known_fusions,
        ch_repeatmasker_annotations,
        ch_targets_bed,
        ch_gridss_config,
        genome_version
    )

    //
    // Collate and save software versions
    //
    def topic_versions = Channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
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
    ch_multiqc_files = ch_multiqc_files.mix(ALIGNMENT_QC.out.multiple_metrics.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ALIGNMENT_QC.out.hs_metrics.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ALIGNMENT_QC.out.flagstat.collect{it[1]}.ifEmpty([]))


    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
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

    // Prepare final output channel

    autoseq_output = Channel
        .empty()
        .mix(
            ch_aligned_bam.map { meta, bam, _bai -> [ meta + [file: "bam"], bam] },
            ch_aligned_bam.map { meta, _bam, bai -> [ meta + [file: "bai"], bai] },
            ALIGNMENT_QC.out.flagstat.map { meta, flagstat -> [ meta + [file: "flagstat"], flagstat] },
            ALIGNMENT_QC.out.hs_metrics.map { meta, hs_metrics -> [ meta + [file: "hs_metrics"], hs_metrics] },
            ALIGNMENT_QC.out.multiple_metrics.map { meta, multiple_metrics -> [ meta + [file: "multiple_metrics"], multiple_metrics] },
            CNV_CALLING.out.jumble_cns.map { meta, cns -> [ meta + [file: "jumble_cns"], cns] },
            CNV_CALLING.out.cnr.map { meta, cnr -> [ meta + [file: "cnr"], cnr] },
            CNV_CALLING.out.seg.map { meta, seg -> [ meta + [file: "seg"], seg] },
            CNV_CALLING.out.profile_bedgraph.map { meta, profile_bedgraph -> [ meta + [file: "profile_bedgraph"], profile_bedgraph] },
            CNV_CALLING.out.segments_bedgraph.map { meta, segments_bedgraph -> [ meta + [file: "segments_bedgraph"], segments_bedgraph] },
            CNV_CALLING.out.png.map { meta, png -> [ meta + [file: "cnv_plot_png"], png] },
            CNV_CALLING.out.cns.map { meta, annotated_cns -> [ meta + [file: "annotated_cns"], annotated_cns] },
            SOMATIC_SNV_CALLING.out.contamination_table.map { meta, table -> [ meta + [file: "contamination_table"], table] },
            SOMATIC_SNV_CALLING.out.mutect2_stats.map { meta, stats -> [ meta + [file: "mutect2_stats"], stats] },
            SOMATIC_SNV_CALLING.out.mutect2_tbi.map { meta, tbi -> [ meta + [file: "mutect2_tbi"], tbi] },
            SOMATIC_SNV_CALLING.out.mutect2_vcf.map { meta, vcf -> [ meta + [file: "mutect2_vcf"], vcf] },
            SOMATIC_SNV_CALLING.out.sage_vcf.map { meta, vcf -> [ meta + [file: "sage_vcf"], vcf] },
            SOMATIC_SNV_CALLING.out.sage_tbi.map { meta, tbi -> [ meta + [file: "sage_tbi"], tbi] },
            SOMATIC_SNV_CALLING.out.somatic_vcf.map { meta, somatic_vcf -> [ meta + [file: "somatic_vcf"], somatic_vcf] },
            SOMATIC_SNV_CALLING.out.somatic_tbi.map { meta, somatic_tbi -> [ meta + [file: "somatic_tbi"], somatic_tbi] },
            SOMATIC_SNV_CALLING.out.vep_vcf.map { meta, vep_vcf -> [ meta + [file: "vep_vcf"], vep_vcf] },
            SOMATIC_SNV_CALLING.out.vep_tbi.map { meta, vep_tbi -> [ meta + [file: "vep_tbi"], vep_tbi] },
            SVS_CALLING.out.gripss_somatic_filtered_vcf.map { meta, vcf, tbi -> [ meta + [file: "gripss_somatic_filtered_vcf"], [vcf, tbi]] },
            SVS_CALLING.out.gripss_somatic_unfiltered_vcf.map { meta, vcf, tbi -> [ meta + [file: "gripss_somatic_unfiltered_vcf"], [vcf, tbi]] },
            SVS_CALLING.out.gripss_germline_filtered_vcf.map { meta, vcf, tbi -> [ meta + [file: "gripss_germline_filtered_vcf"], [vcf, tbi]] },
            SVS_CALLING.out.gripss_germline_unfiltered_vcf.map { meta, vcf, tbi -> [ meta + [file: "gripss_germline_unfiltered_vcf"], [vcf, tbi]] }
        )

    emit:
    autoseq_output = autoseq_output               // channel: [ val(meta + [file: description]), path(file) ]
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]


}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

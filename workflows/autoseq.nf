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

include { ALIGNMENT                                           } from '../subworkflows/local/alignment/main.nf'
include { FASTQ_CREATE_UMI_CONSENSUS_FGBIO as UMI_PROCESSING  } from '../subworkflows/nf-core/fastq_create_umi_consensus_fgbio/main'
include { QC_ALIGNMENT                                        } from '../subworkflows/local/qc_alignment/main.nf'
include { CALL_SOMATIC_SNVS                                   } from '../subworkflows/local/call_somatic_snvs/main.nf'
include { CALL_GERMLINE_SNVS                                  } from '../subworkflows/local/call_germline_snvs/main.nf'
include { CALL_CNVS                                           } from '../subworkflows/local/call_cnvs/main.nf'
include { CALL_SVS                                            } from '../subworkflows/local/call_svs/main.nf'
include { PROFILE_TUMOR_BIOMARKERS                            } from '../subworkflows/local/profile_tumor_biomarkers/main.nf'

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
    ch_dbsnp_vcf      // channel: optional dbSNP VCF for SNV annotation
    ch_dbsnp_vcf_tbi  // channel: optional dbSNP VCF

    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

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
    // SUBWORKFLOW: ALIGNMENT
    //

    if (params.umi_structure) {

        ch_input_reads
            .map { meta, reads ->
                def id = "${meta.case_id}.${meta.sample_name}".toString()
                meta   = meta + [id: id]
                return [meta.sample_name, [meta , reads]]
            }
            .groupTuple()
            .map { _sample_name, grouped_reads ->
                def metas = grouped_reads.collect{it -> it[0]}
                def files = grouped_reads.collect{it -> it[1]}.flatten()
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

        ALIGNMENT(
            ch_input_reads,
            ch_genome_fasta,
            ch_genome_fai,
            ch_bwamem2_index
        )

        ch_multiqc_files = ch_multiqc_files.mix(ALIGNMENT.out.dedup_metrics.collect{it -> it[1]}.ifEmpty([]))
        ch_versions = ch_versions.mix(ALIGNMENT.out.versions)
        ch_aligned_bam = ALIGNMENT.out.dedup_bam
            .join(ALIGNMENT.out.dedup_bai)
    }

    //
    // SUBWORKFLOW: QC of aligned BAM files
    //

    QC_ALIGNMENT(
        ch_aligned_bam,
        ch_genome_fasta,
        ch_genome_fai,
        ch_dict,
        ch_interval_list_slopped20
    )

    ch_versions = ch_versions.mix(QC_ALIGNMENT.out.versions)

    //
    // SUBWORKFLOW: Somatic SNV and INDELs Calling
    //

    // Branch samples by tumor/normal
    ch_aligned_bam
        .branch { meta, _bam, _bai ->
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
        .map { case_id, tumor_meta, tumor_bam, tumor_bai, normal_meta, normal_bam, normal_bai, _meta_intervals, intervals_file ->
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

    CALL_SOMATIC_SNVS (
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

    ch_versions = ch_versions.mix(CALL_SOMATIC_SNVS.out.versions)

    //
    // SUBWORKFLOW: GERMLINE Variant Calling
    //

    ch_germline_input = normal_ch
        .combine(ch_interval_list_slopped20)
        .map { _case_id, meta, bam, bai, _meta_intervals, intervals ->
            [meta, bam, bai, intervals, []]
        }


    CALL_GERMLINE_SNVS(
        ch_germline_input,  // normal samples only
        ch_genome_fasta,
        ch_genome_fai,
        ch_dict,
        ch_ensembl_vep_cache,
        ch_dbsnp_vcf,
        ch_dbsnp_vcf_tbi
    )

    //
    // SUBWORKFLOW: CNV Calling
    //

    CALL_CNVS(
        ch_aligned_bam,
        ch_jumble_ref,
        ch_curation_ann
    )

    //
    // SUBWORKFLOW: SV Calling
    //

    CALL_SVS(
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
    // SUBWORKFLOW: Tumor Biomarker Profiling (e.g. purity/ploidy, MSI, TMB, etc.)
    //

    ch_tumor_cnr = CALL_CNVS.out.cnr
        .filter { meta, _cnr -> meta.sample_type == "tumor" }

    ch_tumor_seg = CALL_CNVS.out.seg
        .filter { meta, _seg -> meta.sample_type == "tumor" }

    PROFILE_TUMOR_BIOMARKERS(
        ch_tumor_cnr,
        ch_tumor_seg,
        CALL_SOMATIC_SNVS.out.mutect2_unfiltered_vcf,
        ch_aligned_bam
    )

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
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
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it -> it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it -> it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_ALIGNMENT.out.multiple_metrics.collect{it -> it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_ALIGNMENT.out.hs_metrics.collect{it -> it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_ALIGNMENT.out.flagstat.collect{it -> it[1]}.ifEmpty([]))


    ch_multiqc_config        = channel.fromPath(
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

    autoseq_output = channel
        .empty()
        .mix(
            ch_aligned_bam.map { meta, bam, _bai -> [ meta + [file: "bam"], bam] },
            ch_aligned_bam.map { meta, _bam, bai -> [ meta + [file: "bai"], bai] },
            QC_ALIGNMENT.out.flagstat.map { meta, flagstat -> [ meta + [file: "flagstat"], flagstat] },
            QC_ALIGNMENT.out.hs_metrics.map { meta, hs_metrics -> [ meta + [file: "hs_metrics"], hs_metrics] },
            QC_ALIGNMENT.out.multiple_metrics.map { meta, multiple_metrics -> [ meta + [file: "multiple_metrics"], multiple_metrics] },
            CALL_CNVS.out.jumble_cns.map { meta, cns -> [ meta + [file: "jumble_cns"], cns] },
            CALL_CNVS.out.cnr.map { meta, cnr -> [ meta + [file: "cnr"], cnr] },
            CALL_CNVS.out.seg.map { meta, seg -> [ meta + [file: "seg"], seg] },
            CALL_CNVS.out.profile_bedgraph.map { meta, profile_bedgraph -> [ meta + [file: "profile_bedgraph"], profile_bedgraph] },
            CALL_CNVS.out.segments_bedgraph.map { meta, segments_bedgraph -> [ meta + [file: "segments_bedgraph"], segments_bedgraph] },
            CALL_CNVS.out.png.map { meta, png -> [ meta + [file: "cnv_plot_png"], png] },
            CALL_CNVS.out.cns.map { meta, annotated_cns -> [ meta + [file: "annotated_cns"], annotated_cns] },
            CALL_SOMATIC_SNVS.out.contamination_table.map { meta, table -> [ meta + [file: "contamination_table"], table] },
            CALL_SOMATIC_SNVS.out.mutect2_stats.map { meta, stats -> [ meta + [file: "mutect2_stats"], stats] },
            CALL_SOMATIC_SNVS.out.mutect2_tbi.map { meta, tbi -> [ meta + [file: "mutect2_tbi"], tbi] },
            CALL_SOMATIC_SNVS.out.mutect2_vcf.map { meta, vcf -> [ meta + [file: "mutect2_vcf"], vcf] },
            CALL_SOMATIC_SNVS.out.sage_vcf.map { meta, vcf -> [ meta + [file: "sage_vcf"], vcf] },
            CALL_SOMATIC_SNVS.out.sage_tbi.map { meta, tbi -> [ meta + [file: "sage_tbi"], tbi] },
            CALL_SOMATIC_SNVS.out.somatic_vcf.map { meta, somatic_vcf -> [ meta + [file: "somatic_vcf"], somatic_vcf] },
            CALL_SOMATIC_SNVS.out.somatic_tbi.map { meta, somatic_tbi -> [ meta + [file: "somatic_tbi"], somatic_tbi] },
            CALL_SOMATIC_SNVS.out.vep_vcf.map { meta, vep_vcf -> [ meta + [file: "vep_vcf"], vep_vcf] },
            CALL_SOMATIC_SNVS.out.vep_tbi.map { meta, vep_tbi -> [ meta + [file: "vep_tbi"], vep_tbi] },
            CALL_GERMLINE_SNVS.out.vcf.map { meta, vcf -> [ meta + [file: "germline_vcf"], vcf] },
            CALL_GERMLINE_SNVS.out.tbi.map { meta, tbi -> [ meta + [file: "germline_tbi"], tbi] },
            CALL_GERMLINE_SNVS.out.vep_vcf.map { meta, vep_vcf -> [ meta + [file: "germline_vep_vcf"], vep_vcf] },
            CALL_GERMLINE_SNVS.out.vep_tbi.map { meta, vep_tbi -> [ meta + [file: "germline_vep_tbi"], vep_tbi] },
            CALL_SVS.out.gripss_somatic_filtered_vcf.map { meta, vcf, tbi -> [ meta + [file: "gripss_somatic_filtered_vcf"], [vcf, tbi]] },
            CALL_SVS.out.gripss_somatic_unfiltered_vcf.map { meta, vcf, tbi -> [ meta + [file: "gripss_somatic_unfiltered_vcf"], [vcf, tbi]] },
            CALL_SVS.out.gripss_germline_filtered_vcf.map { meta, vcf, tbi -> [ meta + [file: "gripss_germline_filtered_vcf"], [vcf, tbi]] },
            CALL_SVS.out.gripss_germline_unfiltered_vcf.map { meta, vcf, tbi -> [ meta + [file: "gripss_germline_unfiltered_vcf"], [vcf, tbi]] },
            PROFILE_TUMOR_BIOMARKERS.out.dpyd_csv.map  { meta, csv  -> [ meta + [file: "dpyd_csv"],  csv  ] },
            PROFILE_TUMOR_BIOMARKERS.out.dpyd_json.map { meta, json -> [ meta + [file: "dpyd_json"], json ] }
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

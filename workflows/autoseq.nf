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
    val_multiqc_config                  // val: /path/to/multiqc_config.yaml
    val_multiqc_logo                    // val: /path/to/multiqc_logo.png
    val_multiqc_methods_description     // val: /path/to/multiqc_methods_description.md
    val_outdir                          // val: /path/to/output/directory

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )

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
                // groupTuple is unordered. Sort before concatenating to
                // ensure consistent FASTQ read order and reproducible
                // downstream results (alignments, UMIs, variants).
                def sorted = grouped_reads.sort(false) { it -> it[1].name }
                def metas = sorted.collect{it -> it[0]}
                def files = sorted.collect{it -> it[1]}.flatten()
                return [metas[0], files]
            }
            .set { ch_input_reads }

        CAT_FASTQ (
            ch_input_reads
        )

        // `.collect()` restores the value-channel semantics lost by `combine`, so that the
        // reference can be reused by every sample inside UMI_PROCESSING
        def ch_fasta_fai_dict = ch_genome_fasta
            .combine(ch_genome_fai)
            .combine(ch_dict)
            .map { meta, fasta, _meta_fai, fai, _meta_dict, dict ->
                [meta, fasta, fai, dict]
            }
            .collect()

        UMI_PROCESSING(
            CAT_FASTQ.out.reads,
            ch_fasta_fai_dict,
            ch_bwamem2_index,
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

    } else {

        ALIGNMENT(
            ch_input_reads,
            ch_genome_fasta,
            ch_genome_fai,
            ch_bwamem2_index
        )

        ch_multiqc_files = ch_multiqc_files.mix(ALIGNMENT.out.dedup_metrics.collect{it -> it[1]}.ifEmpty([]))
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

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${val_outdir}/pipeline_info",
            name:  'autoseq_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it -> it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it -> it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_ALIGNMENT.out.multiple_metrics.collect{it -> it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_ALIGNMENT.out.hs_metrics.collect{it -> it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(QC_ALIGNMENT.out.flagstat.collect{it -> it[1]}.ifEmpty([]))


    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = val_multiqc_methods_description
        ? file(val_multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))

    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'autoseq'],
                files,
                val_multiqc_config
                    ? file(val_multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                val_multiqc_logo ? file(val_multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )


    emit:
    annotated_cns                  = CALL_CNVS.out.cns                                              // channel: [ val(meta), path(cns) ]
    bai                            = ch_aligned_bam.map { meta, _bam, bai -> [ meta, bai ] }        // channel: [ val(meta), path(bai) ]
    bam                            = ch_aligned_bam.map { meta, bam, _bai -> [ meta, bam ] }        // channel: [ val(meta), path(bam) ]
    cnr                            = CALL_CNVS.out.cnr                                              // channel: [ val(meta), path(cnr) ]
    cnv_plot_png                   = CALL_CNVS.out.png                                              // channel: [ val(meta), path(png) ]
    contamination_table            = CALL_SOMATIC_SNVS.out.contamination_table                      // channel: [ val(meta), path(table) ]
    dpyd_csv                       = PROFILE_TUMOR_BIOMARKERS.out.dpyd_csv                          // channel: [ val(meta), path(csv) ]
    dpyd_json                      = PROFILE_TUMOR_BIOMARKERS.out.dpyd_json                         // channel: [ val(meta), path(json) ]
    flagstat                       = QC_ALIGNMENT.out.flagstat                                      // channel: [ val(meta), path(flagstat) ]
    germline_tbi                   = CALL_GERMLINE_SNVS.out.tbi                                     // channel: [ val(meta), path(tbi) ]
    germline_vcf                   = CALL_GERMLINE_SNVS.out.vcf                                     // channel: [ val(meta), path(vcf) ]
    germline_vep_tbi               = CALL_GERMLINE_SNVS.out.vep_tbi                                 // channel: [ val(meta), path(tbi) ]
    germline_vep_vcf               = CALL_GERMLINE_SNVS.out.vep_vcf                                 // channel: [ val(meta), path(vcf) ]
    gripss_germline_filtered_vcf   = CALL_SVS.out.gripss_germline_filtered_vcf                      // channel: [ val(meta), path(vcf), path(tbi) ]
    gripss_germline_unfiltered_vcf = CALL_SVS.out.gripss_germline_unfiltered_vcf                    // channel: [ val(meta), path(vcf), path(tbi) ]
    gripss_somatic_filtered_vcf    = CALL_SVS.out.gripss_somatic_filtered_vcf                       // channel: [ val(meta), path(vcf), path(tbi) ]
    gripss_somatic_unfiltered_vcf  = CALL_SVS.out.gripss_somatic_unfiltered_vcf                     // channel: [ val(meta), path(vcf), path(tbi) ]
    hs_metrics                     = QC_ALIGNMENT.out.hs_metrics                                    // channel: [ val(meta), path(metrics) ]
    jumble_cns                     = CALL_CNVS.out.jumble_cns                                       // channel: [ val(meta), path(cns) ]
    multiple_metrics               = QC_ALIGNMENT.out.multiple_metrics                              // channel: [ val(meta), path(metrics) ]
    multiqc_report                 = MULTIQC.out.report.map { _meta, report -> [report] }.toList()  // channel: /path/to/multiqc_report.html
    mutect2_stats                  = CALL_SOMATIC_SNVS.out.mutect2_stats                            // channel: [ val(meta), path(stats) ]
    mutect2_tbi                    = CALL_SOMATIC_SNVS.out.mutect2_tbi                              // channel: [ val(meta), path(tbi) ]
    mutect2_vcf                    = CALL_SOMATIC_SNVS.out.mutect2_vcf                              // channel: [ val(meta), path(vcf) ]
    profile_bedgraph               = CALL_CNVS.out.profile_bedgraph                                 // channel: [ val(meta), path(bedgraph) ]
    purecn_csv                     = PROFILE_TUMOR_BIOMARKERS.out.purecn_csv                        // channel: [ val(meta), path(csv) ]
    purecn_pdf                     = PROFILE_TUMOR_BIOMARKERS.out.purecn_pdf                        // channel: [ val(meta), path(pdf) ]
    sage_tbi                       = CALL_SOMATIC_SNVS.out.sage_tbi                                 // channel: [ val(meta), path(tbi) ]
    sage_vcf                       = CALL_SOMATIC_SNVS.out.sage_vcf                                 // channel: [ val(meta), path(vcf) ]
    seg                            = CALL_CNVS.out.seg                                              // channel: [ val(meta), path(seg) ]
    segments_bedgraph              = CALL_CNVS.out.segments_bedgraph                                // channel: [ val(meta), path(bedgraph) ]
    somatic_tbi                    = CALL_SOMATIC_SNVS.out.somatic_tbi                              // channel: [ val(meta), path(tbi) ]
    somatic_vcf                    = CALL_SOMATIC_SNVS.out.somatic_vcf                              // channel: [ val(meta), path(vcf) ]
    vep_tbi                        = CALL_SOMATIC_SNVS.out.vep_tbi                                  // channel: [ val(meta), path(tbi) ]
    vep_vcf                        = CALL_SOMATIC_SNVS.out.vep_vcf                                  // channel: [ val(meta), path(vcf) ]
    versions                       = ch_versions                                                    // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    genomic-medicine-sweden/autoseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/genomic-medicine-sweden/autoseq
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { AUTOSEQ                 } from './workflows/autoseq'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_autoseq_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_autoseq_pipeline'
include { PREPARE_REFERENCES      } from './subworkflows/local/prepare_references/main'
include { channelFromPathWithMeta } from './subworkflows/local/utils_nfcore_autoseq_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow GENOMICMEDICINESWEDEN_AUTOSEQ {

    take:
    samplesheet                         // channel: [mandatory] samplesheet read in from --input
    val_bwamem2_index                   // string:  [mandatory] path to reference genome bwa-mem2 index
    val_curation_ann                    // string:  [mandatory] path to CNV curation annotation CSV
    val_dbsnp_vcf                       // string:  [optional]  path to dbSNP VCF
    val_dbsnp_vcf_tbi                   // string:  [optional]  path to dbSNP VCF index
    val_ensembl_vep_cache               // string:  [mandatory] path to Ensembl VEP cache directory
    val_ensembl_vep_cache_tar           // string:  [optional]  path to Ensembl VEP cache tarball
    val_genome_dict                     // string:  [mandatory] path to reference genome sequence dictionary
    val_genome_fai                      // string:  [mandatory] path to reference genome FASTA index
    val_genome_fasta                    // string:  [mandatory] path to reference genome FASTA
    val_germline_resource               // string:  [mandatory] path to germline resource VCF
    val_germline_resource_tbi           // string:  [mandatory] path to germline resource VCF index
    val_gridss_config                   // string:  [mandatory] path to GRIDSS configuration file
    val_gridss_index                    // string:  [mandatory] path to GRIDSS genome index directory
    val_gridss_index_tar                // string:  [optional]  path to GRIDSS genome index tarball
    val_gridss_known_fusions            // string:  [mandatory] path to GRIDSS known fusions
    val_gridss_pon_breakends            // string:  [mandatory] path to GRIDSS PON breakends
    val_gridss_pon_breakpoints          // string:  [mandatory] path to GRIDSS PON breakpoints
    val_gridss_repeatmasker_annotations // string:  [mandatory] path to GRIDSS RepeatMasker annotations
    val_hmf_ensembl_data                // string:  [mandatory] path to HMF Ensembl data resources directory
    val_hmf_ensembl_data_tar            // string:  [optional]  path to HMF Ensembl data resources tarball
    val_interval_list                   // string:  [mandatory] path to target regions interval list
    val_jumble_ref                      // string:  [mandatory] path to Jumble reference RDS
    val_multiqc_config                  // string:  [optional]  path to MultiQC config
    val_multiqc_logo                    // string:  [optional]  path to MultiQC logo
    val_multiqc_methods_description     // string:  [optional]  path to MultiQC methods description
    val_outdir                          // string:  [mandatory] path to output directory
    val_sage_highconf_regions           // string:  [mandatory] path to SAGE high-confidence regions
    val_sage_known_hotspots_somatic     // string:  [mandatory] path to SAGE somatic hotspots
    val_sage_pon                        // string:  [mandatory] path to SAGE panel of normals
    val_targets_bed                     // string:  [mandatory] path to target regions BED

    main:

    // Minimal reference preparation workflow
    def ch_references = PREPARE_REFERENCES (
                            val_genome_fasta,
                            val_bwamem2_index,
                            val_ensembl_vep_cache,
                            val_ensembl_vep_cache_tar,
                            val_gridss_index,
                            val_gridss_index_tar,
                            val_hmf_ensembl_data,
                            val_hmf_ensembl_data_tar
                        )


    //
    // Initialise channels for reference genome
    //
    ch_genome_fasta            = ch_references.genome_fasta
    ch_bwamem2_index           = ch_references.bwamem2_index
    ch_ensembl_vep_cache       = ch_references.vep_cache
    ch_hmf_ensembl_data        = ch_references.hmf_ensembl_data
    ch_gridss_index            = ch_references.gridss_index

    // channelFromPathWithMeta() builds a [[id:simpleName], file] channel from a path,
    // or channel.empty() when the path is null

    ch_genome_fai                      = channelFromPathWithMeta(val_genome_fai)
    ch_dict                            = channelFromPathWithMeta(val_genome_dict)
    ch_dbsnp_vcf                       = channelFromPathWithMeta(val_dbsnp_vcf)
    ch_dbsnp_vcf_tbi                   = channelFromPathWithMeta(val_dbsnp_vcf_tbi)
    ch_germline_resource               = channelFromPathWithMeta(val_germline_resource)
    ch_germline_resource_tbi           = channelFromPathWithMeta(val_germline_resource_tbi)

    ch_targets_bed                     = channelFromPathWithMeta(val_targets_bed)
    ch_interval_list                   = channelFromPathWithMeta(val_interval_list)
    ch_jumble_ref                      = channelFromPathWithMeta(val_jumble_ref)

    ch_sage_known_hotspots_somatic     = channelFromPathWithMeta(val_sage_known_hotspots_somatic)
    ch_sage_highconf_regions           = channelFromPathWithMeta(val_sage_highconf_regions)
    ch_sage_pon                        = channelFromPathWithMeta(val_sage_pon)
    ch_curation_ann                    = channelFromPathWithMeta(val_curation_ann)

    // GRIDSS-specific channels for SV calling
    ch_gridss_pon_breakends            = channelFromPathWithMeta(val_gridss_pon_breakends)
    ch_gridss_pon_breakpoints          = channelFromPathWithMeta(val_gridss_pon_breakpoints)
    ch_gridss_known_fusions            = channelFromPathWithMeta(val_gridss_known_fusions)
    ch_gridss_repeatmasker_annotations = channelFromPathWithMeta(val_gridss_repeatmasker_annotations)
    ch_gridss_config                   = channelFromPathWithMeta(val_gridss_config)


    //
    // WORKFLOW: Run pipeline
    //
    AUTOSEQ (
        samplesheet,
        ch_genome_fasta,
        ch_genome_fai,
        ch_dict,
        ch_bwamem2_index,
        ch_targets_bed,
        ch_interval_list,
        ch_jumble_ref,
        ch_sage_known_hotspots_somatic,
        ch_sage_highconf_regions,
        ch_sage_pon,
        ch_ensembl_vep_cache,
        ch_hmf_ensembl_data,
        ch_curation_ann,
        ch_germline_resource,
        ch_germline_resource_tbi,
        ch_gridss_index,
        ch_gridss_pon_breakends,
        ch_gridss_pon_breakpoints,
        ch_gridss_known_fusions,
        ch_gridss_repeatmasker_annotations,
        ch_gridss_config,
        ch_dbsnp_vcf,
        ch_dbsnp_vcf_tbi,
        val_multiqc_config,
        val_multiqc_logo,
        val_multiqc_methods_description,
        val_outdir,
    )

    emit:
    annotated_cns                  = AUTOSEQ.out.annotated_cns                  // channel: [ val(meta), path(cns) ]
    bai                            = AUTOSEQ.out.bai                            // channel: [ val(meta), path(bai) ]
    bam                            = AUTOSEQ.out.bam                            // channel: [ val(meta), path(bam) ]
    cnr                            = AUTOSEQ.out.cnr                            // channel: [ val(meta), path(cnr) ]
    cnv_plot_png                   = AUTOSEQ.out.cnv_plot_png                   // channel: [ val(meta), path(png) ]
    contamination_table            = AUTOSEQ.out.contamination_table            // channel: [ val(meta), path(table) ]
    dpyd_csv                       = AUTOSEQ.out.dpyd_csv                       // channel: [ val(meta), path(csv) ]
    dpyd_json                      = AUTOSEQ.out.dpyd_json                      // channel: [ val(meta), path(json) ]
    flagstat                       = AUTOSEQ.out.flagstat                       // channel: [ val(meta), path(flagstat) ]
    germline_taf_tbi               = AUTOSEQ.out.germline_taf_tbi               // channel: [ val(meta), path(tbi) ]
    germline_taf_vcf               = AUTOSEQ.out.germline_taf_vcf               // channel: [ val(meta), path(vcf) ]
    germline_tbi                   = AUTOSEQ.out.germline_tbi                   // channel: [ val(meta), path(tbi) ]
    germline_vcf                   = AUTOSEQ.out.germline_vcf                   // channel: [ val(meta), path(vcf) ]
    germline_vep_tbi               = AUTOSEQ.out.germline_vep_tbi               // channel: [ val(meta), path(tbi) ]
    germline_vep_vcf               = AUTOSEQ.out.germline_vep_vcf               // channel: [ val(meta), path(vcf) ]
    gripss_germline_filtered_vcf   = AUTOSEQ.out.gripss_germline_filtered_vcf   // channel: [ val(meta), path(vcf), path(tbi) ]
    gripss_germline_unfiltered_vcf = AUTOSEQ.out.gripss_germline_unfiltered_vcf // channel: [ val(meta), path(vcf), path(tbi) ]
    gripss_somatic_filtered_vcf    = AUTOSEQ.out.gripss_somatic_filtered_vcf    // channel: [ val(meta), path(vcf), path(tbi) ]
    gripss_somatic_unfiltered_vcf  = AUTOSEQ.out.gripss_somatic_unfiltered_vcf  // channel: [ val(meta), path(vcf), path(tbi) ]
    hetsnps_tbi                    = AUTOSEQ.out.hetsnps_tbi                    // channel: [ val(meta), path(tbi) ]
    hetsnps_vcf                    = AUTOSEQ.out.hetsnps_vcf                    // channel: [ val(meta), path(vcf) ]
    hs_metrics                     = AUTOSEQ.out.hs_metrics                     // channel: [ val(meta), path(metrics) ]
    jumble_cns                     = AUTOSEQ.out.jumble_cns                     // channel: [ val(meta), path(cns) ]
    multiple_metrics               = AUTOSEQ.out.multiple_metrics               // channel: [ val(meta), path(metrics) ]
    multiqc_report                 = AUTOSEQ.out.multiqc_report                 // channel: /path/to/multiqc_report.html
    mutect2_stats                  = AUTOSEQ.out.mutect2_stats                  // channel: [ val(meta), path(stats) ]
    mutect2_tbi                    = AUTOSEQ.out.mutect2_tbi                    // channel: [ val(meta), path(tbi) ]
    mutect2_vcf                    = AUTOSEQ.out.mutect2_vcf                    // channel: [ val(meta), path(vcf) ]
    profile_bedgraph               = AUTOSEQ.out.profile_bedgraph               // channel: [ val(meta), path(bedgraph) ]
    purecn_csv                     = AUTOSEQ.out.purecn_csv                     // channel: [ val(meta), path(csv) ]
    purecn_pdf                     = AUTOSEQ.out.purecn_pdf                     // channel: [ val(meta), path(pdf) ]
    sage_tbi                       = AUTOSEQ.out.sage_tbi                       // channel: [ val(meta), path(tbi) ]
    sage_vcf                       = AUTOSEQ.out.sage_vcf                       // channel: [ val(meta), path(vcf) ]
    seg                            = AUTOSEQ.out.seg                            // channel: [ val(meta), path(seg) ]
    segments_bedgraph              = AUTOSEQ.out.segments_bedgraph              // channel: [ val(meta), path(bedgraph) ]
    somatic_tbi                    = AUTOSEQ.out.somatic_tbi                    // channel: [ val(meta), path(tbi) ]
    somatic_vcf                    = AUTOSEQ.out.somatic_vcf                    // channel: [ val(meta), path(vcf) ]
    vep_tbi                        = AUTOSEQ.out.vep_tbi                        // channel: [ val(meta), path(tbi) ]
    vep_vcf                        = AUTOSEQ.out.vep_vcf                        // channel: [ val(meta), path(vcf) ]
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )


    //
    // WORKFLOW: Run main workflow
    //
    GENOMICMEDICINESWEDEN_AUTOSEQ (
        PIPELINE_INITIALISATION.out.samplesheet,
        params.bwamem2_index,
        params.curation_ann,
        params.dbsnp_vcf,
        params.dbsnp_vcf_tbi,
        params.ensembl_vep_cache,
        params.ensembl_vep_cache_tar,
        params.genome_dict,
        params.genome_fai,
        params.genome_fasta,
        params.germline_resource,
        params.germline_resource_tbi,
        params.gridss_config,
        params.gridss_index,
        params.gridss_index_tar,
        params.gridss_known_fusions,
        params.gridss_pon_breakends,
        params.gridss_pon_breakpoints,
        params.gridss_repeatmasker_annotations,
        params.hmf_ensembl_data,
        params.hmf_ensembl_data_tar,
        params.interval_list,
        params.jumble_ref,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir,
        params.sage_highconf_regions,
        params.sage_known_hotspots_somatic,
        params.sage_pon,
        params.targets_bed
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.multiqc_report
    )

    //
    // Mix together the outputs that share an output directory — one item per output directory
    //
    def ch_alignment = channel.empty().mix(
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.bai,
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.bam
    )

    def ch_cnv = channel.empty().mix(
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.annotated_cns,
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.cnr,
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.cnv_plot_png,
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.jumble_cns,
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.profile_bedgraph,
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.seg,
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.segments_bedgraph
    )

    def ch_dpyd = channel.empty().mix(
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.dpyd_csv,
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.dpyd_json
    )

    def ch_purecn = channel.empty().mix(
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.purecn_csv,
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.purecn_pdf
    )

    def ch_qc_picard = channel.empty().mix(
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.hs_metrics,
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.multiple_metrics
    )

    def ch_snvs_germline = channel.empty().mix(
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.germline_taf_tbi,
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.germline_taf_vcf,
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.germline_tbi,
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.germline_vcf,
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.germline_vep_tbi,
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.germline_vep_vcf,
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.hetsnps_tbi,
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.hetsnps_vcf
    )

    def ch_snvs_somatic = channel.empty().mix(
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.mutect2_stats,
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.mutect2_tbi,
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.mutect2_vcf,
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.sage_tbi,
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.sage_vcf,
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.somatic_tbi,
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.somatic_vcf,
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.vep_tbi,
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.vep_vcf
    )

    def ch_svs_germline = channel.empty().mix(
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.gripss_germline_filtered_vcf,
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.gripss_germline_unfiltered_vcf
    )

    def ch_svs_somatic = channel.empty().mix(
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.gripss_somatic_filtered_vcf,
        GENOMICMEDICINESWEDEN_AUTOSEQ.out.gripss_somatic_unfiltered_vcf
    )

    publish:
    alignment        = ch_alignment                                          // channel: [ val(meta), path(file) ]
    cnv              = ch_cnv                                                // channel: [ val(meta), path(file) ]
    dpyd             = ch_dpyd                                               // channel: [ val(meta), path(file) ]
    multiqc_report   = GENOMICMEDICINESWEDEN_AUTOSEQ.out.multiqc_report      // channel: /path/to/multiqc_report.html
    purecn           = ch_purecn                                             // channel: [ val(meta), path(file) ]
    qc_contamination = GENOMICMEDICINESWEDEN_AUTOSEQ.out.contamination_table // channel: [ val(meta), path(table) ]
    qc_picard        = ch_qc_picard                                          // channel: [ val(meta), path(metrics) ]
    qc_samtools      = GENOMICMEDICINESWEDEN_AUTOSEQ.out.flagstat            // channel: [ val(meta), path(flagstat) ]
    snvs_germline    = ch_snvs_germline                                      // channel: [ val(meta), path(file) ]
    snvs_somatic     = ch_snvs_somatic                                       // channel: [ val(meta), path(file) ]
    svs_germline     = ch_svs_germline                                       // channel: [ val(meta), path(vcf), path(tbi) ]
    svs_somatic      = ch_svs_somatic                                        // channel: [ val(meta), path(vcf), path(tbi) ]
}


output {
    alignment {
        path { "alignment" }
    }
    cnv {
        path { "cnv" }
    }
    dpyd {
        path { "dpyd" }
    }
    multiqc_report {
        path { "multiqc" }
    }
    purecn {
        path { "purecn" }
    }
    qc_contamination {
        path { "qc/contamination" }
    }
    qc_picard {
        path { "qc/picard" }
    }
    qc_samtools {
        path { "qc/samtools" }
    }
    snvs_germline {
        path { "snvs/germline" }
    }
    snvs_somatic {
        path { "snvs/somatic" }
    }
    svs_germline {
        path { "svs/germline" }
    }
    svs_somatic {
        path { "svs/somatic" }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

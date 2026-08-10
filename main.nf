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

    // Using channelFromPathWithMeta helper (with simpleName as meta id).
    // If filepath is null, returns, channel.empty())

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
    ch_hmf_ensembl_data                = channelFromPathWithMeta(val_hmf_ensembl_data)
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
    autoseq_output = AUTOSEQ.out.autoseq_output      // channel: [ val(meta + [file: description]), path(file) ]
    multiqc_report = AUTOSEQ.out.multiqc_report     // channel: /path/to/multiqc_report.html
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

    publish:
    autoseq_output  = GENOMICMEDICINESWEDEN_AUTOSEQ.out.autoseq_output  // channel: [ val(meta + [file: description]), path(file) ]
    multiqc_report  = GENOMICMEDICINESWEDEN_AUTOSEQ.out.multiqc_report   // channel: /path/to/multiqc_report.html
}


output {
    multiqc_report {
        path { "multiqc" }
    }
    autoseq_output {
        path { meta, _file ->
            if (meta.file == 'bam' || meta.file == 'bai') {
                return 'alignment'
            } else if (meta.file == 'flagstat') {
                return 'qc/samtools'
            } else if (meta.file == 'contamination_table') {
                return 'qc/contamination'
            } else if (meta.file == 'hs_metrics' || meta.file == 'multiple_metrics') {
                return 'qc/picard'
            } else if (meta.file == 'jumble_cns' || meta.file == 'cnr' || meta.file == 'seg' ||
                    meta.file == 'profile_bedgraph' || meta.file == 'segments_bedgraph' ||
                    meta.file == 'annotated_cns' || meta.file == 'cnv_plot_png') {
                return 'cnv'
            } else if (meta.file == 'mutect2_stats' || meta.file == 'mutect2_vcf' || meta.file == 'mutect2_tbi') {
                return 'variants/somatic/mutect2'
            } else if (meta.file == 'sage_vcf' || meta.file == 'sage_tbi') {
                return 'variants/somatic/sage'
            } else if (meta.file == 'somatic_vcf' || meta.file == 'somatic_tbi') {
                return 'variants/somatic/merged'
            } else if (meta.file == 'vep_vcf' || meta.file == 'vep_tbi') {
                return 'variants/somatic'
            } else if (meta.file == 'germline_vcf' || meta.file == 'germline_tbi') {
                return 'variants/germline/haplotypecaller'
            } else if (meta.file == 'germline_vep_vcf' || meta.file == 'germline_vep_tbi') {
                return 'variants/germline/'
            } else if (meta.file == 'gripss_somatic_filtered_vcf' || meta.file == 'gripss_somatic_unfiltered_vcf') {
                return 'svs/somatic/'
            } else if (meta.file == 'gripss_germline_filtered_vcf' || meta.file == 'gripss_germline_unfiltered_vcf') {
                return 'svs/germline/'
            } else if (meta.file == "dpyd_csv" || meta.file == "dpyd_json") {
                return 'dpyd'
            } else if (meta.file == "purecn_csv" || meta.file == "purecn_pdf") {
                return 'purecn'
            } else {
                return ''
            }
        }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

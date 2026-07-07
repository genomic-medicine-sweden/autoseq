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
include { channelFromPathWithMeta } from './subworkflows/local/utils_nfcore_autoseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NXF_AUTOSEQ {

    take:
    samplesheet                         // channel: samplesheet read in from --input
    val_genome_fasta                    // val: /path/to/genome.fasta
    val_genome_fai                      // val: /path/to/genome.fasta.fai
    val_genome_dict                     // val: /path/to/genome.dict
    val_bwamem2_index                   // val: /path/to/bwamem2_index
    val_dbsnp_vcf                       // val: /path/to/dbsnp.vcf
    val_dbsnp_vcf_tbi                   // val: /path/to/dbsnp.vcf.tbi
    val_targets_bed                     // val: /path/to/targets.bed
    val_interval_list                   // val: /path/to/interval_list
    val_jumble_ref                      // val: /path/to/jumble_ref
    val_ensembl_vep_cache               // val: /path/to/vep
    val_hmf_ensembl_data                // val: /path/to/ensembl_data
    val_curation_ann                    // val: /path/to/curation_ann
    val_germline_resource               // val: /path/to/germline_resource.vcf
    val_germline_resource_tbi           // val: /path/to/germline_resource.vcf.tbi
    val_sage_known_hotspots_somatic     // val: /path/to/sage_known_hotspots_somatic
    val_sage_highconf_regions           // val: /path/to/sage_highconf_regions
    val_sage_pon                        // val: /path/to/sage_pon
    val_gridss_index                    // val: /path/to/gridss_index
    val_gridss_pon_breakends            // val: /path/to/pon_breakends
    val_gridss_pon_breakpoints          // val: /path/to/pon_breakpoints
    val_gridss_known_fusions            // val: /path/to/known_fusions
    val_gridss_repeatmasker_annotations // val: /path/to/repeatmasker_annotations
    val_gridss_config                   // val: /path/to/gridss_config
    val_multiqc_config                  // val: /path/to/multiqc_config.yaml
    val_multiqc_logo                    // val: /path/to/multiqc_logo.png
    val_multiqc_methods_description     // val: /path/to/multiqc_methods_description.md
    val_outdir                          // val: /path/to/output/directory

    main:

    //
    // Initialise channels for reference genome
    //
    // Using channelFromPathWithMeta helper (with simpleName as meta id).
    // If filepath is null, returns, channel.empty())
    ch_genome_fasta                = channelFromPathWithMeta(val_genome_fasta)
    ch_genome_fai                  = channelFromPathWithMeta(val_genome_fai)
    ch_dict                        = channelFromPathWithMeta(val_genome_dict)
    ch_bwamem2_index               = channelFromPathWithMeta(val_bwamem2_index)
    ch_dbsnp_vcf                   = channelFromPathWithMeta(val_dbsnp_vcf)
    ch_dbsnp_vcf_tbi               = channelFromPathWithMeta(val_dbsnp_vcf_tbi)
    ch_germline_resource           = channelFromPathWithMeta(val_germline_resource)
    ch_germline_resource_tbi       = channelFromPathWithMeta(val_germline_resource_tbi)

    ch_targets_bed                 = channelFromPathWithMeta(val_targets_bed)
    ch_interval_list               = channelFromPathWithMeta(val_interval_list)
    ch_jumble_ref                  = channelFromPathWithMeta(val_jumble_ref)
    ch_ensembl_vep_cache           = channelFromPathWithMeta(val_ensembl_vep_cache)

    ch_sage_known_hotspots_somatic = channelFromPathWithMeta(val_sage_known_hotspots_somatic)
    ch_sage_highconf_regions       = channelFromPathWithMeta(val_sage_highconf_regions)
    ch_sage_pon                    = channelFromPathWithMeta(val_sage_pon)
    ch_hmf_ensembl_data            = channelFromPathWithMeta(val_hmf_ensembl_data)
    ch_curation_ann                = channelFromPathWithMeta(val_curation_ann)

    // GRIDSS-specific channels for SV calling
    ch_gridss_index                    = channelFromPathWithMeta(val_gridss_index)
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
    NXF_AUTOSEQ (
        PIPELINE_INITIALISATION.out.samplesheet,
        params.genome_fasta,
        params.genome_fai,
        params.genome_dict,
        params.bwamem2_index,
        params.targets_bed,
        params.interval_list,
        params.jumble_ref,
        params.sage_known_hotspots_somatic,
        params.sage_highconf_regions,
        params.sage_pon,
        params.ensembl_vep_cache,
        params.hmf_ensembl_data,
        params.curation_ann,
        params.germline_resource,
        params.germline_resource_tbi,
        params.gridss_index,
        params.gridss_pon_breakends,
        params.gridss_pon_breakpoints,
        params.gridss_known_fusions,
        params.gridss_repeatmasker_annotations,
        params.gridss_config,
        params.dbsnp_vcf,
        params.dbsnp_vcf_tbi,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir
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
        NXF_AUTOSEQ.out.multiqc_report
    )

    publish:
    autoseq_output  = NXF_AUTOSEQ.out.autoseq_output  // channel: [ val(meta + [file: description]), path(file) ]
    multiqc_report  = NXF_AUTOSEQ.out.multiqc_report   // channel: /path/to/multiqc_report.html
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

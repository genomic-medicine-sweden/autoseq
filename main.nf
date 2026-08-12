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
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_autoseq_pipeline'
include { getPanelsAttribute      } from './subworkflows/local/utils_nfcore_autoseq_pipeline'
include { channelFromPathWithMeta } from './subworkflows/local/utils_nfcore_autoseq_pipeline'
include { PREPARE_REFERENCES      } from './subworkflows/local/prepare_references/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`
params.ref_genome_fasta                 = getGenomeAttribute('fasta')
params.ref_genome_fai                   = getGenomeAttribute('fai')
params.ref_genome_dict                  = getGenomeAttribute('dict')
params.bwamem2_index                    = getGenomeAttribute('bwamem2_index')
params.dbsnp_vcf                        = getGenomeAttribute('dbsnp_vcf')
params.dbsnp_vcf_tbi                    = getGenomeAttribute('dbsnp_vcf_tbi')
params.germline_resource                = getGenomeAttribute('germline_resource')
params.germline_resource_tbi            = getGenomeAttribute('germline_resource_tbi')
params.sage_known_hotspots_somatic      = getGenomeAttribute('sage_known_hotspots_somatic')
params.sage_highconf_regions            = getGenomeAttribute('sage_highconf_regions')
params.sage_pon                         = getGenomeAttribute('sage_pon')
params.ensembl_vep_cache                = getGenomeAttribute('ensembl_vep_cache')
params.ensembl_data_resources           = getGenomeAttribute('ensembl_data_resources')
params.curation_ann                     = getGenomeAttribute('curation_annotations')
params.genome_gridss_index              = getGenomeAttribute('gridss_index')
params.gridss_config                    = getGenomeAttribute('gridss_config')
params.gridss_pon_breakends             = getGenomeAttribute('gridss_pon_breakends')
params.gridss_pon_breakpoints           = getGenomeAttribute('gridss_pon_breakpoints')
params.gridss_known_fusions             = getGenomeAttribute('gridss_known_fusions')
params.gridss_repeatmasker_annotations  = getGenomeAttribute('gridss_repeatmasker_annotations')


params.targets_bed             = getPanelsAttribute('targets_bed_slopped20')
params.interval_list_slopped20 = getPanelsAttribute('targets_interval_list_slopped20')
params.jumble_ref              = getPanelsAttribute('jumble_ref')


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
    val_multiqc_config                  // string:  [optional]  path to MultiQC config
    val_multiqc_logo                    // string:  [optional]  path to MultiQC logo
    val_multiqc_methods_description     // string:  [optional]  path to MultiQC methods description
    val_outdir                          // string:  [mandatory] path to output directory

    main:

    // Minimal reference preparation workflow
    def ch_references = PREPARE_REFERENCES (
                            params.ref_genome_fasta,
                            params.bwamem2_index,
                            params.ensembl_vep_cache,
                            params.ensembl_vep_cache_tar,
                            params.genome_gridss_index,
                            params.gridss_index_tar,
                            params.ensembl_data_resources,
                            params.hmf_ensembl_data_tar
                        )


    //
    // Initialise channels for reference genome
    //
    ch_genome_fasta            = ch_references.genome_fasta
    ch_bwamem2_index           = ch_references.bwamem2_index
    ch_ensembl_vep_cache       = ch_references.vep_cache
    ch_ensembl_data_resources  = ch_references.hmf_ensembl_data
    ch_genome_gridss_index     = ch_references.gridss_index

    // channelFromPathWithMeta() builds a [[id:simpleName], file] channel from a path,
    // or channel.empty() when the path is null
    ch_genome_fai    = channelFromPathWithMeta(params.ref_genome_fai)
    ch_dict          = channelFromPathWithMeta(params.ref_genome_dict)
    ch_dbsnp_vcf     = channelFromPathWithMeta(params.dbsnp_vcf)
    ch_dbsnp_vcf_tbi = channelFromPathWithMeta(params.dbsnp_vcf_tbi)

    //
    ch_targets_bed             = channelFromPathWithMeta(params.targets_bed)
    ch_interval_list_slopped20 = channelFromPathWithMeta(params.interval_list_slopped20)
    ch_jumble_ref              = channelFromPathWithMeta(params.jumble_ref)

    //
    ch_sage_known_hotspots_somatic = channelFromPathWithMeta(params.sage_known_hotspots_somatic)
    ch_sage_highconf_regions       = channelFromPathWithMeta(params.sage_highconf_regions)
    ch_sage_pon                    = channelFromPathWithMeta(params.sage_pon)
    ch_curation_ann                = channelFromPathWithMeta(params.curation_ann)
    ch_germline_resource           = channelFromPathWithMeta(params.germline_resource)
    ch_germline_resource_tbi       = channelFromPathWithMeta(params.germline_resource_tbi)

    // GRIDSS-specific channels for SV calling
    ch_pon_breakends            = channelFromPathWithMeta(params.gridss_pon_breakends)
    ch_pon_breakpoints          = channelFromPathWithMeta(params.gridss_pon_breakpoints)
    ch_known_fusions            = channelFromPathWithMeta(params.gridss_known_fusions)
    ch_repeatmasker_annotations = channelFromPathWithMeta(params.gridss_repeatmasker_annotations)
    ch_gridss_config            = channelFromPathWithMeta(params.gridss_config)

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
        ch_interval_list_slopped20,
        ch_jumble_ref,
        ch_sage_known_hotspots_somatic,
        ch_sage_highconf_regions,
        ch_sage_pon,
        ch_ensembl_vep_cache,
        ch_ensembl_data_resources,
        ch_curation_ann,
        ch_germline_resource,
        ch_germline_resource_tbi,
        ch_genome_gridss_index,
        ch_pon_breakends,
        ch_pon_breakpoints,
        ch_known_fusions,
        ch_repeatmasker_annotations,
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

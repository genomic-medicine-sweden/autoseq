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
workflow NXF_AUTOSEQ {

    take:
    samplesheet                         // channel: samplesheet read in from --input
    val_multiqc_config                  // val: /path/to/multiqc_config.yaml
    val_multiqc_logo                    // val: /path/to/multiqc_logo.png
    val_multiqc_methods_description     // val: /path/to/multiqc_methods_description.md
    val_outdir                          // val: /path/to/output/directory

    main:

    //
    // Initialise channels for reference genome
    //
    ch_genome_fasta  = params.ref_genome_fasta  ? channel.fromPath(params.ref_genome_fasta).map{ it -> [[id:'genome_fasta'], it]}.collect() : channel.empty()
    ch_genome_fai    = params.ref_genome_fai    ? channel.fromPath(params.ref_genome_fai).map{ it -> [[id:'genome_fai'], it]}.collect() : channel.empty()
    ch_dict          = params.ref_genome_dict   ? channel.fromPath(params.ref_genome_dict).map{ it -> [[id:'genome_dict'], it]}.collect() : channel.empty()
    ch_bwamem2_index = params.bwamem2_index     ? channel.fromPath(params.bwamem2_index).map{ it -> [[id:'bwamem2_index'], it]}.collect() : channel.empty()
    ch_dbsnp_vcf     = params.dbsnp_vcf        ? channel.fromPath(params.dbsnp_vcf).map{ it -> [[id:'dbsnp_vcf'], it]}.collect() : channel.empty()
    ch_dbsnp_vcf_tbi = params.dbsnp_vcf_tbi  ? channel.fromPath(params.dbsnp_vcf_tbi).map{ it -> [[id:'dbsnp_vcf_tbi'], it]}.collect() : channel.empty()

    //
    ch_targets_bed             = params.targets_bed ? channel.fromPath(params.targets_bed).map{ it -> [[id:'targets_bed'], it]}.collect() : channel.empty()
    ch_interval_list_slopped20 = params.interval_list_slopped20 ? channel.fromPath(params.interval_list_slopped20).map{ it -> [[id:'interval_list_slopped20'], it]}.collect() : channel.empty()
    ch_jumble_ref              = params.jumble_ref ? channel.fromPath(params.jumble_ref).map{ it -> [[id:'jumble_ref'], it]}.collect() : channel.empty()
    ch_ensembl_vep_cache       = params.ensembl_vep_cache ? channel.fromPath(params.ensembl_vep_cache).map{ it -> [[id:'ensembl_vep_cache'], it]}.collect() : channel.empty()

    //
    ch_sage_known_hotspots_somatic = params.sage_known_hotspots_somatic ? channel.fromPath(params.sage_known_hotspots_somatic).map{ it -> [[id:'sage_known_hotspots_somatic'], it]}.collect() : channel.empty()
    ch_sage_highconf_regions       = params.sage_highconf_regions ? channel.fromPath(params.sage_highconf_regions).map{ it -> [[id:'sage_highconf_regions'], it]}.collect() : channel.empty()
    ch_sage_pon                    = params.sage_pon ? channel.fromPath(params.sage_pon).map{ it -> [[id:'sage_pon'], it]}.collect() : channel.empty()
    ch_ensembl_data_resources      = params.ensembl_data_resources ? channel.fromPath(params.ensembl_data_resources).map{ it -> [[id:'ensembl_data_resources'], it]}.collect() : channel.empty()
    ch_curation_ann                = params.curation_ann ? channel.fromPath(params.curation_ann).map{ it -> [[id:'curation_ann'], it]}.collect() : channel.empty()
    ch_germline_resource           = params.germline_resource ? channel.fromPath(params.germline_resource).map{ it -> [[id:'germline_resource'], it]}.collect() : channel.empty()
    ch_germline_resource_tbi       = params.germline_resource_tbi ? channel.fromPath(params.germline_resource_tbi).map{ it -> [[id:'germline_resource_tbi'], it]}.collect() : channel.empty()

    // GRIDSS-specific channels for SV calling
    ch_genome_gridss_index      = params.genome_gridss_index ? channel.fromPath(params.genome_gridss_index).map{ it -> [[id:'genome_gridss_index'], it]}.collect() : channel.empty()
    ch_pon_breakends            = params.gridss_pon_breakends ? channel.fromPath(params.gridss_pon_breakends).map{ it -> [[id:'pon_breakends'], it]}.collect() : channel.empty()
    ch_pon_breakpoints          = params.gridss_pon_breakpoints ? channel.fromPath(params.gridss_pon_breakpoints).map{ it -> [[id:'pon_breakpoints'], it]}.collect() : channel.empty()
    ch_known_fusions            = params.gridss_known_fusions ? channel.fromPath(params.gridss_known_fusions).map{ it -> [[id:'known_fusions'], it]}.collect() : channel.empty()
    ch_repeatmasker_annotations = params.gridss_repeatmasker_annotations ? channel.fromPath(params.gridss_repeatmasker_annotations).map{ it -> [[id:'repeatmasker_annotations'], it]}.collect() : channel.empty()
    ch_gridss_config            = params.gridss_config ? channel.fromPath(params.gridss_config).map{ it -> [[id: 'gridss_config'], it]}.collect() : channel.empty()

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
    NXF_AUTOSEQ (
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

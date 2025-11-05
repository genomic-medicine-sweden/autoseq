#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/autoseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/autoseq
    Website: https://nf-co.re/autoseq
    Slack  : https://nfcore.slack.com/channels/autoseq
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

// TODO nf-core: Remove this line if you don't need a FASTA file
//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`
params.ref_genome_fasta            = getGenomeAttribute('fasta')
params.ref_genome_fai              = getGenomeAttribute('fai')
params.ref_genome_dict             = getGenomeAttribute('dict')
params.bwamem2_index               = getGenomeAttribute('bwamem2_index')
params.sage_known_hotspots_somatic = getGenomeAttribute('sage_known_hotspots_somatic')
params.sage_highconf_regions       = getGenomeAttribute('sage_highconf_regions')
params.sage_pon                    = getGenomeAttribute('sage_pon')
params.ensembl_data_resources      = getGenomeAttribute('ensembl_data_resources')


params.targets_bed             = getPanelsAttribute('targets_bed_slopped20')
params.targets_bed_gz          = getPanelsAttribute('targets_bed_slopped20_gz')
params.interval_list           = getPanelsAttribute('targets_interval_list')
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
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // Initialise channels for reference genome
    //
    ch_genome_fasta  = params.ref_genome_fasta  ? Channel.fromPath(params.ref_genome_fasta).map{ it -> [[id:'genome_fasta'], it]}.collect() : Channel.empty()
    ch_genome_fai    = params.ref_genome_fai    ? Channel.fromPath(params.ref_genome_fai).map{ it -> [[id:'genome_fai'], it]}.collect() : Channel.empty()
    ch_dict          = params.ref_genome_dict   ? Channel.fromPath(params.ref_genome_dict).map{ it -> [[id:'genome_dict'], it]}.collect() : Channel.empty()
    ch_bwamem2_index = params.bwamem2_index     ? Channel.fromPath(params.bwamem2_index).map{ it -> [[id:'bwamem2_index'], it]}.collect() : Channel.empty()

    //
    ch_targets_bed             = params.targets_bed ? Channel.fromPath(params.targets_bed).map{ it -> [[id:'targets_bed'], it]}.collect() : Channel.empty()
    ch_targets_bed_gz          = params.targets_bed_gz ? Channel.fromPath(params.targets_bed_gz).map{ it -> [[id:'targets_bed_gz'], it]}.collect() : Channel.empty()
    ch_interval_list           = params.interval_list ? Channel.fromPath(params.interval_list).map{ it -> [[id:'interval_list'], it]}.collect() : Channel.empty()
    ch_interval_list_slopped20 = params.interval_list_slopped20 ? Channel.fromPath(params.interval_list_slopped20).map{ it -> [[id:'interval_list_slopped20'], it]}.collect() : Channel.empty()
    ch_jumble_ref              = params.jumble_ref ? Channel.fromPath(params.jumble_ref).map{ it -> [[id:'jumble_ref'], it]}.collect() : Channel.empty()

    //
    ch_sage_known_hotspots_somatic = params.sage_known_hotspots_somatic ? Channel.fromPath(params.sage_known_hotspots_somatic).map{ it -> [[id:'sage_known_hotspots_somatic'], it]}.collect() : Channel.empty()
    ch_sage_highconf_regions       = params.sage_highconf_regions ? Channel.fromPath(params.sage_highconf_regions).map{ it -> [[id:'sage_highconf_regions'], it]}.collect() : Channel.empty()
    ch_sage_pon                    = params.sage_pon ? Channel.fromPath(params.sage_pon).map{ it -> [[id:'sage_pon'], it]}.collect() : Channel.empty()
    ch_ensembl_data_resources      = params.ensembl_data_resources ? Channel.fromPath(params.ensembl_data_resources).map{ it -> [[id:'ensembl_data_resources'], it]}.collect() : Channel.empty()

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
        ch_targets_bed_gz,
        ch_interval_list,
        ch_interval_list_slopped20,
        ch_jumble_ref,
        ch_sage_known_hotspots_somatic,
        ch_sage_highconf_regions,
        ch_sage_pon,
        ch_ensembl_data_resources
    )
    emit:
    multiqc_report = AUTOSEQ.out.multiqc_report // channel: /path/to/multiqc_report.html
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
        PIPELINE_INITIALISATION.out.samplesheet
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
        params.hook_url,
        NXF_AUTOSEQ.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { SAGE_SOMATIC } from '../../../modules/local/sage/somatic/main'
include { BAM_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING_GATK as CALL_GATK_MUTECT2 } from '../../../subworkflows/nf-core/bam_tumor_normal_somatic_variant_calling_gatk/main'



workflow CALL_SOMATIC_SNVS {
    take:
    ch_input                 // channel: [ val(meta), path(input), path(input_index), val(which_norm) ]
    ch_fasta                 // channel: [ val(meta), path(fasta) ]
    ch_fai                   // channel: [ val(meta), path(fai) ]
    ch_dict                  // channel: [ val(meta), path(dict) ]
    ch_germline_resource     // channel: /path/to/germline/resource
    ch_germline_resource_tbi // channel: /path/to/germline/index
    ch_panel_of_normals      // channel: /path/to/panel/of/normals
    ch_panel_of_normals_tbi  // channel: /path/to/panel/of/normals/index
    ch_interval_file         // channel: /path/to/interval/file
    ch_sage_known_hotspots_somatic
    ch_sage_highconf_regions
    ch_sage_pon
    ch_ensembl_data_resources

    main:
    versions = Channel.empty()

    // Call somatic SNVs using GATK Mutect2 in tumor-normal mode
    CALL_GATK_MUTECT2(
        ch_input,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_germline_resource,
        ch_germline_resource_tbi,
        ch_panel_of_normals,
        ch_panel_of_normals_tbi,
        ch_interval_file
    )

    versions = versions.mix(CALL_GATK_MUTECT2.out.versions)

    genome_version  = '37' // hardcoded for now, could be passed in as a param if needed
    // Call somatic SNVs using SAGE in tumor-normal mode
    SAGE_SOMATIC(
        ch_input,
        ch_fasta,
        ch_fai,
        ch_dict,
        genome_version, // genome version
        ch_panel_of_normals,
        ch_sage_known_hotspots_somatic.collect{it[1]},
        ch_sage_highconf_regions.collect{it[1]},
        ch_ensembl_data_resources.collect{it[1]},
        true // targeted_mode
    )

    versions = versions.mix(SAGE_SOMATIC.out.versions)


    emit:
    contamination_table = CALL_GATK_MUTECT2.out.contamination_table    // channel: [ val(meta), path(table) ]
    mutect2_stats       = CALL_GATK_MUTECT2.out.filtered_stats                 // channel: [ val(meta), path(stats) ]
    mutect2_tbi         = CALL_GATK_MUTECT2.out.filtered_tbi                   // channel: [ val(meta), path(tbi) ]
    mutect2_vcf         = CALL_GATK_MUTECT2.out.filtered_vcf                   // channel: [ val(meta), path(vcf) ]
    pileup_table_normal = CALL_GATK_MUTECT2.out.pileup_table_normal         // channel: [ val(meta), path(table) ]
    pileup_table_tumor  = CALL_GATK_MUTECT2.out.pileup_table_tumor          // channel: [ val(meta), path(table) ]
    sage_vcf            = SAGE_SOMATIC.out.vcf                                 // channel: [ val(meta), path(vcf) ]
    sage_tbi            = SAGE_SOMATIC.out.tbi                                 // channel: [
    versions
}

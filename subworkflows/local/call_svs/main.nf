//
// Subworkflow for structural variant calling using GRIDSS and filtering with GRIPSS
//

include { GRIDSS_EXTRACT_OVERLAPPING_FRAGMENTS  } from '../../../modules/local/gridss/extract_overlapping_fragments/main'
include { GRIDSS_PREPROCESS                     } from '../../../modules/local/gridss/preprocess/main'
include { GRIDSS_ASSEMBLE                       } from '../../../modules/local/gridss/assemble/main'
include { GRIDSS_CALL                           } from '../../../modules/local/gridss/call/main'
include { GRIPSS_SOMATIC                        } from '../../../modules/local/gripss/somatic/main'
include { GRIPSS_GERMLINE                       } from '../../../modules/local/gripss/germline/main'

workflow SVS_CALLING {
    take:
    ch_aligned_bam
    ch_target_bed
    ch_genome_fasta
    ch_genome_fasta_fai
    ch_genome_gridss_index
    ch_genome_dict
    ch_blacklist
    ch_pon_breakends
    ch_pon_breakpoints
    ch_known_fusions
    ch_repeatmasker_annotations
    ch_target_region_bed
    gridss_config
    val_genome_version

    main:
    ch_versions = Channel.empty()

    //
    // GRIDSS: Extract overlapping fragments from tumor BAM
    //
    GRIDSS_EXTRACT_OVERLAPPING_FRAGMENTS (
        ch_aligned_bam,
        ch_target_bed
    )
    ch_versions = ch_versions.mix(GRIDSS_EXTRACT_OVERLAPPING_FRAGMENTS.out.versions)

    //
    // GRIDSS: Preprocess step
    //
    GRIDSS_PREPROCESS (
        GRIDSS_EXTRACT_OVERLAPPING_FRAGMENTS.out.gridss_targeted_bam,
        ch_genome_fasta,
        ch_genome_gridss_index,
        ch_genome_fasta_fai,
        ch_genome_dict,
        gridss_config
    )
    ch_versions = ch_versions.mix(GRIDSS_PREPROCESS.out.versions)


    ch_assemble_input = GRIDSS_EXTRACT_OVERLAPPING_FRAGMENTS.out.gridss_targeted_bam
        .join(GRIDSS_EXTRACT_OVERLAPPING_FRAGMENTS.out.gridss_targeted_bai)
        .join(GRIDSS_PREPROCESS.out.preprocess_dir)
        .view()

    //
    // GRIDSS: Assemble step
    //
    GRIDSS_ASSEMBLE (
        ch_assemble_input,
        ch_genome_fasta,
        ch_genome_gridss_index,
        ch_genome_fasta_fai,
        ch_genome_dict,
        ch_blacklist,
        gridss_config
    )
    ch_versions = ch_versions.mix(GRIDSS_ASSEMBLE.out.versions)


    ch_call_input = GRIDSS_EXTRACT_OVERLAPPING_FRAGMENTS.out.gridss_targeted_bam
        .join(GRIDSS_EXTRACT_OVERLAPPING_FRAGMENTS.out.gridss_targeted_bai)
        .join(GRIDSS_ASSEMBLE.out.assemble_dir)
        .view()
    //
    // GRIDSS: Call step to generate SV VCF
    //
    GRIDSS_CALL (
        ch_call_input,
        ch_genome_fasta,
        ch_genome_gridss_index,
        ch_genome_fasta_fai,
        ch_genome_dict,
        ch_blacklist,
        gridss_config
    )
    ch_versions = ch_versions.mix(GRIDSS_CALL.out.versions)

    //
    // GRIPSS: Somatic variant filtering
    //
    GRIPSS_SOMATIC (
        GRIDSS_CALL.out.vcf,
        ch_genome_fasta,
        ch_genome_fasta_fai,
        ch_pon_breakends,
        ch_pon_breakpoints,
        ch_known_fusions,
        ch_repeatmasker_annotations,
        ch_target_region_bed,
        val_genome_version
    )
    ch_versions = ch_versions.mix(GRIPSS_SOMATIC.out.versions)

    //
    // GRIPSS: Germline variant filtering
    //
    GRIPSS_GERMLINE (
        GRIDSS_CALL.out.vcf,
        ch_genome_fasta,
        ch_genome_fasta_fai,
        ch_pon_breakends,
        ch_pon_breakpoints,
        ch_known_fusions,
        ch_repeatmasker_annotations,
        ch_target_region_bed,
        val_genome_version
    )
    ch_versions = ch_versions.mix(GRIPSS_GERMLINE.out.versions)

    emit:
    gridss_vcf = GRIDSS_CALL.out.vcf_file
    gripss_somatic_filtered_vcf = GRIPSS_SOMATIC.out.filtered_vcf
    gripss_somatic_unfiltered_vcf = GRIPSS_SOMATIC.out.unfiltered_vcf
    gripss_germline_filtered_vcf = GRIPSS_GERMLINE.out.filtered_vcf
    gripss_germline_unfiltered_vcf = GRIPSS_GERMLINE.out.unfiltered_vcf
    versions = ch_versions
}

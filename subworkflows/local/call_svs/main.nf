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
    ch_genome_fasta
    ch_genome_fai
    ch_genome_gridss_index
    ch_genome_dict
    ch_pon_breakends
    ch_pon_breakpoints
    ch_known_fusions
    ch_repeatmasker_annotations
    ch_target_region_bed
    ch_gridss_config
    genome_version

    main:

    //
    // GRIDSS: Extract overlapping fragments from tumor BAM
    //
    GRIDSS_EXTRACT_OVERLAPPING_FRAGMENTS (
        ch_aligned_bam,
        ch_target_region_bed
    )

    //
    // GRIDSS: Preprocess step
    //
    GRIDSS_PREPROCESS (
        GRIDSS_EXTRACT_OVERLAPPING_FRAGMENTS.out.gridss_targeted_bam,
        ch_genome_fasta,
        ch_genome_gridss_index,
        ch_genome_fai,
        ch_genome_dict,
        ch_gridss_config  // Convert to value channel
    )

    ch_assemble_input = GRIDSS_EXTRACT_OVERLAPPING_FRAGMENTS.out.gridss_targeted_bam
        .join(GRIDSS_PREPROCESS.out.preprocess_dir)
        .map { meta, bam, bai, preprocess_dir ->
            [meta.case_id, [meta, bam, bai, preprocess_dir]]
        }
        .groupTuple()
        .map { case_id, samples ->
            // Separate normal and tumor samples
            def normal_sample = samples.find { sample -> sample[0].sample_type == 'normal' }
            def tumor_sample = samples.find { sample -> sample[0].sample_type == 'tumor' }

            // Extract components
            def meta = [
                id: case_id,
                normal_id: normal_sample[0].id,
                tumor_id: tumor_sample[0].id
            ]
            def normal_bam = normal_sample[1]
            def tumor_bam = tumor_sample[1]
            def normal_bai = normal_sample[2]
            def tumor_bai = tumor_sample[2]
            def normal_process_dir = normal_sample[3]
            def tumor_process_dir = tumor_sample[3]

            // Return in desired format
            [meta, [normal_bam, tumor_bam], [normal_bai, tumor_bai], [normal_process_dir, tumor_process_dir]]
        }

    //
    // GRIDSS: Assemble step
    //
    GRIDSS_ASSEMBLE (
        ch_assemble_input,
        ch_genome_fasta,
        ch_genome_gridss_index,
        ch_genome_fai,
        ch_genome_dict,
        ch_gridss_config
    )


    ch_call_input = ch_assemble_input
        .join(GRIDSS_ASSEMBLE.out.assemble_dir)

    //
    // GRIDSS: Call step to generate SV VCF
    //
    GRIDSS_CALL (
        ch_call_input,
        ch_genome_fasta,
        ch_genome_gridss_index,
        ch_genome_fai,
        ch_genome_dict,
        ch_gridss_config
    )

    //
    // GRIPSS: Somatic variant filtering
    //
    GRIPSS_SOMATIC (
        GRIDSS_CALL.out.vcf,
        ch_genome_fasta,
        ch_genome_fai,
        ch_pon_breakends,
        ch_pon_breakpoints,
        ch_known_fusions,
        ch_repeatmasker_annotations,
        ch_target_region_bed,
        genome_version
    )

    //
    // GRIPSS: Germline variant filtering
    //
    GRIPSS_GERMLINE (
        GRIDSS_CALL.out.vcf,
        ch_genome_fasta,
        ch_genome_fai,
        ch_pon_breakends,
        ch_pon_breakpoints,
        ch_known_fusions,
        ch_repeatmasker_annotations,
        ch_target_region_bed,
        genome_version
    )

    emit:
    gridss_vcf                      = GRIDSS_CALL.out.vcf
    gripss_somatic_filtered_vcf     = GRIPSS_SOMATIC.out.filtered_vcf
    gripss_somatic_unfiltered_vcf   = GRIPSS_SOMATIC.out.unfiltered_vcf
    gripss_germline_filtered_vcf    = GRIPSS_GERMLINE.out.filtered_vcf
    gripss_germline_unfiltered_vcf  = GRIPSS_GERMLINE.out.unfiltered_vcf
}

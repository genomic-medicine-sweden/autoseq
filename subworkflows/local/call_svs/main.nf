//
// Subworkflow for structural variant calling using GRIDSS and filtering with GRIPSS
//

include { GRIDSS_EXTRACTOVERLAPPINGFRAGMENTS } from '../../../modules/nf-core/gridss/extractoverlappingfragments/main'
include { GRIDSS_PREPROCESS                  } from '../../../modules/nf-core/gridss/preprocess/main'
include { GRIDSS_ASSEMBLE                    } from '../../../modules/nf-core/gridss/assemble/main'
include { GRIDSS_CALL                        } from '../../../modules/nf-core/gridss/call/main'
include { SAMTOOLS_INDEX                     } from '../../../modules/nf-core/samtools/index/main'
include { GRIPSS_SOMATIC                     } from '../../../modules/local/gripss/somatic/main'
include { GRIPSS_GERMLINE                    } from '../../../modules/local/gripss/germline/main'

workflow CALL_SVS {
    take:
    ch_aligned_bam
    ch_genome_fasta
    ch_genome_fai
    ch_genome_gridss_index
    ch_pon_breakends
    ch_pon_breakpoints
    ch_known_fusions
    ch_repeatmasker_annotations
    ch_target_region_bed
    ch_gridss_config
    genome_version

    main:

    // The nf-core GRIDSS modules take the reference as a single tuple. `.collect()`
    // restores the value-channel semantics lost by `combine`, so that the reference
    // can be reused by every sample.
    def ch_fasta_fai_gridss_index = ch_genome_fasta
        .combine(ch_genome_fai)
        .combine(ch_genome_gridss_index)
        .map { meta, fasta, _meta_fai, fai, _meta_gridss_index, gridss_index ->
            [meta, fasta, fai, gridss_index]
        }
        .collect()

    //
    // GRIDSS: Extract overlapping fragments from the aligned BAM
    //
    GRIDSS_EXTRACTOVERLAPPINGFRAGMENTS (
        ch_aligned_bam,
        ch_target_region_bed
    )

    //
    // SAMTOOLS: Index the targeted BAM, the GRIDSS steps require an indexed input
    //
    SAMTOOLS_INDEX (
        GRIDSS_EXTRACTOVERLAPPINGFRAGMENTS.out.bam
    )

    ch_targeted_bam = GRIDSS_EXTRACTOVERLAPPINGFRAGMENTS.out.bam
        .join(SAMTOOLS_INDEX.out.bai)

    //
    // GRIDSS: Preprocess step
    //
    GRIDSS_PREPROCESS (
        ch_targeted_bam,
        ch_fasta_fai_gridss_index
    )

    ch_assemble_input = ch_targeted_bam
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
        ch_fasta_fai_gridss_index,
        ch_gridss_config
    )

    ch_call_input = ch_assemble_input
        .join(GRIDSS_ASSEMBLE.out.assemble_dir)

    //
    // GRIDSS: Call step to generate SV VCF
    //
    GRIDSS_CALL (
        ch_call_input,
        ch_fasta_fai_gridss_index,
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

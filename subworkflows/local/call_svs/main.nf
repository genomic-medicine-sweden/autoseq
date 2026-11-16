//
// Subworkflow for structural variant calling using GRIDSS and filtering with GRIPSS
//

include { GRIDSS_ASSEMBLE                    } from '../../../modules/nf-core/gridss/assemble/main'
include { GRIDSS_CALL                        } from '../../../modules/nf-core/gridss/call/main'
include { GRIDSS_EXTRACTOVERLAPPINGFRAGMENTS } from '../../../modules/nf-core/gridss/extractoverlappingfragments/main'
include { GRIDSS_PREPROCESS                  } from '../../../modules/nf-core/gridss/preprocess/main'
include { GRIPSS_GERMLINE                    } from '../../../modules/local/gripss/germline/main'
include { GRIPSS_SOMATIC                     } from '../../../modules/local/gripss/somatic/main'
include { SAMTOOLS_INDEX                     } from '../../../modules/nf-core/samtools/index/main'

workflow CALL_SVS {
    take:
    ch_bam_bai                  // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_genome_fai               // channel: [mandatory] [ val(meta), path(fai) ]
    ch_genome_fasta             // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_gridss_config            // channel: [optional]  [ val(meta), path(properties) ]
    ch_gridss_index             // channel: [mandatory] [ val(meta), path(index) ]
    ch_known_fusions            // channel: [mandatory] [ val(meta), path(bedpe) ]
    ch_pon_breakends            // channel: [mandatory] [ val(meta), path(bed) ]
    ch_pon_breakpoints          // channel: [mandatory] [ val(meta), path(bedpe) ]
    ch_repeatmasker_annotations // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_target_region_bed        // channel: [mandatory] [ val(meta), path(bed) ]
    val_genome_version          // string:  [mandatory] numeric genome version for GRIPSS ("37" or "38")

    main:

    // The nf-core GRIDSS modules take the reference as a single tuple. `.collect()`
    // restores the value-channel semantics lost by `combine`, so that the reference
    // can be reused by every sample.
    def ch_fasta_fai_gridss_index = ch_genome_fasta
        .combine(ch_genome_fai)
        .combine(ch_gridss_index)
        .map { meta, fasta, _meta_fai, fai, _meta_gridss_index, gridss_index ->
            [meta, fasta, fai, gridss_index]
        }
        .collect()

    //
    // MODULE: Extract the fragments overlapping the target regions
    //

    GRIDSS_EXTRACTOVERLAPPINGFRAGMENTS(
        ch_bam_bai,
        ch_target_region_bed,
    )

    //
    // MODULE: Index the targeted BAM, the GRIDSS steps require an indexed input
    //

    SAMTOOLS_INDEX(
        GRIDSS_EXTRACTOVERLAPPINGFRAGMENTS.out.bam
    )

    def ch_samtools_index_for_gridss_preprocess = GRIDSS_EXTRACTOVERLAPPINGFRAGMENTS.out.bam
        .join(SAMTOOLS_INDEX.out.index)

    //
    // MODULE: GRIDSS preprocess
    //

    GRIDSS_PREPROCESS(
        ch_samtools_index_for_gridss_preprocess,
        ch_fasta_fai_gridss_index,
    )

    //
    // MODULE: GRIDSS assemble, tumor and normal samples are paired on `case_id`
    //

    def ch_gridss_preprocess_for_gridss_assemble = ch_samtools_index_for_gridss_preprocess
        .join(GRIDSS_PREPROCESS.out.preprocess_dir)
        .map { meta, bam, bai, preprocess_dir ->
            [meta.case_id, [meta, bam, bai, preprocess_dir]]
        }
        .groupTuple()
        .map { case_id, samples ->
            def normal_sample = samples.find { sample -> sample[0].sample_type == 'normal' }
            def tumor_sample  = samples.find { sample -> sample[0].sample_type == 'tumor' }

            def meta = [
                id: case_id,
                normal_id: normal_sample[0].id,
                tumor_id: tumor_sample[0].id,
            ]

            [
                meta,
                [normal_sample[1], tumor_sample[1]], // bam
                [normal_sample[2], tumor_sample[2]], // bai
                [normal_sample[3], tumor_sample[3]], // preprocess_dir
            ]
        }

    GRIDSS_ASSEMBLE(
        ch_gridss_preprocess_for_gridss_assemble,
        ch_fasta_fai_gridss_index,
        ch_gridss_config,
    )

    //
    // MODULE: GRIDSS call to generate the SV VCF
    //

    def ch_gridss_assemble_for_gridss_call = ch_gridss_preprocess_for_gridss_assemble
        .join(GRIDSS_ASSEMBLE.out.assemble_dir)

    GRIDSS_CALL(
        ch_gridss_assemble_for_gridss_call,
        ch_fasta_fai_gridss_index,
        ch_gridss_config,
    )

    //
    // MODULE: GRIPSS somatic SV filtering
    //

    GRIPSS_SOMATIC(
        GRIDSS_CALL.out.vcf,
        ch_genome_fasta,
        ch_genome_fai,
        ch_pon_breakends,
        ch_pon_breakpoints,
        ch_known_fusions,
        ch_repeatmasker_annotations,
        ch_target_region_bed,
        val_genome_version,
    )

    //
    // MODULE: GRIPSS germline SV filtering
    //

    GRIPSS_GERMLINE(
        GRIDSS_CALL.out.vcf,
        ch_genome_fasta,
        ch_genome_fai,
        ch_pon_breakends,
        ch_pon_breakpoints,
        ch_known_fusions,
        ch_repeatmasker_annotations,
        ch_target_region_bed,
        val_genome_version,
    )

    emit:
    gridss_vcf                     = GRIDSS_CALL.out.vcf                // channel: [ val(meta), path(vcf) ]
    gripss_germline_filtered_vcf   = GRIPSS_GERMLINE.out.filtered_vcf   // channel: [ val(meta), path(vcf), path(tbi) ]
    gripss_germline_unfiltered_vcf = GRIPSS_GERMLINE.out.unfiltered_vcf // channel: [ val(meta), path(vcf), path(tbi) ]
    gripss_somatic_filtered_vcf    = GRIPSS_SOMATIC.out.filtered_vcf    // channel: [ val(meta), path(vcf), path(tbi) ]
    gripss_somatic_unfiltered_vcf  = GRIPSS_SOMATIC.out.unfiltered_vcf  // channel: [ val(meta), path(vcf), path(tbi) ]
}

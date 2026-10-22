//
// Subworkflow for alignment QC metrics
//

include { PICARD_COLLECTHSMETRICS       } from '../../../modules/nf-core/picard/collecthsmetrics/main'
include { PICARD_COLLECTMULTIPLEMETRICS } from '../../../modules/nf-core/picard/collectmultiplemetrics/main'
include { SAMTOOLS_FLAGSTAT             } from '../../../modules/nf-core/samtools/flagstat/main'

workflow QC_ALIGNMENT {
    take:
    ch_bam_bai       // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_genome_dict   // channel: [mandatory] [ val(meta), path(dict) ]
    ch_genome_fai    // channel: [mandatory] [ val(meta), path(fai) ]
    ch_genome_fasta  // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_interval_list // channel: [mandatory] [ val(meta), path(interval_list) ]

    main:

    //
    // MODULE: Collect BAM metrics with Picard
    //

    PICARD_COLLECTMULTIPLEMETRICS(
        ch_bam_bai,
        ch_genome_fasta,
        ch_genome_fai
    )

    ch_bam_for_picard_collecthsmetrics = ch_bam_bai
        .combine(ch_interval_list)
        .map { meta, bam, bai, _meta2, interval_list ->
            [meta, bam, bai, interval_list, interval_list]
        }


    //
    // MODULE: Collect HS metrics with Picard
    //

    PICARD_COLLECTHSMETRICS(
        ch_bam_for_picard_collecthsmetrics,
        ch_genome_fasta,
        ch_genome_fai,
        ch_genome_dict,
        [[], []]
    )

    //
    // MODULE: Collect BAM metrics with Samtools
    //

    SAMTOOLS_FLAGSTAT(
        ch_bam_bai
    )

    emit:
    flagstat         = SAMTOOLS_FLAGSTAT.out.flagstat            // channel: [ val(meta), path(flagstat) ]
    hs_metrics       = PICARD_COLLECTHSMETRICS.out.metrics       // channel: [ val(meta), path(metrics) ]
    multiple_metrics = PICARD_COLLECTMULTIPLEMETRICS.out.metrics // channel: [ val(meta), [ path(metrics) ] ]

}

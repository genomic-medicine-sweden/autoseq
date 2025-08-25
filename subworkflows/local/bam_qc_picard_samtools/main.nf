

include { PICARD_COLLECTMULTIPLEMETRICS  } from '../../../modules/nf-core/picard/collectmultiplemetrics/main'
include { PICARD_COLLECTHSMETRICS        } from '../../../modules/nf-core/picard/collecthsmetrics/main'
include { SAMTOOLS_FLAGSTAT              } from '../../../modules/nf-core/samtools/flagstat/main'

workflow BAM_QC_PICARD_SAMTOOLS {
    take:
    ch_input_bam // channel: input BAM files from alignment workflow
    ch_genome_fasta
    ch_genome_fai
    
    main:

    ch_versions = Channel.empty()
    
    //
    // MODULE: Collect BAM metrics with Picard
    //

    PICARD_COLLECTMULTIPLEMETRICS(
        ch_input_bam,
        ch_genome_fasta,
        ch_bam_meta
    )

    ch_versions = ch_versions.mix(PICARD_COLLECTMULTIPLEMETRICS.out.versions.first())

    //
    // MODULE: Collect HS metrics with Picard
    //

    PICARD_COLLECTHSMETRICS(
        ch_input_bam,
        ch_genome_fasta,
        ch_bam_meta
    )

    ch_versions = ch_versions.mix(PICARD_COLLECTHSMETRICS.out.versions.first())
    ch_versions = ch_versions.mix(PICARD_COLLECTHSMETRICS.out.versions.second())
    
    //
    // MODULE: Collect BAM metrics with Samtools
    //

    SAMTOOLS_FLAGSTAT(
        ch_input_bam,
        ch_bam_meta
    )

    emit:
    multiple_metrics = PICARD_COLLECTMULTIPLEMETRICS.out.metrics
    hs_metrics       = PICARD_COLLECTHSMETRICS.out.metrics
    flagstat         = SAMTOOLS_FLAGSTAT.out.flagstat 
    versions = ch_versions.mix(PICARD_COLLECTHSMETRICS.out.versions.first())

} 
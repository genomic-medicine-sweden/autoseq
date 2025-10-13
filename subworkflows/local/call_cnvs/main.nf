

include { JUMBLE_RUN         } from '../../../modules/local/jumble/main'
include { ANNOTATE_CNVS      } from '../../../modules/local/annotate_cnvs/main'


workflow CALL_CNVS_JUMBLE {
    take:
    ch_input_bam            // channel: input BAM files from alignment workflow
    ch_jumble_ref           // channel: Jumble reference files
    ch_curation_ann         // channel: CNV curation/annotation file

    main:

    ch_versions = Channel.empty()

    //
    // MODULE: JUMBLE_RUN to call CNVs
    //

    JUMBLE_RUN(
        ch_input_bam,
        ch_jumble_ref
    )

    ch_versions = ch_versions.mix(JUMBLE_RUN.out.versions.first())

    //
    // MODULE: Annotate CNVs with cancer relevant genes
    //

    ANNOTATE_CNVS(
        JUMBLE_RUN.out.cns,
        ch_curation_ann
    )

    ch_versions = ch_versions.mix(JUMBLE_RUN.out.versions.first())

    emit:
    jumble_cns          = JUMBLE_RUN.out.cns
    cnr                 = JUMBLE_RUN.out.cnr
    seg                 = JUMBLE_RUN.out.seg
    profile_bedgraph    = JUMBLE_RUN.out.profile_bedgraph
    segments_bedgraph   = JUMBLE_RUN.out.segments_bedgraph
    png                 = JUMBLE_RUN.out.png
    cns                 = ANNOTATE_CNVS.out.cns
    versions            = ch_versions

}

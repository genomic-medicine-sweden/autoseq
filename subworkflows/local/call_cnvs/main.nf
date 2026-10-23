//
// Subworkflow for Copy Number Variants calling
//

include { ANNOTATE_CNVS } from '../../../modules/local/annotate_cnvs/main'
include { JUMBLE_RUN    } from '../../../modules/local/jumble/run/main'


workflow CALL_CNVS {
    take:
    ch_bam_bai      // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_curation_ann // channel: [mandatory] [ val(meta), path(csv) ]
    ch_jumble_ref   // channel: [mandatory] [ val(meta), path(rds) ]

    main:

    //
    // MODULE: JUMBLE_RUN to call CNVs
    //

    JUMBLE_RUN(
        ch_bam_bai,
        ch_jumble_ref,
    )

    //
    // MODULE: Annotate CNVs with cancer relevant genes
    //

    def ch_jumble_run_for_annotate_cnvs = JUMBLE_RUN.out.cns
        .map { meta, cns ->
            def sample_type = meta.sample_type == "tumor" ? 'somatic' : 'germline'
            [meta, cns, sample_type]
        }

    ANNOTATE_CNVS(
        ch_jumble_run_for_annotate_cnvs,
        ch_curation_ann,
    )

    emit:
    annotated_cns     = ANNOTATE_CNVS.out.cns            // channel: [ val(meta), path(cns) ]
    bins_csv          = JUMBLE_RUN.out.bins_csv          // channel: [ val(meta), path(csv) ]
    cnr               = JUMBLE_RUN.out.cnr               // channel: [ val(meta), path(cnr) ]
    cnv_plot_png      = JUMBLE_RUN.out.png               // channel: [ val(meta), path(png) ]
    jumble_cns        = JUMBLE_RUN.out.cns               // channel: [ val(meta), path(cns) ]
    profile_bedgraph  = JUMBLE_RUN.out.profile_bedgraph  // channel: [ val(meta), path(bedgraph) ]
    seg               = JUMBLE_RUN.out.seg               // channel: [ val(meta), path(seg) ]
    segments_bedgraph = JUMBLE_RUN.out.segments_bedgraph // channel: [ val(meta), path(bedgraph) ]

}

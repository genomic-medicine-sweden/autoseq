

include { JUMBLE_RUN                                  } from '../../../modules/local/jumble/main'
include { ANNOTATE_CNVS as ANNOTATE_CNVS_SOMATIC      } from '../../../modules/local/annotate_cnvs/main'
include { ANNOTATE_CNVS as ANNOTATE_CNVS_GERMLINE     } from '../../../modules/local/annotate_cnvs/main'


workflow CALL_CNVS {
    take:
    ch_input_bam            // channel: input BAM files from alignment workflow
    ch_jumble_ref           // channel: Jumble reference files
    ch_curation_ann         // channel: CNV curation/annotation file

    main:

    //
    // MODULE: JUMBLE_RUN to call CNVs
    //

    JUMBLE_RUN(
        ch_input_bam,
        ch_jumble_ref
    )

    //
    // Split CNS by sample type so each aliased ANNOTATE_CNVS run applies
    // the appropriate somatic/germline curation profile via ext.args.
    //

    def ch_jumble_cns_branched = JUMBLE_RUN.out.cns
        .branch { meta, _cns ->
            somatic:  meta.sample_type == "tumor"
            germline: true
        }

    //
    // MODULE: Annotate CNVs with cancer relevant genes
    //

    ANNOTATE_CNVS_SOMATIC(
        ch_jumble_cns_branched.somatic,
        ch_curation_ann
    )

    ANNOTATE_CNVS_GERMLINE(
        ch_jumble_cns_branched.germline,
        ch_curation_ann
    )

    ch_annotated_cns = ANNOTATE_CNVS_SOMATIC.out.cns.mix(ANNOTATE_CNVS_GERMLINE.out.cns)

    emit:
    jumble_cns          = JUMBLE_RUN.out.cns
    cnr                 = JUMBLE_RUN.out.cnr
    seg                 = JUMBLE_RUN.out.seg
    profile_bedgraph    = JUMBLE_RUN.out.profile_bedgraph
    segments_bedgraph   = JUMBLE_RUN.out.segments_bedgraph
    png                 = JUMBLE_RUN.out.png
    cns                 = ch_annotated_cns

}

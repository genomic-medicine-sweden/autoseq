
include { PURECN_RUN } from '../../../modules/local/purecn/run/main'
include { TYPEDPYD   } from '../../../modules/local/typeDPYD/main'


workflow PROFILE_TUMOR_BIOMARKERS {
    take:
    ch_cnr                  // channel: [mandatory] [ val(meta), path(cnr) ] jumble cnr
    ch_seg                  // channel: [mandatory] [ val(meta), path(seg) ] jumble seg
    ch_mutect2_vcf          // channel: [optional]  [ val(meta), path(vcf) ] unfiltered mutect2 vcf
    ch_bam_bai              // channel: [mandatory] [ val(meta), path(bam), path(bai) ]

    main:

    //
    // MODULE: purecn
    //

    // The VCF carries the case meta, so it is joined to the tumor CNR/SEG on `case_id`.
    // `remainder` keeps the tumor when no VCF is given, and the VCF-only entries are dropped
    def ch_purecn_input = ch_cnr
        .join(ch_seg)
        .map { meta, cnr, seg ->
            [meta.case_id, meta, cnr, seg]
        }
        .join(ch_mutect2_vcf.map { meta, vcf -> [meta.case_id, vcf] }, remainder: true)
        .filter { row -> row[1] != null }
        .map { _case_id, meta, cnr, seg, vcf ->
            [meta, cnr, seg, vcf ?: []]
        }

    def purecn_genome = params.genome.equals("GRCh37") ? 'hg19' : 'hg38'

    PURECN_RUN(
        ch_purecn_input,
        purecn_genome
    )

    //
    // MODULE: DPYD status
    //
    TYPEDPYD(ch_bam_bai)

    //
    // MODULE: MSI status
    //

    //
    // MODULE: FRANKENPLOT (Genomic Overview, HRD, etc.)
    //

    emit:
    purecn_csv          = PURECN_RUN.out.csv            // channel: [ val(meta), path(purecn_csv) ]
    purecn_genes_csv    = PURECN_RUN.out.genes_csv      // channel: [ val(meta), path(purecn_genes_csv) ]
    purecn_variants_csv = PURECN_RUN.out.variants_csv   // channel: [ val(meta), path(purecn_variants_csv) ]
    purecn_loh_csv      = PURECN_RUN.out.loh_csv        // channel: [ val(meta), path(purecn_loh_csv) ]
    purecn_pdf          = PURECN_RUN.out.pdf            // channel: [ val(meta), path(purecn_pdf) ]
    dpyd_csv            = TYPEDPYD.out.csv              // channel: [ val(meta), path(dpyd_csv) ]
    dpyd_json           = TYPEDPYD.out.json             // channel: [ val(meta), path(dpyd_json) ]

}


include { PURECN_RUN } from '../../../modules/local/purecn/run/main'
include { TYPEDPYD   } from '../../../modules/local/typeDPYD/main'


workflow PROFILE_TUMOR_BIOMARKERS {
    take:
    ch_cnr                  // channel: [ val(meta), path(cnr) ] jumble cnr
    ch_seg                  // channel: [ val(meta), path(seg) ] jumble seg
    ch_mutect2_vcf          // channel: [ val(meta), path(vcf) ] unfiltered mutect2 vcf
    ch_bam_bai              // channel: [ val(meta), path(bam), path(bai) ]

    main:

    //
    // MODULE: purecn
    //
    ch_purecn_input = ch_cnr
        .combine(ch_seg)
        .combine(ch_mutect2_vcf)
        .map { cnr_meta, cnr, _seg_meta, seg, _vcf_meta, vcf ->
            def meta = cnr_meta
            def tcnr = cnr
            def tseg = seg
            def tvcf = vcf

            return tuple(meta, tcnr, tseg, tvcf)
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

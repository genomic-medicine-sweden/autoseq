
include { PURECN_RUN } from '../../../modules/local/purecn/run/main'
include { TYPEDPYD   } from '../../../modules/local/typeDPYD/main'


workflow PROFILE_TUMOR_BIOMARKERS {
    take:
    ch_cnr
    ch_seg
    ch_mutect2_vcf
    ch_bam

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
    TYPEDPYD(ch_bam)

    //
    // MODULE: MSI status
    //

    //
    // MODULE: FRANKENPLOT (Genomic Overview, HRD, etc.)
    //

    emit:
    purecn_csv          = PURECN_RUN.out.csv
    purecn_genes_csv    = PURECN_RUN.out.genes_csv
    purecn_variants_csv = PURECN_RUN.out.variants_csv
    purecn_loh_csv      = PURECN_RUN.out.loh_csv
    purecn_pdf          = PURECN_RUN.out.pdf
    dpyd_csv            = TYPEDPYD.out.csv
    dpyd_json           = TYPEDPYD.out.json

}

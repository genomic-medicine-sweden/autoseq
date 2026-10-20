
include { JUMBLE_FRANKENPLOT } from '../../../modules/local/jumble/frankenplot/main'
include { PURECN_RUN         } from '../../../modules/local/purecn/run/main'
include { TYPEDPYD           } from '../../../modules/local/typeDPYD/main'
include { keyByCaseId        } from '../utils_nfcore_autoseq_pipeline'


workflow PROFILE_TUMOR_BIOMARKERS {
    take:
    ch_cnr                  // channel: [mandatory] [ val(meta), path(cnr) ] jumble cnr
    ch_seg                  // channel: [mandatory] [ val(meta), path(seg) ] jumble seg
    ch_mutect2_vcf          // channel: [mandatory] [ val(meta), path(vcf) ] unfiltered mutect2 vcf
    ch_bam_bai              // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_germline_vcf         // channel: [optional]  [ val(meta), path(vcf) ] VEP annotated germline vcf with tumor allele fraction
    ch_hetsnps_vcf          // channel: [optional]  [ val(meta), path(vcf) ] heterozygous SNPs with tumor allelic depths
    ch_jumble_cns           // channel: [mandatory] [ val(meta), path(cns) ] jumble cns
    ch_jumble_csv           // channel: [mandatory] [ val(meta), path(csv) ] jumble bin-level table
    ch_somatic_vcf          // channel: [mandatory]  [ val(meta), path(vcf) ] VEP annotated somatic vcf

    main:

    //
    // MODULE: purecn
    //

    // The meta from `ch_mutect2_vcf` carries the `case_id`, so this can be used to join it to the tumor [meta, cnr, seg] to ensure `ch_purecn_input` is separated per case.
    // The tumor meta is kept so PureCN outputs are named per tumor.
    def ch_purecn_input = ch_cnr
        .join(ch_seg)
        .map { meta, cnr, seg ->
            [meta.case_id, meta, cnr, seg]
        }
        .join(ch_mutect2_vcf.map { meta, vcf -> [meta.case_id, vcf] })
        .map { _case_id, meta, cnr, seg, vcf ->
            [meta, cnr, seg, vcf]
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

    // Every input is keyed on `case_id` and its files are packed into a single element, so that
    // `join(remainder: true)` pads a missing optional input with one `null`. The tumor meta is
    // carried downstream because the report describes the tumor
    def ch_jumble = ch_jumble_cns.join(ch_jumble_csv)
    def ch_dpyd   = TYPEDPYD.out.json.join(TYPEDPYD.out.csv)

    def ch_tumor_jumble  = ch_jumble.filter { meta, _cns, _csv -> meta.sample_type == "tumor" }
    def ch_normal_jumble = ch_jumble.filter { meta, _cns, _csv -> meta.sample_type == "normal" }
    // DPYD is a germline marker, so the normal typing is reported
    def ch_normal_dpyd   = ch_dpyd.filter { meta, _json, _csv -> meta.sample_type == "normal" }

    // The somatic VCF is always produced (Mutect2 runs for every case). The normal Jumble output,
    // germline VCF, het SNPs and normal DPYD come from the normal sample, so they are missing in
    // tumor-only runs
    def ch_frankenplot_input = keyByCaseId(ch_tumor_jumble, true)
        .join(keyByCaseId(ch_normal_jumble), remainder: true)
        .join(keyByCaseId(ch_somatic_vcf))
        .join(keyByCaseId(ch_germline_vcf), remainder: true)
        .join(keyByCaseId(ch_hetsnps_vcf), remainder: true)
        .join(keyByCaseId(ch_normal_dpyd), remainder: true)
        .filter { row -> row[1] != null }
        .map { _case_id, tumor, normal, somatic_vcf, germline_vcf, hetsnps_vcf, dpyd ->
            def (meta, tumor_cns, tumor_csv) = tumor
            def (normal_cns, normal_csv)     = normal ?: [[], []]
            def (dpyd_json, dpyd_csv)        = dpyd ?: [[], []]
            [
                meta,
                tumor_cns,
                normal_cns,
                tumor_csv,
                normal_csv,
                somatic_vcf,
                germline_vcf ?: [],
                hetsnps_vcf ?: [],
                dpyd_json,
                dpyd_csv
            ]
        }

    JUMBLE_FRANKENPLOT(ch_frankenplot_input)

    emit:
    purecn_csv          = PURECN_RUN.out.csv            // channel: [ val(meta), path(purecn_csv) ]
    purecn_genes_csv    = PURECN_RUN.out.genes_csv      // channel: [ val(meta), path(purecn_genes_csv) ]
    purecn_variants_csv = PURECN_RUN.out.variants_csv   // channel: [ val(meta), path(purecn_variants_csv) ]
    purecn_loh_csv      = PURECN_RUN.out.loh_csv        // channel: [ val(meta), path(purecn_loh_csv) ]
    purecn_pdf          = PURECN_RUN.out.pdf            // channel: [ val(meta), path(purecn_pdf) ]
    dpyd_csv            = TYPEDPYD.out.csv              // channel: [ val(meta), path(dpyd_csv) ]
    dpyd_json           = TYPEDPYD.out.json             // channel: [ val(meta), path(dpyd_json) ]
    frankenplot_html    = JUMBLE_FRANKENPLOT.out.html   // channel: [ val(meta), path(html) ]

}

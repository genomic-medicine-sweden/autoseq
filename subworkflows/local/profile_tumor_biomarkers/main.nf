
include { JUMBLE_FRANKENPLOT } from '../../../modules/local/jumble/frankenplot/main'
include { PURECN_RUN         } from '../../../modules/local/purecn/run/main'
include { TYPEDPYD           } from '../../../modules/local/typeDPYD/main'


workflow PROFILE_TUMOR_BIOMARKERS {
    take:
    ch_cnr                  // channel: [mandatory] [ val(meta), path(cnr) ] jumble cnr
    ch_seg                  // channel: [mandatory] [ val(meta), path(seg) ] jumble seg
    ch_mutect2_vcf          // channel: [optional]  [ val(meta), path(vcf) ] unfiltered mutect2 vcf
    ch_bam_bai              // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_jumble_cns           // channel: [mandatory] [ val(meta), path(cns) ] jumble cns
    ch_jumble_csv           // channel: [mandatory] [ val(meta), path(csv) ] jumble bin-level table of tumor and normal samples
    ch_somatic_vcf          // channel: [optional]  [ val(meta), path(vcf) ] VEP annotated somatic vcf
    ch_germline_vcf         // channel: [optional]  [ val(meta), path(vcf) ] VEP annotated germline vcf with tumor allele fraction
    ch_hetsnps_vcf          // channel: [optional]  [ val(meta), path(vcf) ] heterozygous SNPs with tumor allelic depths

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

    // Every input is keyed on `case_id` and its files are packed into a single element, so that
    // `join(remainder: true)` pads a missing optional input with one `null`. The tumor meta is
    // carried downstream because the report describes the tumor
    def ch_jumble = ch_jumble_cns
        .join(ch_jumble_csv)
        .branch { meta, _cns, _csv ->
            tumor: meta.sample_type == "tumor"
            normal: meta.sample_type == "normal"
        }

    def ch_tumor_jumble_by_case = ch_jumble.tumor
        .map { meta, cns, csv ->
            return [meta.case_id, [meta, cns, csv]]
        }

    def ch_normal_jumble_by_case = ch_jumble.normal
        .map { meta, cns, csv ->
            return [meta.case_id, [cns, csv]]
        }

    // DPYD is a germline marker, so the normal typing is reported
    def ch_normal_dpyd_by_case = TYPEDPYD.out.json
        .join(TYPEDPYD.out.csv)
        .filter { meta, _json, _csv -> meta.sample_type == "normal" }
        .map { meta, json, csv ->
            return [meta.case_id, [json, csv]]
        }

    def ch_frankenplot_input = ch_tumor_jumble_by_case
        .join(ch_normal_jumble_by_case, remainder: true)
        .join(ch_somatic_vcf.map { meta, vcf -> [meta.case_id, vcf] }, remainder: true)
        .join(ch_germline_vcf.map { meta, vcf -> [meta.case_id, vcf] }, remainder: true)
        .join(ch_hetsnps_vcf.map { meta, vcf -> [meta.case_id, vcf] }, remainder: true)
        .join(ch_normal_dpyd_by_case, remainder: true)
        .filter { row -> row[1] != null }
        .map { _case_id, tumor, normal, somatic_vcf, germline_vcf, hetsnps_vcf, dpyd ->
            def (meta, tumor_cns, tumor_csv) = tumor
            def (normal_cns, normal_csv)     = normal ?: [[], []]
            def (dpyd_json, dpyd_csv)        = dpyd ?: [[], []]
            return [
                meta,
                tumor_cns,
                normal_cns,
                tumor_csv,
                normal_csv,
                somatic_vcf ?: [],
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

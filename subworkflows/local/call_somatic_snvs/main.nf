//
// Subworkflow for somatic SNV calling
//

include { BAM_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING_GATK as CALL_GATK_MUTECT2 } from '../../../subworkflows/nf-core/bam_tumor_normal_somatic_variant_calling_gatk/main'
include { BCFTOOLS_FILTER as PASSFILTER_FOR_MUTECT2                         } from '../../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_FILTER as PASSFILTER_FOR_SAGE                            } from '../../../modules/nf-core/bcftools/filter/main'
include { ENSEMBLVEP_VEP as ANNOTATE_VEP                                    } from '../../../modules/nf-core/ensemblvep/vep/main'
include { SAGE_SOMATIC                                                      } from '../../../modules/local/sage/somatic/main'
include { SOMATIC_VCFMERGE                                                  } from '../../../modules/local/vcfmerge/main'
include { TABIX_TABIX as INDEX_VCF                                          } from '../../../modules/nf-core/tabix/tabix/main'
include { VT_DECOMPOSE                                                      } from '../../../modules/nf-core/vt/decompose/main'
include { VT_NORMALIZE                                                      } from '../../../modules/nf-core/vt/normalize/main'

workflow CALL_SOMATIC_SNVS {
    take:
    ch_bam_bai                     // channel: [mandatory] [ val(meta), [ path(tumor_bam), path(normal_bam) ], [ path(tumor_bai), path(normal_bai) ], path(intervals) ]
    ch_dict                        // channel: [mandatory] [ val(meta), path(dict) ]
    ch_ensembl_data_resources      // channel: [mandatory] [ val(meta), path(dir) ]
    ch_ensembl_vep_cache           // channel: [mandatory] [ val(meta), path(cache) ]
    ch_fai                         // channel: [mandatory] [ val(meta), path(fai) ]
    ch_fasta                       // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_germline_resource           // channel: [mandatory] [ path(vcf) ]
    ch_germline_resource_tbi       // channel: [mandatory] [ path(tbi) ]
    ch_interval_file               // channel: [mandatory] [ path(interval_list) ]
    ch_panel_of_normals            // channel: [optional]  [ path(vcf) ]
    ch_panel_of_normals_tbi        // channel: [optional]  [ path(tbi) ]
    ch_sage_highconf_regions       // channel: [mandatory] [ val(meta), path(bed) ]
    ch_sage_known_hotspots_somatic // channel: [mandatory] [ val(meta), path(vcf) ]
    ch_sage_pon                    // channel: [mandatory] [ val(meta), path(tsv) ]
    val_ensemblvep_version         // string:  [mandatory] Ensembl VEP cache version
    val_genome                     // string:  [mandatory] genome assembly (e.g. "GRCh37")
    val_genome_version             // string:  [mandatory] numeric genome version for SAGE ("37" or "38")
    val_min_baseq                  // integer: [mandatory] minimum average base quality for SAGE
    val_sage_hotspot_tumor_qual    // integer: [mandatory] SAGE minimum tumor quality for hotspot variants
    val_sage_min_map_quality       // integer: [mandatory] SAGE minimum mapping quality
    val_sage_min_tumor_qual        // integer: [mandatory] SAGE minimum tumor quality for panel variants
    val_sage_min_tumor_vaf         // float:   [mandatory] SAGE minimum tumor VAF
    val_vep_species                // string:  [mandatory] Ensembl VEP species (e.g. "homo_sapiens")

    main:
    ch_versions = channel.empty()

    //
    // SUBWORKFLOW: Call somatic SNVs with GATK Mutect2 in tumor-normal mode
    //

    // GATK4_MUTECT2 expects the index tuple as [ val(meta), path(fai), path(gzi) ]
    ch_fai_for_call_gatk_mutect2 = ch_fai.map { meta, fai ->
        return tuple(meta, fai, [])
    }

    CALL_GATK_MUTECT2(
        ch_bam_bai,
        ch_fasta,
        ch_fai_for_call_gatk_mutect2,
        ch_dict,
        [], // alleles
        [], // alleles_tbi
        ch_germline_resource,
        ch_germline_resource_tbi,
        ch_panel_of_normals,
        ch_panel_of_normals_tbi,
        ch_interval_file
    )

    //
    // MODULE: Decompose and normalize the filtered Mutect2 calls
    //

    ch_call_gatk_mutect2_for_vt_decompose = CALL_GATK_MUTECT2.out.filtered_vcf
        .map { meta, vcf ->
            return tuple(meta, vcf, []) // empty intervals
        }

    VT_DECOMPOSE(ch_call_gatk_mutect2_for_vt_decompose)

    ch_vt_decompose_for_vt_normalize = VT_DECOMPOSE.out.vcf
        .map { meta, vcf ->
            return tuple(meta, vcf, [], [])
        }

    VT_NORMALIZE(ch_vt_decompose_for_vt_normalize, ch_fasta, ch_fai)
    INDEX_VCF(VT_NORMALIZE.out.vcf)

    ch_vt_normalize_for_passfilter_for_mutect2 = VT_NORMALIZE.out.vcf
        .join(INDEX_VCF.out.tbi)

    //
    // MODULE: Call somatic SNVs with SAGE in tumor-normal mode
    //

    SAGE_SOMATIC(
        ch_bam_bai,
        ch_fasta,
        ch_fai,
        ch_dict,
        val_genome_version,
        ch_sage_pon.collect{ it -> it[1] },
        ch_sage_known_hotspots_somatic.collect{ it -> it[1] },
        ch_sage_highconf_regions.collect{ it -> it[1] },
        ch_ensembl_data_resources.collect{ it -> it[1] },
        true, // targeted_mode
        val_sage_min_tumor_vaf,
        val_sage_min_tumor_qual,
        val_sage_min_map_quality,
        val_sage_hotspot_tumor_qual,
        val_min_baseq
    )

    ch_sage_somatic_for_passfilter_for_sage = SAGE_SOMATIC.out.vcf
        .join(SAGE_SOMATIC.out.tbi)

    //
    // MODULE: Keep PASS variants and merge the Mutect2 and SAGE call sets
    //

    PASSFILTER_FOR_MUTECT2(ch_vt_normalize_for_passfilter_for_mutect2)
    PASSFILTER_FOR_SAGE(ch_sage_somatic_for_passfilter_for_sage)

    SOMATIC_VCFMERGE(
        PASSFILTER_FOR_MUTECT2.out.vcf.join(PASSFILTER_FOR_MUTECT2.out.index),
        PASSFILTER_FOR_SAGE.out.vcf
    )

    ch_versions = ch_versions.mix(SOMATIC_VCFMERGE.out.versions)

    //
    // MODULE: Annotate somatic variants with Ensembl VEP
    //

    ch_somatic_vcfmerge_for_annotate_vep = SOMATIC_VCFMERGE.out.vcf
        .map { meta, vcf ->
            return tuple(meta, vcf, []) // empty channel for custom files
        }

    ANNOTATE_VEP(
        ch_somatic_vcfmerge_for_annotate_vep,
        val_genome,
        val_vep_species,
        val_ensemblvep_version,
        ch_ensembl_vep_cache.collect(),
        ch_fasta,
        [],
        [[], []] // no GTF, annotate from the cache
    )

    emit:
    contamination_table    = CALL_GATK_MUTECT2.out.contamination_table // channel: [ val(meta), path(table) ]
    mutect2_filtered_tbi   = CALL_GATK_MUTECT2.out.filtered_tbi        // channel: [ val(meta), path(tbi) ]
    mutect2_filtered_vcf   = CALL_GATK_MUTECT2.out.filtered_vcf        // channel: [ val(meta), path(vcf) ]
    mutect2_stats          = CALL_GATK_MUTECT2.out.filtered_stats      // channel: [ val(meta), path(stats) ]
    mutect2_unfiltered_tbi = CALL_GATK_MUTECT2.out.mutect2_tbi         // channel: [ val(meta), path(tbi) ]
    mutect2_unfiltered_vcf = CALL_GATK_MUTECT2.out.mutect2_vcf         // channel: [ val(meta), path(vcf) ]
    pileup_table_normal    = CALL_GATK_MUTECT2.out.pileup_table_normal // channel: [ val(meta), path(table) ]
    pileup_table_tumor     = CALL_GATK_MUTECT2.out.pileup_table_tumor  // channel: [ val(meta), path(table) ]
    sage_tbi               = SAGE_SOMATIC.out.tbi                      // channel: [ val(meta), path(tbi) ]
    sage_vcf               = SAGE_SOMATIC.out.vcf                      // channel: [ val(meta), path(vcf) ]
    segmentation_table     = CALL_GATK_MUTECT2.out.segmentation_table  // channel: [ val(meta), path(table) ]
    somatic_tbi            = SOMATIC_VCFMERGE.out.tbi                  // channel: [ val(meta), path(tbi) ]
    somatic_vcf            = SOMATIC_VCFMERGE.out.vcf                  // channel: [ val(meta), path(vcf) ] merged Mutect2 + SAGE calls
    vep_tbi                = ANNOTATE_VEP.out.tbi                      // channel: [ val(meta), path(tbi) ]
    vep_vcf                = ANNOTATE_VEP.out.vcf                      // channel: [ val(meta), path(vcf) ]
    versions               = ch_versions                               // channel: [ path(versions.yml) ]
}

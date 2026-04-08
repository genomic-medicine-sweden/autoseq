
include { BAM_TUMOR_NORMAL_SOMATIC_VARIANT_CALLING_GATK as CALL_GATK_MUTECT2 } from '../../../subworkflows/nf-core/bam_tumor_normal_somatic_variant_calling_gatk/main'

include { VT_DECOMPOSE                                  } from '../../../modules/nf-core/vt/decompose/main'
include { VT_NORMALIZE                                  } from '../../../modules/nf-core/vt/normalize/main'
include { BCFTOOLS_FILTER as PASSFILTER_FOR_MUTECT2     } from '../../../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_FILTER as PASSFILTER_FOR_SAGE        } from '../../../modules/nf-core/bcftools/filter/main'
include { SAGE_SOMATIC                                  } from '../../../modules/local/sage/somatic/main'
include { SOMATIC_VCFMERGE                              } from '../../../modules/local/vcfmerge/main'
include { TABIX_TABIX  as INDEX_VCF                     } from '../../../modules/nf-core/tabix/tabix/main'
include { ENSEMBLVEP_VEP as ANNOTATE_VEP                } from '../../../modules/nf-core/ensemblvep/vep'


workflow CALL_SOMATIC_SNVS {
    take:
    ch_input                 // channel: [ val(meta), path(input), path(input_index), val(which_norm) ]
    ch_fasta                 // channel: [ val(meta), path(fasta) ]
    ch_fai                   // channel: [ val(meta), path(fai) ]
    ch_dict                  // channel: [ val(meta), path(dict) ]
    ch_germline_resource     // channel: /path/to/germline/resource
    ch_germline_resource_tbi // channel: /path/to/germline/index
    ch_panel_of_normals      // channel: /path/to/panel/of/normals
    ch_panel_of_normals_tbi  // channel: /path/to/panel/of/normals/index
    ch_interval_file         // channel: /path/to/interval/file
    ch_sage_known_hotspots_somatic
    ch_sage_highconf_regions
    ch_sage_pon
    ch_ensembl_data_resources
    ch_ensembl_vep_cache
    genome_version

    main:
    versions = Channel.empty()

    // Call somatic SNVs using GATK Mutect2 in tumor-normal mode
    CALL_GATK_MUTECT2(
        ch_input,
        ch_fasta,
        ch_fai,
        ch_dict,
        ch_germline_resource,
        ch_germline_resource_tbi,
        ch_panel_of_normals,
        ch_panel_of_normals_tbi,
        ch_interval_file
    )

    ch_mutect2_vcf = CALL_GATK_MUTECT2.out.filtered_vcf
        .map { meta, vcf ->
            return tuple(meta, vcf, []) // empty intervals
        }

    VT_DECOMPOSE (ch_mutect2_vcf)

    ch_vtdecompose_vcf = VT_DECOMPOSE.out.vcf
        .map { meta, vcf ->
            return tuple(meta, vcf, [], [])
        }

    VT_NORMALIZE (ch_vtdecompose_vcf, ch_fasta, ch_fai)
    INDEX_VCF (VT_NORMALIZE.out.vcf)

    ch_mutect2_vcf = VT_NORMALIZE.out.vcf
        .join(INDEX_VCF.out.tbi)


    // Call somatic SNVs using SAGE in tumor-normal mode
    SAGE_SOMATIC(
        ch_input,
        ch_fasta,
        ch_fai,
        ch_dict,
        genome_version, // genome version
        ch_sage_pon.collect{it[1]},
        ch_sage_known_hotspots_somatic.collect{it[1]},
        ch_sage_highconf_regions.collect{it[1]},
        ch_ensembl_data_resources.collect{it[1]},
        true, // targeted_mode,
        params.sage_min_tumor_vaf,
        params.sage_min_tumor_qual,
        params.sage_min_map_quality,
        params.sage_hotspot_tumor_qual,
        params.min_baseq
    )

    ch_sage_vcf = SAGE_SOMATIC.out.vcf
            .join(SAGE_SOMATIC.out.tbi)

    PASSFILTER_FOR_MUTECT2 (ch_mutect2_vcf)
    PASSFILTER_FOR_SAGE (ch_sage_vcf)

    SOMATIC_VCFMERGE (
        PASSFILTER_FOR_MUTECT2.out.vcf,
        PASSFILTER_FOR_SAGE.out.vcf
    )

    ch_input_vcf = SOMATIC_VCFMERGE.out.vcf
        .map { meta, vcf ->
            return tuple(meta, vcf, []) // empty channel for custom files
        }

    // VEP Annotation
    ANNOTATE_VEP (
        ch_input_vcf,
        params.genome,
        "homo_sapiens",
        params.ensemblvep_version,
        ch_ensembl_vep_cache.collect{it -> it[1]},
        ch_fasta,
        []
    )


    versions = versions.mix(CALL_GATK_MUTECT2.out.versions)
    versions = versions.mix(VT_DECOMPOSE.out.versions)
    versions = versions.mix(VT_NORMALIZE.out.versions)
    versions = versions.mix(SAGE_SOMATIC.out.versions)
    versions = versions.mix(PASSFILTER_FOR_MUTECT2.out.versions)
    versions = versions.mix(PASSFILTER_FOR_SAGE.out.versions)
    versions = versions.mix(SOMATIC_VCFMERGE.out.versions)
    versions = versions.mix(ANNOTATE_VEP.out.versions)


    emit:
    contamination_table = CALL_GATK_MUTECT2.out.contamination_table    // channel: [ val(meta), path(table) ]
    mutect2_stats       = CALL_GATK_MUTECT2.out.filtered_stats                 // channel: [ val(meta), path(stats) ]
    mutect2_tbi         = CALL_GATK_MUTECT2.out.filtered_tbi                   // channel: [ val(meta), path(tbi) ]
    mutect2_vcf         = CALL_GATK_MUTECT2.out.filtered_vcf                   // channel: [ val(meta), path(vcf) ]
    pileup_table_normal = CALL_GATK_MUTECT2.out.pileup_table_normal         // channel: [ val(meta), path(table) ]
    pileup_table_tumor  = CALL_GATK_MUTECT2.out.pileup_table_tumor          // channel: [ val(meta), path(table) ]
    sage_vcf            = SAGE_SOMATIC.out.vcf                                 // channel: [ val(meta), path(vcf) ]
    sage_tbi            = SAGE_SOMATIC.out.tbi
    somatic_vcf         = SOMATIC_VCFMERGE.out.vcf
    somatic_tbi         = SOMATIC_VCFMERGE.out.tbi
    vep_vcf             = ANNOTATE_VEP.out.vcf
    vep_tbi             = ANNOTATE_VEP.out.tbi
    versions
}

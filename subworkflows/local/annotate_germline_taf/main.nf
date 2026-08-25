include { GATK4_HAPLOTYPECALLER } from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_GENOTYPEGVCFS   } from '../../../modules/nf-core/gatk4/genotypegvcfs/main'
include { BCFTOOLS_MERGE        } from '../../../modules/nf-core/bcftools/merge/main'


workflow ANNOTATE_GERMLINE_TAF {
    take:
    ch_tbam_gvcf        // channel: [mandatory] [ val(meta), path(tbam), path(tbai), path(gvcf), path(gvcf_tbi) ]
    ch_genome_fasta     // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_genome_fai       // channel: [mandatory] [ val(meta), path(fai) ]
    ch_genome_dict      // channel: [mandatory] [ val(meta), path(dict) ]

    main:

    // The germline VCF is passed through the `intervals` slot so that the tumour is
    // genotyped only at the sites where a germline variant was called
    def ch_tbam_gvcf_for_haplotypecaller = ch_tbam_gvcf
        .map { meta, tbam, tbai, gvcf, _gvcf_tbi ->
            return [meta, tbam, tbai, gvcf, []]
        }

    // GATK refuses to read a block-compressed VCF as intervals unless its index sits next to
    // it, and the module has no slot for one. `dbsnp_tbi` is only ever staged, never
    // referenced in the command, so it is used to bring the index into the work directory
    def ch_gvcf_tbi_for_haplotypecaller = ch_tbam_gvcf
        .map { meta, _tbam, _tbai, _gvcf, gvcf_tbi ->
            return [meta, gvcf_tbi]
        }

    GATK4_HAPLOTYPECALLER (
        ch_tbam_gvcf_for_haplotypecaller,
        ch_genome_fasta,
        ch_genome_fai,
        ch_genome_dict,
        [[], []],
        ch_gvcf_tbi_for_haplotypecaller
    )

    def ch_haplotypecaller_for_genotypegvcfs = ch_tbam_gvcf
        .join(GATK4_HAPLOTYPECALLER.out.vcf)
        .join(GATK4_HAPLOTYPECALLER.out.tbi)
        .map { meta, _tbam, _tbai, gvcf, gvcf_tbi, vcf, tbi ->
            return [meta, vcf, tbi, gvcf, gvcf_tbi]
        }

    GATK4_GENOTYPEGVCFS (
        ch_haplotypecaller_for_genotypegvcfs,
        ch_genome_fasta,
        ch_genome_fai,
        ch_genome_dict,
        [[], []],
        [[], []]
    )

    // Merging the tumour VCF with the germline VCF puts the tumour genotype next to the
    // germline one, so every germline variant carries the tumour allele fraction
    def ch_genotypegvcfs_for_bcftools_merge = ch_tbam_gvcf
        .join(GATK4_GENOTYPEGVCFS.out.vcf)
        .join(GATK4_GENOTYPEGVCFS.out.tbi)
        .map { meta, _tbam, _tbai, gvcf, gvcf_tbi, vcf, tbi ->
            return [meta, [gvcf, vcf], [gvcf_tbi, tbi], []]
        }

    // `.collect()` restores the value-channel semantics lost by `combine`, so that the
    // reference can be reused by every sample
    def ch_fasta_fai_for_bcftools_merge = ch_genome_fasta
        .combine(ch_genome_fai)
        .map { meta, fasta, _meta_fai, fai ->
            return [meta, fasta, fai]
        }
        .collect()

    BCFTOOLS_MERGE (
        ch_genotypegvcfs_for_bcftools_merge,
        ch_fasta_fai_for_bcftools_merge
    )

    emit:
    germline_taf_vcf = BCFTOOLS_MERGE.out.vcf   // channel: [ val(meta), path(vcf) ]
    germline_taf_tbi = BCFTOOLS_MERGE.out.index // channel: [ val(meta), path(tbi) ]

}

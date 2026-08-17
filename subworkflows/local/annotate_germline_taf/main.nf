include { GATK4_HAPLOTYPECALLER } from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_GENOTYPEGVCFS   } from '../../../modules/nf-core/gatk4/genotypegvcfs/main'
include { BCFTOOLS_MERGE        } from '../../../modules/nf-core/bcftools/merge/main'


workflow ANNOTATE_GERMLINE_TAF {
    take:
    ch_tbam_gvcf        // channel: [mandatory] input channel of tumor bam and germline vep vcf files
    ch_genome_fasta     // channel: [mandatory] input channel of reference genome fasta file
    ch_genome_fai       // channel: [mandatory] input channel of reference genome fai file
    ch_genome_dict      // channel: [mandatory] input channel of reference genome dict file

    main:

    def ch_input_haplotypecaller = ch_tbam_gvcf
        .map { meta, tbam, tbai, gvcf, gvcf_tbi ->
            return [meta, tbam, tbai, gvcf, []]
        }

    GATK4_HAPLOTYPECALLER (
        ch_input_haplotypecaller,
        ch_genome_fasta,
        ch_genome_fai,
        ch_genome_dict,
        [],
        []
    )

    def ch_input_genotypegvcfs = ch_tbam_gvcf
        .join(GATK4_HAPLOTYPECALLER.out.vcf)
        .join(GATK4_HAPLOTYPECALLER.out.tbi)
        .map { meta, tbam, tbai, gvcf, gvcf_tbi, dragstr_model, vcf, tbi ->
            return [meta, vcf, tbi, gvcf, []]
        }

    GATK4_GENOTYPEGVCFS (
        ch_input_genotypegvcfs,
        ch_genome_fasta,
        ch_genome_fai,
        ch_genome_dict,
        [],
        []
    )

    def ch_input_bcftools_merge = ch_tbam_gvcf
        .join(GATK4_GENOTYPEGVCFS.out.vcf)
        .join(GATK4_GENOTYPEGVCFS.out.tbi)
        .map { meta, tbam, tbai, gvcf, gvcf_tbi, dragstr_model, vcf, tbi ->
            return [meta, [gvcf, vcf], [gvcf_tbi, tbi], []]
        }

    def ch_fasta_fai = ch_genome_fasta
        .join(ch_genome_fai)
        .map { fasta, fai ->
            return [fasta, fai]
        }

    BCFTOOLS_MERGE (
        ch_input_bcftools_merge,
        ch_fasta_fai,
    )



    emit:
    germline_taf_vcf = BCFTOOLS_MERGE.out.vcf
    germline_taf_tbi = BCFTOOLS_MERGE.out.index

}


include { GATK4_HAPLOTYPECALLER     } from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { ENSEMBLVEP_VEP            } from '../../../modules/nf-core/ensemblvep/vep/main'

workflow CALL_GERMLINE_SNVS {
    take:
        ch_input                // [meta, bam, bai, intervals, dragstr_model]
        ch_genome_fasta         // [meta_ref, fasta]
        ch_genome_fai           // [meta_ref, fai]
        ch_dict                 // [meta_ref, dict]
        ch_vep_cache            // [meta_cache, cache_dir]
        ch_dbsnp_vcf
        ch_dbsnp_vcf_tbi

    main:

        GATK4_HAPLOTYPECALLER(
            ch_input,
            ch_genome_fasta,
            ch_genome_fai,
            ch_dict,
            ch_dbsnp_vcf,
            ch_dbsnp_vcf_tbi
        )

        haplotypecaller_vcf = GATK4_HAPLOTYPECALLER.out.vcf
            .map { meta, vcf ->
                [meta, vcf, []]  // Add empty custom_extra_files for VEP
            }

        ENSEMBLVEP_VEP(
            haplotypecaller_vcf,
            params.genome,
            params.vep_species,  // Assuming human genome; adjust as needed
            params.ensemblvep_version,
            ch_vep_cache.collect{it -> it[1]},
            ch_genome_fasta,
            []
        )

    emit:
    vcf     = GATK4_HAPLOTYPECALLER.out.vcf
    tbi     = GATK4_HAPLOTYPECALLER.out.tbi
    vep_vcf = ENSEMBLVEP_VEP.out.vcf
    vep_tbi = ENSEMBLVEP_VEP.out.tbi

}


include { GATK4_HAPLOTYPECALLER     } from '../../../modules/nf-core/gatk4/haplotypecaller/main'
include { ENSEMBLVEP_VEP            } from '../../../modules/nf-core/ensemblvep/vep/main'

workflow CALL_GERMLINE_SNVS {
    take:
        ch_bam                  // [meta, bam, bai]
        ch_genome_fasta         // [meta_ref, fasta]
        ch_genome_fai           // [meta_ref, fai]
        ch_dict                 // [meta_ref, dict]
        ch_vep_cache            // [meta_cache, cache_dir]

    main:

        GATK4_HAPLOTYPECALLER(
            ch_bam,
            ch_genome_fasta,
            ch_genome_fai,
            ch_dict
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
    vep_vcf = ENSEMBLVEP_VEP.out.vcf
}

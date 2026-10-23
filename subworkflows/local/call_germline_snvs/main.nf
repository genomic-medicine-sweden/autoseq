//
// Subworkflow for germline SNV calling
//

include { ENSEMBLVEP_VEP        } from '../../../modules/nf-core/ensemblvep/vep/main'
include { GATK4_HAPLOTYPECALLER } from '../../../modules/nf-core/gatk4/haplotypecaller/main'

workflow CALL_GERMLINE_SNVS {
    take:
    ch_bam_bai             // channel: [mandatory] [ val(meta), path(bam), path(bai), path(intervals), path(dragstr_model) ]
    ch_dbsnp_vcf           // channel: [optional]  [ val(meta), path(vcf) ]
    ch_dbsnp_vcf_tbi       // channel: [optional]  [ val(meta), path(tbi) ]
    ch_dict                // channel: [mandatory] [ val(meta), path(dict) ]
    ch_genome_fai          // channel: [mandatory] [ val(meta), path(fai) ]
    ch_genome_fasta        // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_vep_cache           // channel: [mandatory] [ val(meta), path(cache) ]
    val_ensemblvep_version // string:  [mandatory] Ensembl VEP cache version
    val_genome             // string:  [mandatory] genome assembly (e.g. "GRCh37")
    val_vep_species        // string:  [mandatory] Ensembl VEP species (e.g. "homo_sapiens")

    main:

    //
    // MODULE: Call germline SNVs and indels with GATK HaplotypeCaller
    //

    GATK4_HAPLOTYPECALLER(
        ch_bam_bai,
        ch_genome_fasta,
        ch_genome_fai,
        ch_dict,
        ch_dbsnp_vcf,
        ch_dbsnp_vcf_tbi,
    )

    ch_haplotypecaller_for_ensemblvep_vep = GATK4_HAPLOTYPECALLER.out.vcf
        .map { meta, vcf ->
            [meta, vcf, []]  // Add empty custom_extra_files for VEP
        }

    //
    // MODULE: Annotate germline variants with Ensembl VEP
    //

    ENSEMBLVEP_VEP(
        ch_haplotypecaller_for_ensemblvep_vep,
        val_genome,
        val_vep_species,
        val_ensemblvep_version,
        ch_vep_cache.collect(),
        ch_genome_fasta,
        [],
        [[], []] // no GTF, annotate from the cache
    )

    emit:
    haplotypecaller_tbi = GATK4_HAPLOTYPECALLER.out.tbi // channel: [ val(meta), path(tbi) ]
    haplotypecaller_vcf = GATK4_HAPLOTYPECALLER.out.vcf // channel: [ val(meta), path(vcf) ]
    vep_tbi             = ENSEMBLVEP_VEP.out.tbi        // channel: [ val(meta), path(tbi) ]
    vep_vcf             = ENSEMBLVEP_VEP.out.vcf        // channel: [ val(meta), path(vcf) ]

}

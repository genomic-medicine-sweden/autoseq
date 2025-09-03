//
//
//

include { BWAMEM2_MEM } from '../../../modules/nf-core/bwamem2/mem/main'
include { GATK4_MARKDUPLICATES as MARKDUPLICATES } from '../../../modules/nf-core/gatk4/markduplicates/main'


workflow ALIGNMENT {
    take:
    ch_input_reads // channel: input reads from samplesheet
    ch_genome_fasta
    ch_genome_fai
    ch_bwamem2_index

    main:

    ch_versions  = Channel.empty()

    //
    // MODULE: Run BWA-MEM2 alignment
    //

    ch_fastq_reads = ch_input_reads
        .map { meta, reads ->
            meta  = meta + [read_group: "${meta.case_id}.${meta.sample_name}.${meta.lane}".toString(), split: null]

            return [meta, reads[0], reads[1]]
        }

    BWAMEM2_MEM(
        ch_fastq_reads,
        ch_bwamem2_index,
        ch_genome_fasta,
        true // sort_bam
    )

    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())


    ch_input_bam = BWAMEM2_MEM.out.bam
        .map { meta, bam, bai ->
            def id = "${meta.case_id}_${meta.sample_name}".toString()
            tuple(id, bam)
        }
        .groupTuple()


    ch_bam_meta = BWAMEM2_MEM.out.bam
        .map { meta, bam, bai ->
            def id = "${meta.case_id}_${meta.sample_name}".toString()
            tuple(id, meta.case_id, meta.sample_name, meta.sample_type)
        }
        .groupTuple()
        .map { id, case_id, sample_name, sample_type ->
            def meta = [
                id         : id,
                case_id    : case_id.unique()[0],
                sample_name: sample_name.unique()[0],
                sample_type: sample_type.unique()[0]
            ]

            tuple(id, meta)
        }
        .combine(ch_input_bam, by: 0)

    ch_input_bam = ch_bam_meta
        .map { id, meta, bam ->
            tuple(meta, bam)
        }

    MARKDUPLICATES (
        ch_input_bam,
        ch_genome_fasta.collect{it[1]},
        ch_genome_fai.collect{it[1]}
    )

    ch_versions = ch_versions.mix(MARKDUPLICATES.out.versions.first())


    emit:
    bam        = MARKDUPLICATES.out.bam
    metrics    = MARKDUPLICATES.out.metrics
    versions   = ch_versions            // channel: versions.yml
}

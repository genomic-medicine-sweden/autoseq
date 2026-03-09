//
//
//

include { BWAMEM2_MEM } from '../../../modules/nf-core/bwamem2/mem/main'
include { GATK4_MARKDUPLICATES as MARKDUPLICATES } from '../../../modules/nf-core/gatk4/markduplicates/main'


workflow READ_ALIGNMENT {
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

    ch_input_fastqs = ch_input_reads
        .map { meta, reads ->
            tuple(meta.sample_name, meta, reads)
        }
        .groupTuple()
        .flatMap { _sample_name, meta_list, reads_list ->
            def processed_reads = []
            meta_list.eachWithIndex { meta, index ->

                def seq_index = index + 1
                def new_meta = meta + [
                    id         : "${meta.sample_name}.${seq_index}".toString(),
                    read_group : "\"@RG\\tID:${meta.sample_name}_${seq_index}\\tSM:${meta.sample_name}\\tLB:${meta.sample_name}\\tPL:ILLUMINA\"".toString(),
                    split      : null
                ]

                processed_reads.add( tuple(new_meta, reads_list[index]) )
            }

            return processed_reads
        }


    sort_bam = true
    BWAMEM2_MEM(
        ch_input_fastqs,
        ch_bwamem2_index,
        ch_genome_fasta,
        sort_bam
    )

    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())


    ch_input_bam = BWAMEM2_MEM.out.bam
        .map { meta, bam ->
            def id = "${meta.case_id}_${meta.sample_name}".toString()
            tuple(id, bam)
        }
        .groupTuple()


    ch_bam_meta = BWAMEM2_MEM.out.bam
        .map { meta, bam ->
            def id = "${meta.case_id}_${meta.sample_name}".toString()
            tuple(id, meta.case_id, meta.sample_name, meta.sample_type)
        }
        .groupTuple()
        .map { id, case_id, sample_name, sample_type ->
            def meta = [
                id         : sample_name.first(),
                case_id    : case_id.first(),
                sample_name: sample_name.first(),
                sample_type: sample_type.first()
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
    dedup_bam        = MARKDUPLICATES.out.bam
    dedup_bai        = MARKDUPLICATES.out.bai
    dedup_metrics    = MARKDUPLICATES.out.metrics
    versions         = ch_versions            // channel: versions.yml
}

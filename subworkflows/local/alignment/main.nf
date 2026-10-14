//
// Align reads with BWA-MEM2, merge the per-lane alignments of each sample and remove duplicates
//

include { BWAMEM2_MEM                            } from '../../../modules/nf-core/bwamem2/mem/main'
include { GATK4_MARKDUPLICATES as MARKDUPLICATES } from '../../../modules/nf-core/gatk4/markduplicates/main'


workflow ALIGNMENT {
    take:
    ch_bwamem2_index // channel: [mandatory] [ val(meta), path(index) ]
    ch_genome_fai    // channel: [mandatory] [ val(meta), path(fai) ]
    ch_genome_fasta  // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_input_reads   // channel: [mandatory] [ val(meta), [ path(fastq_1), path(fastq_2) ] ]

    main:

    //
    // MODULE: Run BWA-MEM2 alignment
    //

    ch_reads_for_bwamem2_mem = ch_input_reads
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
            processed_reads
        }

    def sort_bam = true
    BWAMEM2_MEM(
        ch_reads_for_bwamem2_mem,
        ch_bwamem2_index,
        ch_genome_fasta,
        sort_bam
    )

    ch_output_from_bwamem2_mem = BWAMEM2_MEM.out.bam

    // Merge the per-lane BAMs of each sample and remove duplicates
    ch_bwamem2_for_markduplicates = ch_output_from_bwamem2_mem
        .map { meta, bam ->
            def sample_meta = [
                id         : meta.sample_name,
                case_id    : meta.case_id,
                sample_name: meta.sample_name,
                sample_type: meta.sample_type
            ]
            tuple(sample_meta, bam)
        }
        .groupTuple()

    MARKDUPLICATES (
        ch_bwamem2_for_markduplicates,
        ch_genome_fasta.map { _meta, fasta -> fasta }.collect(),
        ch_genome_fai.map { _meta, fai -> fai }.collect()
    )

    emit:
    dedup_bai     = MARKDUPLICATES.out.bai     // channel: [ val(meta), path(bai) ]
    dedup_bam     = MARKDUPLICATES.out.bam     // channel: [ val(meta), path(bam) ]
    dedup_metrics = MARKDUPLICATES.out.metrics // channel: [ val(meta), path(metrics) ]
}


process GRIDSS_PREPROCESS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a8/a899fdab3c3954dcaf82fd3c2f624fcb07051043e4ddf1c748f95966be7fbe96/data':
        'community.wave.seqera.io/library/gridss:2.13.2--8dbd752ea64bc813' }"


    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(genome_fasta)
    tuple val(meta3), path(genome_gridss_index)
    tuple val(meta4), path(genome_fai)
    tuple val(meta5), path(genome_dict)
    path gridss_config

    output:
    tuple val(meta), path("gridss_preprocess/*.gridss.targeted.bam.gridss.working/"),  emit: preprocess_dir
    path  "versions.yml"                                                            ,  emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def arg_config = gridss_config ? "-c ${gridss_config}" : ""

    """
    ln -s \$(find -L ${genome_gridss_index} -regex '.*\\.\\(amb\\|ann\\|pac\\|gridsscache\\|sa\\|bwt\\|img\\|alt\\)') ./

    gridss ${args} \\
        --jvmheap ${Math.round(task.memory.bytes * 0.95)} \\
        --steps preprocess \\
        --reference ${genome_fasta} \\
        --workingdir gridss_preprocess/ \\
        --threads ${task.cpus} ${arg_config} ${bam}


    # Capture versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: \$(CallVariants --version 2>&1 | sed 's/-gridss\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p gridss_preprocess/${prefix}.gridss.targeted.bam.gridss.working/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: "stub"
        R: \$( R --version | sed -n '1p' )
    END_VERSIONS
    """
}

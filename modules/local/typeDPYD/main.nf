process TYPEDPYD {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f0/f013455b8093034814c198fdbd7bee7946b81449e3eb13799599df13dc766a6f/data':
        'community.wave.seqera.io/library/bioconductor-genomicranges_bioconductor-rsamtools_r-data.table_r-optparse_pruned:6a23dd737988db92' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.DPYD.csv")  , emit: csv
    tuple val(meta), path("*.DPYD.json") , emit: json
    tuple val("${task.process}"), val('typeDPYD'), eval("type_DPYD.R --version"), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    type_DPYD.R \\
        -b ${bam} \\
        -c ${prefix}.DPYD.csv \\
        -j ${prefix}.DPYD.json \\
        ${args}

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.DPYD.csv
    touch ${prefix}.DPYD.json
    """
}

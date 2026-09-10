process PREPARE_HETSNPS_FOR_FRANKENPLOT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/be/be0fed5a63cca4fd0f70dcf67ef383c8ea1a305697bda9103e7387ca96e27967/data' :
        'community.wave.seqera.io/library/pysam_tabix_gzip:8124a02a31e03aad' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(bam), path(bai)

    output:
    tuple val(meta), path("*.hetsnps.vcf.gz")    , emit: vcf
    tuple val(meta), path("*.hetsnps.vcf.gz.tbi"), emit: tbi
    tuple val("${task.process}"), val('prepare_hetsnps_for_frankenplot'), eval("prepare_hetsnps_for_frankenplot.py --version"), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    prepare_hetsnps_for_frankenplot.py \\
        --input ${vcf} \\
        --bam ${bam} \\
        --output ${prefix}.hetsnps.vcf.gz \\
        --samplename ${meta.id} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "" | gzip > ${prefix}.hetsnps.vcf.gz
    touch ${prefix}.hetsnps.vcf.gz.tbi
    """
}

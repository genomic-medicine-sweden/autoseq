process PREPARE_HETSNPS_FOR_FRANKENPLOT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.24.0--py312hf5ad864_1' :
        'biocontainers/pysam:0.24.0--py312hf5ad864_1' }"

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

    // The added sample defaults to the meta id. `$args` is appended last, so a `--samplename`
    // set through `ext.args` is the one argparse keeps.
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
    touch ${prefix}.hetsnps.vcf.gz
    touch ${prefix}.hetsnps.vcf.gz.tbi
    """
}

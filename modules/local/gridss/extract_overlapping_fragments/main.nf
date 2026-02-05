
process GRIDSS_EXTRACT_OVERLAPPING_FRAGMENTS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a8/a899fdab3c3954dcaf82fd3c2f624fcb07051043e4ddf1c748f95966be7fbe96/data':
        'community.wave.seqera.io/library/gridss:2.13.2--8dbd752ea64bc813' }"


    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(target_bed)

    output:
    tuple val(meta), path("*.gridss.targeted.bam"),  emit: gridss_targeted_bam
    tuple val(meta), path("*.gridss.targeted.bam.bai"),  emit: gridss_targeted_bai
    path  "versions.yml"                          ,  emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    GRIDSS_JAR=\$(find /opt/conda/share/gridss-* -name "gridss*.jar" | head -n 1 )

    # Extract overlapping fragments using GRIDSS utility
    gridss_extract_overlapping_fragments -w '.' \\
        --targetbed  $target_bed -j \$GRIDSS_JAR \\
        -o ${prefix}.gridss.targeted.bam $bam

    samtools index ${prefix}.gridss.targeted.bam

    # Capture versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: \$(CallVariants --version 2>&1 | sed 's/-gridss\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.gridss.targeted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss_extract_overlapping_fragments: "stub"
        R: \$( R --version | sed -n '1p' )
    END_VERSIONS
    """
}

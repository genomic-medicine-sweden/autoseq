
process GRIPSS_SOMATIC {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/60/60c44edaccb69b42b3c92103dbfedadd474f5496b2b0f2f1e2b47bfa499044a5/data':
        'community.wave.seqera.io/library/hmftools-gripss:2.4--2d77b9970bcdd750' }"


    input:
    tuple val(meta), path(gridss_vcf)
    tuple val(meta2), path(genome_fasta)
    tuple val(meta3), path(genome_fai)
    tuple val(meta4), path(pon_breakends)
    tuple val(meta5), path(pon_breakpoints)
    tuple val(meta6), path(known_fusions)
    tuple val(meta7), path(repeatmasker_annotations)
    tuple val(meta8), path(target_region_bed)
    val genome_version

    output:
    tuple val(meta), path("*.gripss.filtered.somatic.vcf.gz"), path("*.gripss.filtered.somatic.vcf.gz.tbi") ,  emit: filtered_vcf
    tuple val(meta), path("*.gripss.somatic.vcf.gz"), path("*.gripss.somatic.vcf.gz.tbi")                   ,  emit: unfiltered_vcf
    path  "versions.yml"                                                                                    ,  emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.tumor_id}"
    def target_regions_bed_arg = target_region_bed ? "-target_regions_bed ${target_region_bed}" : ''


    """
    gripss \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        ${args} \\
        -sample ${prefix} \\
        -reference ${meta.normal_id} \\
        -vcf ${gridss_vcf} \\
        -ref_genome ${genome_fasta} \\
        -ref_genome_version ${genome_version} \\
        -pon_sgl_file ${pon_breakends} \\
        -pon_sv_file ${pon_breakpoints} \\
        -known_hotspot_file ${known_fusions} \\
        -repeat_mask_file ${repeatmasker_annotations} \\
        ${target_regions_bed_arg} \\
        -output_id somatic \\
        -output_dir ./

    # Capture versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gripss: \$(gripss -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.tumor_id}"

    """
    touch ${prefix}.gripss.filtered.somatic.vcf.gz
    touch ${prefix}.gripss.filtered.somatic.vcf.gz.tbi
    touch ${prefix}.gripss.somatic.vcf.gz
    touch ${prefix}.gripss.somatic.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gripss: "stub"
    END_VERSIONS
    """
}

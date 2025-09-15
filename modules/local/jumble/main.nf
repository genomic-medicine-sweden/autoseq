
process JUMBLE_CNV {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/sarathkumar-murugan/jumble:0.2.1'

    input:
    tuple val(meta), path(bam), path(bai)
    val jumble_ref

    output:
    tuple val(meta), path("${meta.case_id}.jumble.cnv_calls.tsv")

    script:
    """
    jumble-run.R  \\
        --tumor_bam ${tumor_bam} \\
        --normal_bam ${normal_bam} \\
        --reference ${genome_fasta} \\
        --genome_version ${genome_ver} \\
        --output ${meta.case_id}.jumble.cnv_calls.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jumble: \$(jumble --version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """


    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}

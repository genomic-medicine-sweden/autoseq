
process ANNOTATE_CNVS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bd/bda1410103a2570321bb830cc0471125fa7579b0fe026dc09e1e482035c57be5/data' :
        'community.wave.seqera.io/library/pip_logging_pandas:2d7d54d059a6ecf2' }"

    input:
    tuple val(meta), path(cns)
    tuple val(meta2), path(curation_ann)

    output:
    tuple val(meta), path("*_ann.cns"), emit: cns
    path  "versions.yml"          , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    annotate_cnvs.py -i ${cns} -c ${curation_ann} -o ${prefix}_ann.cns

    # Capture versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annotate_cnvs.py: \$( annotate_cnvs.py --version 2>&1 | head -n 1 || echo "unknown" )
        python: \$( python --version | sed -n '1p' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_ann.cns

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annotate_cnvs.py: "stub"
        python: \$( python --version | sed -n '1p' )
    END_VERSIONS
    """
}


process JUMBLE_RUN {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/19/19ea92369dd5c2d099cc7ddc13b7a86c901d88b9b55a556505e25f239fccbbde/data':
        'community.wave.seqera.io/library/bioconductor-aroma.light_bioconductor-bamsignals_bioconductor-bsgenome.hsapiens.ucsc.hg19_bioconductor-bsgenome_pruned:f3bfc207eb7292f6' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(jumbleref)

    output:
    tuple val(meta), path("*.cns"), emit: cns
    tuple val(meta), path("*.cnr"), emit: cnr
    tuple val(meta), path("*_dnacopy.seg"), emit: seg
    tuple val(meta), path("*_profile_bedgraph"), emit: profile_bedgraph
    tuple val(meta), path("*_segments_bedgraph"), emit: segments_bedgraph
    tuple val(meta), path("*.png"), emit: png , optional: true
    tuple val(meta), path("*.RDS"), emit: rds , optional: true
    path  "versions.yml"          , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """

    # Run the tool
    jumble-run.R \\
        $args \\
        -r ${jumbleref} \\
        -b ${bam} \\
        -o "./"

    ## Convert to bedgraph for IGV visualization
    awk -F'\\t' -v OFS='\\t' '\$1 != "chromosome" {print \$1"\\t"\$2"\\t"\$3"\\t"\$6}' \\
         ${prefix}.cnr > ${prefix}_profile_bedgraph

    awk -F'\\t' -v OFS='\\t' '\$1 != "chromosome" {print \$1"\\t"\$2"\\t"\$3"\\t"\$5}' \\
         ${prefix}.cns > ${prefix}_segments_bedgraph

    # Capture versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jumble-run.R: \$( jumble-run.R --version 2>&1 | head -n 1 || echo "unknown" )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.cns
    touch ${prefix}.cnr
    touch ${prefix}_dnacopy.seg
    touch ${prefix}_profile_bedgraph
    touch ${prefix}_segments_bedgraph
    touch ${prefix}.RDS
    touch ${prefix}.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jumble-run.R: "stub"
        R: \$( R --version | sed -n '1p' )
    END_VERSIONS
    """
}

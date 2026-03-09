
process JUMBLE_RUN {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e6/e69d214fba6867a67796351f143d77406d37af48a16b41bc7083618462bd0359/data':
        'community.wave.seqera.io/library/bioconductor-aroma.light_bioconductor-bamsignals_bioconductor-bsgenome.hsapiens.ucsc.hg19_bioconductor-bsgenome_pruned:c58b1261dd6946a5' }"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(jumbleref)

    output:
    tuple val(meta), path("*.cns"), emit: cns
    tuple val(meta), path("*.cnr"), emit: cnr
    tuple val(meta), path("*_dnacopy.seg"), emit: seg
    tuple val(meta), path("*_profile.bedgraph"), emit: profile_bedgraph
    tuple val(meta), path("*_segments.bedgraph"), emit: segments_bedgraph
    tuple val(meta), path("*.png"), emit: png , optional: true
    tuple val(meta), path("*.RDS"), emit: rds , optional: true
    tuple val("${task.process}"), val('Jumble'), eval("jumble-run.R --version 2>&1 | head -n 1 " ),  topic: versions,  emit: versions_jumble

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """

    # Run the tool
    jumble-run.R \\
        $args \\
        -r ${jumbleref} \\
        -b ${bam} \\
        -p ${prefix} \\
        -o "./"

    ## Convert to bedgraph for IGV visualization
    awk -F'\\t' -v OFS='\\t' '\$1 != "chromosome" {print \$1"\\t"\$2"\\t"\$3"\\t"\$6}' \\
         ${prefix}.cnr > ${prefix}_profile.bedgraph

    awk -F'\\t' -v OFS='\\t' '\$1 != "chromosome" {print \$1"\\t"\$2"\\t"\$3"\\t"\$5}' \\
         ${prefix}.cns > ${prefix}_segments.bedgraph

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.cns
    touch ${prefix}.cnr
    touch ${prefix}_dnacopy.seg
    touch ${prefix}_profile.bedgraph
    touch ${prefix}_segments.bedgraph
    touch ${prefix}.RDS
    touch ${prefix}.png

    """
}

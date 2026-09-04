
process JUMBLE_RUN {
    tag "${meta.id}"
    label 'process_medium'

    // Jumble is not distributed via CRAN, Bioconda or Seqera Containers, so the image
    // hosted on the Clinical Genomics Docker Hub account.
    container "docker.io/clinicalgenomics/jumble:0.5.5"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(jumbleref)

    output:
    tuple val(meta), path("*.jumble.csv")        , emit: bins_csv
    tuple val(meta), path("*.cnr")               , emit: cnr
    tuple val(meta), path("*.cns")               , emit: cns
    tuple val(meta), path("*.counts.RDS")        , emit: counts
    tuple val(meta), path("*.genes.csv")         , emit: genes_csv
    tuple val(meta), path("*_profile.bedgraph")  , emit: profile_bedgraph
    tuple val(meta), path("*.qc.csv")            , emit: qc_csv
    tuple val(meta), path("*_dnacopy.seg")       , emit: seg
    tuple val(meta), path("*_segments.bedgraph") , emit: segments_bedgraph
    tuple val(meta), path("*.png")               , emit: png , optional: true
    tuple val("${task.process}"), val('Jumble'), eval("Rscript -e 'cat(as.character(packageVersion(\"Jumble\")))'"),  topic: versions,  emit: versions_jumble

    script:
    def args = task.ext.args ?: ''
    // Jumble names every output after the input BAM, so the prefix is not configurable
    def prefix = bam.baseName

    """

    # Run the tool
    jumble-run.R \\
        $args \\
        -r ${jumbleref} \\
        -b ${bam} \\
        -o "./"

    ## Convert to bedgraph for IGV visualization
    awk -F'\\t' -v OFS='\\t' '\$1 != "chromosome" {print \$1, \$2, \$3, \$7}' \\
         ${prefix}.cnr > ${prefix}_profile.bedgraph

    awk -F'\\t' -v OFS='\\t' '\$1 != "chromosome" {print \$1, \$2, \$3, \$6}' \\
         ${prefix}.cns > ${prefix}_segments.bedgraph

    """

    stub:
    def prefix = bam.baseName

    """
    touch ${prefix}.jumble.csv
    touch ${prefix}.cnr
    touch ${prefix}.cns
    touch ${prefix}.bam.counts.RDS
    touch ${prefix}.genes.csv
    touch ${prefix}_profile.bedgraph
    touch ${prefix}.qc.csv
    touch ${prefix}_dnacopy.seg
    touch ${prefix}_segments.bedgraph
    touch ${prefix}.png

    """
}

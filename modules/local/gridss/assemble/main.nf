
process GRIDSS_ASSEMBLE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a8/a899fdab3c3954dcaf82fd3c2f624fcb07051043e4ddf1c748f95966be7fbe96/data':
        'community.wave.seqera.io/library/gridss:2.13.2--8dbd752ea64bc813' }"


    input:
    tuple val(meta), path(bams), path(bais), path(preprocess_dirs)
    tuple val(meta2), path(genome_fasta)
    tuple val(meta3), path(genome_gridss_index)
    tuple val(meta4), path(genome_fai)
    tuple val(meta5), path(genome_dict)
    tuple val(meta7), path(gridss_config)

    output:
    tuple val(meta), path("gridss_assemble/"),  emit: assemble_dir
    tuple val("${task.process}"), val('gridss'), eval("CallVariants --version 2>&1 | sed 's/-gridss\$//'")  , topic: versions, emit: versions_gridss

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def arg_config = gridss_config ? "-c ${gridss_config}" : ""
    def outdir = "gridss_assemble"

    def bams_list = bams instanceof List ? bams : [bams]

    """
    ln -s \$(find -L ${genome_gridss_index} -regex '.*\\.\\(amb\\|ann\\|pac\\|gridsscache\\|sa\\|bwt\\|img\\|alt\\)') ./

    gridss ${args} \\
        --jvmheap ${Math.round(task.memory.bytes * 0.95)} \\
        --steps assemble \\
        --reference ${genome_fasta} \\
        --workingdir ${outdir}/work/ \\
        --assembly ${outdir}/${prefix}.sv.assembly.bam \\
        --threads ${task.cpus} ${arg_config} ${bams_list.join(' ')}


    # Capture versions
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: \$(CallVariants --version 2>&1 | sed 's/-gridss\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p gridss_assemble/work/
    touch gridss_assemble/${prefix}.sv.assembly.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gridss: "stub"
        R: \$( R --version | sed -n '1p' )
    END_VERSIONS
    """
}

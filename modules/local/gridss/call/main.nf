
process GRIDSS_CALL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a8/a899fdab3c3954dcaf82fd3c2f624fcb07051043e4ddf1c748f95966be7fbe96/data':
        'community.wave.seqera.io/library/gridss:2.13.2--8dbd752ea64bc813' }"


    input:
    tuple val(meta), path(bams), path(bais), path(assemble_dir)
    tuple val(meta2), path(genome_fasta)
    tuple val(meta3), path(genome_gridss_index)
    tuple val(meta4), path(genome_fai)
    tuple val(meta5), path(genome_dict)
    tuple val(meta6), path(blacklist)
    tuple val(meta7), path(gridss_config)

    output:
    tuple val(meta), path("gridss_call/*.sv.gridss.vcf.gz"),  emit: vcf
    tuple val("${task.process}"), val('gridss'), eval("CallVariants --version 2>&1 | sed 's/-gridss\$//'")  , topic: versions, emit: versions_gridss

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def arg_config = gridss_config ? "-c ${gridss_config}" : ""
    def outdir = "gridss_call"

    def bams_list = bams instanceof List ? bams : [bams]

    """
    ln -s \$(find -L ${genome_gridss_index} -regex '.*\\.\\(amb\\|ann\\|pac\\|gridsscache\\|sa\\|bwt\\|img\\|alt\\)') ./

    gridss ${args} \\
        --jvmheap ${Math.round(task.memory.bytes * 0.95)} \\
        --steps call \\
        --reference ${genome_fasta} \\
        --blacklist ${blacklist} \\
        --workingdir ${outdir}/work/ \\
        --assembly ${outdir}/${prefix}.sv.assemblies.bam \\
        --output ${outdir}/${prefix}.sv.gridss.vcf.gz \\
        --threads ${task.cpus} ${arg_config} ${bams_list.join(' ')}


    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p gridss_call/work/
    touch gridss_call/${prefix}.sv.assemblies.bam
    touch gridss_call/${prefix}.sv.gridss.vcf.gz

    """
}

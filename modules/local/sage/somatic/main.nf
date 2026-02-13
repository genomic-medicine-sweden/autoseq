// This module is adapted from nf-core/oncoanalyser for use in this workflow
// NOTE(SW): logic that determines BQR outputs assumes '-output_vcf' is a path that includes at least one non-empty directory (e.g. /path/to/results/filename.vcf)

process SAGE_SOMATIC {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hmftools-sage:4.1--hdfd78af_0' :
        'biocontainers/hmftools-sage:4.1--hdfd78af_0' }"

    input:
    tuple val(meta), path(input), path(input_index), path(intervals)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)
    val genome_ver
    path sage_pon
    path sage_known_hotspots_somatic
    path sage_highconf_regions
    path ensembl_data_resources
    val targeted_mode
    val min_tumor_vaf
    val min_tumor_qual
    val min_map_quality
    val hotspot_tumor_qual
    val min_avg_base_qual

    output:
    tuple val(meta), path('somatic/*.sage.somatic.vcf.gz')     , emit: vcf
    tuple val(meta), path('somatic/*.sage.somatic.vcf.gz.tbi') , emit: tbi
    tuple val(meta), path('somatic/')                          , emit: sage_dir
    path 'versions.yml'                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def log_level_arg = task.ext.log_level ? "-log_level ${task.ext.log_level}" : ''
    def reference_arg = meta.normal_id != null ? "-reference ${meta.normal_id}" : ''

    // BAMs
    def tumor_bam = input[0]
    def reference_bam = input[1]
    def reference_bam_arg = reference_bam ? "-reference_bam ${reference_bam}" : ''
    def high_depth_mode_arg = targeted_mode ? '-high_depth_mode': ''
    def panel_bed = intervals ? "-panel_bed ${intervals} -panel_only " : ''

    def hotspot_args = hotspot_tumor_qual ? "-hotspot_min_tumor_qual ${hotspot_tumor_qual}" : ''
    def taf_args = min_tumor_vaf ? " -hard_min_tumor_vaf ${min_tumor_vaf} -hotspot_min_tumor_vaf ${min_tumor_vaf} -panel_min_tumor_vaf ${min_tumor_vaf}" : ''
    def panel_tumor_qual_args = min_tumor_qual ? " -panel_min_tumor_qual ${min_tumor_qual} " : ''
    def map_qual_args = min_map_quality ? " -min_map_quality ${min_map_quality} " : ''
    def base_qual_args = min_avg_base_qual ? " -min_avg_base_qual ${min_avg_base_qual}" : ''

    """
    mkdir -p somatic/

    sage \\
        -Xmx${Math.round(task.memory.bytes * 0.95)} \\
        ${args} \\
        ${reference_arg} \\
        ${reference_bam_arg} \\
        ${high_depth_mode_arg} \\
        ${panel_bed} \\
        -tumor ${meta.tumor_id} \\
        -tumor_bam ${tumor_bam} \\
        -ref_genome ${fasta} \\
        -ref_genome_version ${genome_ver} \\
        -hotspots ${sage_known_hotspots_somatic} \\
        -high_confidence_bed ${sage_highconf_regions} \\
        -ensembl_data_dir ${ensembl_data_resources} \\
        -bqr_write_plot \\
        -skip_msi_jitter \\
        -threads ${task.cpus} \\
        ${hotspot_args} \\
        ${taf_args}     \\
        ${panel_tumor_qual_args}  \\
        ${map_qual_args}  \\
        ${base_qual_args}  \\
        ${log_level_arg} \\
        -output_vcf somatic/${meta.tumor_id}.sage.somatic.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sage: \$(sage -version | sed 's/^.* //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p somatic/

    touch somatic/${meta.tumor_id}.sage.somatic.vcf.gz
    touch somatic/${meta.tumor_id}.sage.somatic.vcf.gz.tbi
    touch somatic/${meta.tumor_id}.gene.coverage.tsv
    touch somatic/${meta.tumor_id}.sage.bqr.png
    touch somatic/${meta.tumor_id}.sage.bqr.tsv

    ${ (meta.normal_id != null) ? "touch somatic/${meta.normal_id}.sage.bqr.{png,tsv}" : '' }

    echo -e '${task.process}:\\n  stub: noversions\\n' > versions.yml
    """
}

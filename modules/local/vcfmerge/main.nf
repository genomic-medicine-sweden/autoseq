
process SOMATIC_VCFMERGE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/60/60c6a151bdbba0a0a37aac6110ddd8c0837074a427a2e9512d2ad30fb3ebc44e/data' :
        'community.wave.seqera.io/library/bcftools_tabix:04336756d7b46b1b' }"

    input:
    tuple val(meta), path(mutect_vcf), path(sage_vcf)


    output:
    tuple val(meta), path("*merged.vcf"), emit: merged_vcf
    path  "versions.yml"          , emit: versions

    script:
    def args = task.ext.args ?: ''
    def targets_bed = targets ? "-t ${targets} " : ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    /*
    This module provides a series of bcftools commands to combine VCF files generated
    by two different somatic variant callers. The process consists of the following steps:

        1. Reorder the sample names in SAGE VCF file to ensure consistency. This step is
        necessary because bcftools requires the sample names to be in the same order across
        all input files for successful merging.

        2. Merge the VCF files using bcftools concat. This command concatenates the VCF
        files, producing a single output file that contains the combined variant calls from
        both sources.

    These steps facilitate the integration of variant data from multiple callers, enabling
    downstream analyses that require a unified VCF file.


    TODO: The current implementation has certain limitation. Need better solution.
    */

    """
    # Re-order the sample names in SAGE vcf file
    bcftools query -l ${sage_vcf} | sort > sample_names.txt
    bcftools view -0z -S sample_names.txt ${sage_vcf} -o ${prefix}_ordered.vcf.gz
    tabix -p vcf ${prefix}_ordered.vcf.gz

    # VCF concatenation
    bcftools concat -a -D ${mutect_vcf} ${prefix}_ordered.vcf.gz \\
        | bgzip > ${prefix}-all.somatic.vcf.gz
    tabix -p vcf ${prefix}-all.somatic.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version 2>&1 | head -n 1 || echo "unknown" )
        tabix:  \$( tabix --version 2>&1 | head -n 1 || echo "unknown" )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}-all.somatic.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version 2>&1 | head -n 1 || echo "unknown" )
        tabix:  \$( tabix --version 2>&1 | head -n 1 || echo "unknown" )
    END_VERSIONS
    """
}

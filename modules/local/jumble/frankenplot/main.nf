process JUMBLE_FRANKENPLOT {
    tag "${meta.id}"
    label 'process_single'

    // Jumble is not distributed via CRAN, Bioconda or Seqera Containers, so use the image
    // hosted on the Clinical Genomics Docker Hub account.
    container "docker.io/clinicalgenomics/jumble:0.5.6"

    input:
    tuple val(meta), path(tcns), path(ncns), path(t_csv), path(n_csv), path(svcf), path(gvcf), path(het_snps), path(dpyd_json), path(dpyd_csv)

    output:
    tuple val(meta), path("*.frankenplot.html"), emit: html
    tuple val("${task.process}"), val('Jumble'), eval("Rscript -e 'cat(as.character(packageVersion(\"Jumble\")))'"),  topic: versions,  emit: versions_jumble

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def normal_cns_arg = ncns       ? "--normal-cns ${ncns}"           : ''
    def normal_csv_arg = n_csv      ? "--normal-jumble-csv ${n_csv}"   : ''
    def het_snps_arg   = het_snps   ? "--tumor-snp-vcf ${het_snps}"    : ''
    def somatic_arg    = svcf       ? "--somatic-vcf ${svcf}"          : ''
    def germline_arg   = gvcf       ? "--germline-vcf ${gvcf}"         : ''
    def dpyd_json_arg  = dpyd_json  ? "--dpyd-json ${dpyd_json}"       : ''
    def dpyd_csv_arg   = dpyd_csv   ? "--dpyd-csv ${dpyd_csv}"         : ''

    """
    jumble-frankenplot.R \\
        $args \\
        -c ${tcns} \\
        -j ${t_csv} \\
        ${normal_cns_arg} \\
        ${normal_csv_arg} \\
        ${het_snps_arg} \\
        ${somatic_arg} \\
        ${germline_arg} \\
        ${dpyd_json_arg} \\
        ${dpyd_csv_arg} \\
        -o ${prefix}.frankenplot.html

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.frankenplot.html
    """
}

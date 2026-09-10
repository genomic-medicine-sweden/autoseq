process JUMBLE_FRANKENPLOT {
    tag "${meta.id}"
    label 'process_single'

    // Jumble is not distributed via CRAN, Bioconda or Seqera Containers, so use the image
    // hosted on the Clinical Genomics Docker Hub account.
    container "docker.io/clinicalgenomics/jumble:0.5.6"

    input:
    tuple val(meta), path(cns), path(jumble_csv), path(vcf), path(het_snps), path(dpyd)

    output:
    tuple val(meta), path("*.frankenplot.html"), emit: html
    tuple val("${task.process}"), val('Jumble'), eval("Rscript -e 'cat(as.character(packageVersion(\"Jumble\")))'"),  topic: versions,  emit: versions_jumble

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // `cns`, `jumble_csv`, `vcf` and `dpyd` are lists holding the tumour entry first and the
    // optional normal / secondary entry second. Nextflow unwraps a single-element list into a
    // bare path, and `size()` on a path returns its byte count rather than an element count,
    // so normalise back to a list before unpacking positionally.
    def cns_files  = cns        instanceof List ? cns        : [cns]
    def csv_files  = jumble_csv instanceof List ? jumble_csv : [jumble_csv]
    def vcf_files  = vcf        instanceof List ? vcf        : [vcf]
    def dpyd_files = dpyd       instanceof List ? dpyd       : [dpyd]

    def normal_cns_arg = cns_files.size()  > 1 ? "--normal-cns ${cns_files[1]}"             : ''
    def normal_csv_arg = csv_files.size()  > 1 ? "--normal-jumble-csv ${csv_files[1]}"      : ''
    def het_snps_arg   = het_snps              ? "--tumor-snp-vcf ${het_snps}"              : ''
    def somatic_arg    = vcf_files.size()  > 0 ? "--somatic-vcf ${vcf_files[0]}"            : ''
    def germline_arg   = vcf_files.size()  > 1 ? "--germline-vcf ${vcf_files[1]}"           : ''
    def dpyd_json_arg  = dpyd_files.size() > 0 ? "--dpyd-json ${dpyd_files[0]}"             : ''
    def dpyd_csv_arg   = dpyd_files.size() > 1 ? "--dpyd-csv ${dpyd_files[1]}"              : ''

    """
    jumble-frankenplot.R \\
        $args \\
        -c ${cns_files[0]} \\
        -j ${csv_files[0]} \\
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

# nf-core/autoseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [date]

Initial release of nf-core/autoseq, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- [ #1 ](https://github.com/imsarath/nf-autoseq/pull/1) Added subworkflow for BAM QC, Reference genome configuration and nf-test integration
- [ #3 ](https://github.com/imsarath/nf-autoseq/pull/3) Added subworkflow for somatic variant callers - GATK4 MuTect2, Hmftools - SAGE somatic
- [ #7 ](https://github.com/imsarath/nf-autoseq/pull/7) Added subworkflow for cnv calling (Jumble) and cancer-specific information annotations.
- [ #13 ](https://github.com/imsarath/nf-autoseq/pull/13) Added subworkflow for structural variant calling - GRIDSS
- [ #16 ](https://github.com/imsarath/nf-autoseq/pull/16) Added subworkflow for germline variant calling - GATK4 haplotypecaller and minor bug fixes
- [ #19 ](https://github.com/imsarath/nf-autoseq/pull/19) Repo setup for Genome Medicine Sweden Org
- [ #21 ](https://github.com/genomic-medicine-sweden/nf-autoseq/pull/21) Added sub-workflow for tumor biomarker profiling - purecn-run and typeDPYD

### `Fixed`

- [ #22 ](https://github.com/genomic-medicine-sweden/nf-autoseq/pull/22) Added PR `write` permissions for `add_pr_checklist_comment` action.
- [ #21 ](https://github.com/genomic-medicine-sweden/nf-autoseq/pull/21) Enabled synchronize in nf-test.yml so that GitHub CI/CD stays active for subsequent updates following a review request.
- [ #21 ](https://github.com/genomic-medicine-sweden/nf-autoseq/pull/21) Modified jumble-run.R to optimize the execution speed of version printing.
- [ #21 ](https://github.com/genomic-medicine-sweden/nf-autoseq/pull/21) Updated GATK Mutect2 parameters to include germline variants in the unfiltered VCF, enabling compatibility with PureCN downstream analysis.
- [ #21 ](https://github.com/genomic-medicine-sweden/nf-autoseq/pull/21) Updated the minimum Nextflow version to fix nf-test failures in the GitHub Actions workflow.
- [ #21 ](https://github.com/genomic-medicine-sweden/nf-autoseq/pull/21) Added the modified nf_core_autoseq logo images to .nf-core.yml to ignore them and resolve nf-core linting failures.

### `Dependencies`

### `Deprecated`

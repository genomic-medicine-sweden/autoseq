# genomic-medicine-sweden/autoseq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [unreleased<!-- TODO nf-core: replace with date on release -->]

Initial release of genomic-medicine-sweden/autoseq, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- [ #1 ](https://github.com/imsarath/nf-autoseq/pull/1) Added subworkflow for BAM QC, Reference genome configuration and nf-test integration
- [ #3 ](https://github.com/imsarath/nf-autoseq/pull/3) Added subworkflow for somatic variant callers - GATK4 MuTect2, Hmftools - SAGE somatic
- [ #7 ](https://github.com/imsarath/nf-autoseq/pull/7) Added subworkflow for cnv calling (Jumble) and cancer-specific information annotations.
- [ #13 ](https://github.com/imsarath/nf-autoseq/pull/13) Added subworkflow for structural variant calling - GRIDSS
- [ #16 ](https://github.com/imsarath/nf-autoseq/pull/16) Added subworkflow for germline variant calling - GATK4 haplotypecaller and minor bug fixes
- [ #19 ](https://github.com/imsarath/nf-autoseq/pull/19) Repo setup for Genome Medicine Sweden Org
- [ #21 ](https://github.com/genomic-medicine-sweden/nf-autoseq/pull/21) Added sub-workflow for tumor biomarker profiling - purecn-run and typeDPYD
- [ #40 ](https://github.com/genomic-medicine-sweden/autoseq/pull/40) Added nf-test for the `ANNOTATE_CNVS` local module covering somatic, germline, and stub cases.
- [ #42 ](https://github.com/genomic-medicine-sweden/autoseq/pull/42) Added nf-test for `gridss/extract_overlapping_fragments` module.
- [ #43 ](https://github.com/genomic-medicine-sweden/autoseq/pull/43) Added nf-test for the `gridss/preprocess` local module covering a targeted BAM scenario and stub case.
- [ #50 ](https://github.com/genomic-medicine-sweden/autoseq/pull/50) Added `PREPARE_REFERENCES` subworkflow to build the BWA-MEM2 index and untar the VEP cache, GRIDSS index and HMF ensembl_data references, with nf-test coverage.
- [ #55 ](https://github.com/genomic-medicine-sweden/autoseq/pull/55) Enabled an end-to-end `-profile test` run of the paired tumor/normal workflow on real GRCh37 test data.
- [ #71 ](https://github.com/genomic-medicine-sweden/autoseq/pull/71) Added a `test_umi` profile and pipeline-level nf-test covering the UMI alignment branch on paired tumor/normal test data.
- [ #85 ](https://github.com/genomic-medicine-sweden/autoseq/pull/85) Added nf-test and `meta.yml` for the `ALIGNMENT` subworkflow covering a multi-lane tumor sample and a stub case.
- [ #86 ](https://github.com/genomic-medicine-sweden/autoseq/pull/86) Added nf-test and `meta.yml` for the `QC_ALIGNMENT` subworkflow covering a tumor BAM with the panel interval list and a stub case.
- [ #87 ](https://github.com/genomic-medicine-sweden/autoseq/pull/87) Added nf-test and `meta.yml` for the `CALL_CNVS` subworkflow covering paired tumor/normal samples and a stub case.
- [ #88 ](https://github.com/genomic-medicine-sweden/autoseq/pull/88) Added nf-test and `meta.yml` for the `CALL_GERMLINE_SNVS` subworkflow covering a normal sample with the panel interval list and a stub case.
- [ #89 ](https://github.com/genomic-medicine-sweden/autoseq/pull/89) Added nf-test and `meta.yml` for the `CALL_SVS` subworkflow covering a paired tumor/normal case and a stub case.
- [ #90 ](https://github.com/genomic-medicine-sweden/autoseq/pull/90) Added nf-test and `meta.yml` for the `CALL_SOMATIC_SNVS` subworkflow covering a paired tumor/normal case and a stub case.
- [ #91 ](https://github.com/genomic-medicine-sweden/autoseq/pull/91) Added the `channelFromPathWithMeta` helper to build `[[id:...], file]` reference channels from a path, or `channel.empty()` when the path is null.
- [ #96 ](https://github.com/genomic-medicine-sweden/autoseq/pull/96) Added the `ANNOTATE_GERMLINE_TAF` subworkflow to annotate germline variants with the tumor allele fraction, with nf-test and `meta.yml`.

### `Changed`

- [ #35 ](https://github.com/genomic-medicine-sweden/nf-autoseq/pull/35) Refactor repository name from `nf-core/autoseq` to `genomic-medicine-sweden/autoseq`.
- [ #35 ](https://github.com/genomic-medicine-sweden/nf-autoseq/pull/35) Expanded `CONTRIBUTING.md` with detailed PR title conventions, review guidelines, and GMS-specific contribution workflow.
- [ #35 ](https://github.com/genomic-medicine-sweden/nf-autoseq/pull/35) Updated bug report and issue templates to remove nf-core Slack references and align with GMS organisation.
- [ #35 ](https://github.com/genomic-medicine-sweden/nf-autoseq/pull/35) Updated `README.md` contributors section.
- [ #36 ](https://github.com/genomic-medicine-sweden/autoseq/pull/36) Template update for nf-core/tools v4.0.2.
- [ #36 ](https://github.com/genomic-medicine-sweden/autoseq/pull/36) Updated multiqc module from 1.32 to 1.34 and fastqc module (included in the template update).
- [ #36 ](https://github.com/genomic-medicine-sweden/autoseq/pull/36) Updated the minimum required nextflow version to 25.10.4 (included in new template).
- [ #57 ](https://github.com/genomic-medicine-sweden/autoseq/pull/57) Regenerated the `PREPARE_REFERENCES` nf-test snapshot to match the updated HMF `ensembl_data` reference CSVs.
- [ #63 ](https://github.com/genomic-medicine-sweden/autoseq/pull/63) Updated the nf-core `fastq_create_umi_consensus_fgbio` subworkflow and its fgbio and samtools modules.
- [ #91 ](https://github.com/genomic-medicine-sweden/autoseq/pull/91) Reference channels in `main.nf` now use `channelFromPathWithMeta`.
- [ #95 ](https://github.com/genomic-medicine-sweden/autoseq/pull/95) Updated the contribution guidelines `docs/CONTRIBUTING.md` with new publishing strategy for pipeline outputs.

### `Fixed`

- [ #35 ](https://github.com/genomic-medicine-sweden/nf-autoseq/pull/35) Rebuilt `nextflow_schema.json` to pass `nf-core schema build` validation after the organisation rename.
- [ #22 ](https://github.com/genomic-medicine-sweden/nf-autoseq/pull/22) Added PR `write` permissions for `add_pr_checklist_comment` action.
- [ #21 ](https://github.com/genomic-medicine-sweden/nf-autoseq/pull/21) Enabled synchronize in nf-test.yml so that GitHub CI/CD stays active for subsequent updates following a review request.
- [ #21 ](https://github.com/genomic-medicine-sweden/nf-autoseq/pull/21) Modified jumble-run.R to optimize the execution speed of version printing.
- [ #21 ](https://github.com/genomic-medicine-sweden/nf-autoseq/pull/21) Updated GATK Mutect2 parameters to include germline variants in the unfiltered VCF, enabling compatibility with PureCN downstream analysis.
- [ #21 ](https://github.com/genomic-medicine-sweden/nf-autoseq/pull/21) Updated the minimum Nextflow version to fix nf-test failures in the GitHub Actions workflow.
- [ #21 ](https://github.com/genomic-medicine-sweden/nf-autoseq/pull/21) Added the modified nf_core_autoseq logo images to .nf-core.yml to ignore them and resolve nf-core linting failures.
- [ #37 ](https://github.com/genomic-medicine-sweden/autoseq/pull/38) Updated documentation reference in `multiqc_config.yml` and disabled the nf-core linting check for multiqc config.
- [ #41 ](https://github.com/genomic-medicine-sweden/autoseq/pull/41) Refactored `annotate_cnvs` local module to accept `sample_type` as part of the input tuple instead of deriving it from `meta.sample_type` inside the module.
- [ #52 ](https://github.com/genomic-medicine-sweden/autoseq/pull/52) Replaced leftover `oncorefiner` references in `docs/CONTRIBUTING.md` with the correct pipeline name `autoseq`.
- [ #53 ](https://github.com/genomic-medicine-sweden/autoseq/pull/53) Added `-tumor-segmentation` argument to `GATK4_CALCULATECONTAMINATION` so the tumor segmentation table is written for downstream filtering.
- [ #54 ](https://github.com/genomic-medicine-sweden/autoseq/pull/54) Emit a tabix index for the Mutect2 pass-filtered VCF and pass it to `SOMATIC_VCFMERGE` so `bcftools concat -a` can load the index.
- [ #56 ](https://github.com/genomic-medicine-sweden/autoseq/pull/56) Use `meta.id` instead of `meta.tumor_id` for the `PURECN_RUN` output prefix.
- [ #58 ](https://github.com/genomic-medicine-sweden/autoseq/pull/58) Replaced the stubbed `PURECN_RUN` and `PROFILE_TUMOR_BIOMARKERS` nf-tests with real-VCF runs, added dedicated stub cases, and regenerated the snapshots.
- [ #59 ](https://github.com/genomic-medicine-sweden/autoseq/pull/59) Added `--normal-sample ${meta.normal_id}` to `GATK4_MUTECT2` so the normal sample is correctly identified in paired tumor/normal calling.
- [ #62 ](https://github.com/genomic-medicine-sweden/autoseq/pull/62) Excluded non-deterministic Jumble CNV outputs from the `tests/default.nf.test` content snapshot to fix md5 failures in CI.
- [ #72 ](https://github.com/genomic-medicine-sweden/autoseq/pull/72) Added `--sample` and `--library` to `FASTQTOBAM` so the read group sample name matches `meta.id` for downstream callers.
- [ #73 ](https://github.com/genomic-medicine-sweden/autoseq/pull/73) Moved the `ZIPPERBAMS_(PRE|POST)` tag options from `ext.args` to `ext.args2` so they are passed to `ZipperBams` instead of to the fgbio wrapper.
- [ #74 ](https://github.com/genomic-medicine-sweden/autoseq/pull/74) Collected the reference channel passed to `UMI_PROCESSING` so it is reusable across all samples instead of being consumed by the first one.
- [ #94 ](https://github.com/genomic-medicine-sweden/autoseq/pull/94) Sorted grouped lane FASTQs by filename before `CAT_FASTQ` in the UMI branch, since `groupTuple` does not guarantee ordering and multi-lane samples could be concatenated in a non-reproducible order.

### `Dependencies`

### `Deprecated`

# nf-core/autoseq: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/autoseq/usage](https://nf-co.re/autoseq/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

**nf-autoseq** is an automated Nextflow-based bioinformatics pipeline developed for the comprehensive analysis of cancer genomics data, specifically optimized for deep targeted sequencing and whole-exome sequencing (WES). The pipeline facilitates a reproducible, end-to-end workflow transitioning from raw sequence reads to refined, annotated variant calls suitable for clinical and research interpretation.

**Note**: Currently, it is in active development and not yet ready for production use.

## Prerequisites

To ensure portability and ease of deployment, nf-autoseq requires the following software to be installed on your host system:

- Nextflow
- Container/Package Manager: The pipeline supports multiple environments to manage software dependencies:
  `Docker`, `Singularity`, `Apptainer`, `Conda`, etc. (see [profiles](#profile) for more details).

## Running the pipeline

To run the nf-autoseq pipeline, you just need to get two things ready:

- Samplesheet
- Reference Datasets

It is important to format these files correctly so the pipeline can run automatically without any errors. Please follow the instructions below to download your references and set up your samplesheet.

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-autoseq/main.nf \
        --input ./samplesheet.csv \
        --outdir ./results \
        --genome GRCh37  \
        -profile docker \
        --ref_genomes_base references \
        --panel probio_comprehensive4
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

#### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

````csv title="samplesheet.csv"
```csv
case_id,sample_name,sample_type,lane,fastq_1,fastq_2,bam
PATIENT_ID,TUMOR_ID,tumor,L2,/path/to/SAMPLE_L2_R1_001.fastq.gz,/path/to/SAMPLE_L2_R2_001.fastq.gz,
PATIENT_ID,TUMOR_ID,tumor,L3,/path/to/SAMPLE_L3_R1_001.fastq.gz,/path/to/SAMPLE_L3_R2_001.fastq.gz,
PATIENT_ID,TUMOR_ID,tumor,L4,/path/to/SAMPLE_L4_R1_001.fastq.gz,/path/to/SAMPLE_L4_R2_001.fastq.gz,
PATIENT_ID,TUMOR_ID,tumor,L5,/path/to/SAMPLE_L5_R1_001.fastq.gz,/path/to/SAMPLE_L5_R2_001.fastq.gz,
````

> [!IMPORTANT]
> Note: current version of the pipeline does not support BAM files as input, so the `bam` column should be left empty. Support for BAM input is planned for a future release.

#### Full samplesheet

The samplesheet can also contain multiple samples, for example:

```csv title="samplesheet.csv"
case_id,sample_name,sample_type,lane,fastq_1,fastq_2,bam
PATIENT_ID,TUMOR_ID,tumor,L2,/path/to/SAMPLE_L2_R1_001.fastq.gz,/path/to/SAMPLE_L2_R2_001.fastq.gz,
PATIENT_ID,TUMOR_ID,tumor,L3,/path/to/SAMPLE_L3_R1_001.fastq.gz,/path/to/SAMPLE_L3_R2_001.fastq.gz,
PATIENT_ID,TUMOR_ID,tumor,L4,/path/to/SAMPLE_L4_R1_001.fastq.gz,/path/to/SAMPLE_L4_R2_001.fastq.gz,
PATIENT_ID,TUMOR_ID,tumor,L5,/path/to/SAMPLE_L5_R1_001.fastq.gz,/path/to/SAMPLE_L5_R2_001.fastq.gz,
PATIENT_ID,NORMAL_ID,normal,L4,/path/to/SAMPLE_L4_R1_001.fastq.gz,/path/to/SAMPLE_L4_R2_001.fastq.gz,
PATIENT_ID,NORMAL_ID,normal,L5,/path/to/SAMPLE_L5_R1_001.fastq.gz,/path/to/SAMPLE_L5_R2_001.fastq.gz,
PATIENT_ID,NORMAL_ID,normal,L6,/path/to/SAMPLE_L6_R1_001.fastq.gz,/path/to/SAMPLE_L6_R2_001.fastq.gz,
```

| Column        | Description                                                                                                                                                                            |
| ------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `case_id`     | Unique identifier for the patient or case.                                                                                                                                             |
| `sample_name` | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `sample_type` | Type of sample (e.g., tumor, normal).                                                                                                                                                  |
| `lane`        | Lane number of the sequencing run.                                                                                                                                                     |
| `fastq_1`     | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `fastq_2`     | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `bam`         | Full path to BAM file. This column is optional and should be left empty if not used.                                                                                                   |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

### Reference Files and Panels

The `nf-autoseq` pipeline relies on reference genomes and gene panels to perform its analyses. By default, the pipeline will download these files from public repositories (eg. Ensembl) as needed during the run. However, you can also specify a local directory containing these files if you have them available.

#### Reference Directory Structure

```
references/
└── GRCh37/
    ├── genome/             # FASTA, FAI, and DICT files
    ├── bwamem2_index/      # BWA-mem2 alignment indices
    ├── annotations/        # Germline resources and curation CSVs
    ├── hmfdata/            # SAGE blocklists, hotspots, and driver panels
    ├── vep/                # Ensembl VEP cache
    ├── gridss_index        # GRIDSS cache and bwa-index
    └── gridss/             # GRIDSS configuration, PONs, and fusions
└── GRCh38/
    ├── genome/             # FASTA, FAI, and DICT files
    ...
```

> [!IMPORTANT]
> Automation Note: Automated downloading and generation of these reference files is not yet implemented. Users must ensure the reference path is correctly populated before running the pipeline. Full automation of reference setup is planned for a future release.

### Panel Support

The pipeline uses the `--panel` parameter to load specific genomic coordinates and bait information required for targeted sequencing analysis.

#### Currently Supported Panels

The following panels are pre-configured in the pipeline:

| Panel ID              | Description             | Components Included                                     |
| --------------------- | ----------------------- | ------------------------------------------------------- |
| probio_comprehensive3 | ProBio Comprehensive v3 | BED (slopped), Interval lists, and Jumble RDS reference |
| gmck_v3               | GMCK v3                 | BED (slopped), Interval lists, and Jumble RDS reference |

#### Custom Panels

If you have a custom panel that is not included in the above list, you can specify it using the `--panel_bed` parameter. The pipeline expects a BED file with the genomic coordinates of the panel targets. Pipeline will automatically generate the necessary interval lists. This option is not implemented in the current version of the pipeline.

#### Parameter Usage

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-autoseq/main.nf -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/autoseq releases page](https://github.com/nf-core/autoseq/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```

#!/usr/bin/env Rscript

# Markus Mayrhofer 2022-2025
# Dependencies ------------------------------------------------------------
{
    suppressPackageStartupMessages(library(optparse))
    suppressPackageStartupMessages(library(stringr))
    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(GenomicRanges))
    suppressPackageStartupMessages(library(Rsamtools))
    suppressPackageStartupMessages(library(RJSONIO))
}


# Version -----------------------------------------------------------------

VERSION <- "0.1.0"


# Options -----------------------------------------------------------------

option_list <- list(
    make_option(c("-b", "--bam_file"), action = "store", type = "character", default = NULL,
                help = "A .bam file aligned to hg19 or hg38. A matching bam index (.bai) file is presumed."),
    make_option(c("-c", "--csv_file"), action = "store", type = "character", default = NULL,
                help = "Output csv file."),
    make_option(c("-j", "--json_file"), action = "store", type = "character", default = NULL,
                help = "Output json file."),
    make_option(c("-v", "--version"), action = "store_true", default = FALSE,
                help = "Print version and exit.")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (isTRUE(opt$version)) {
    cat(VERSION, "\n")
    quit(status = 0)
}



# save.image('ws.Rdata')
# stop()
# load('ws.Rdata')



# Variant definition ------------------------------------------------------

{
    # Define a single data.table containing all variant information, including
    # separate columns for hg19 and hg38 positions for easy switching.
    variants_table <- data.table(
        alias = rep(c("*2A", "*13", "c.2846A>T", "HapB3"), each = 2),
        rsid = rep(c("rs3918290", "rs55886062", "rs67376798", "rs56038477"), each = 2),
        chromosome = "1",
        pos_hg19 = rep(c(97915614, 97981343, 97547947, 98039419), each = 2),
        pos_hg38 = rep(c(97341350, 97407079, 96973683, 97465155), each = 2),
        nucleotide = c("C", "T", "A", "C", "T", "A", "C", "T"),
        allele = rep(c("wt", "alt"), 4),
        activity = rep(c(1, 0, 1, 0, 1, 0.5, 1, 0.5), each = 1)
    )

    variants_table[, activity_loss := 1 - activity]
}

# Genotype text -------------------------------------------------------------

genotype_text <- list(
    '*1/*1'=
        'No evidence of the listed SNP alleles associated with deficient DPYD activity was found.
        Normal enzymatic activity can be expected.
        In case of treatment with kapecitabin or fluorouracil, normal start dose is recommended.',
    '*1/HapB3'=
        'This patient is heterozygous for an allele associated with deficient DPYD activity.
        Slightly deficient enzymatic activity can be expected.
        In case of treatment with kapecitabin or fluorouracil, 50-75% of normal start dose is recommended.',
    '*1/c.2846A>T'=
        'This patient is heterozygous for an allele associated with deficient DPYD activity.
        Slightly deficient enzymatic activity can be expected.
        In case of treatment with kapecitabin or fluorouracil, 50-75% of normal start dose is recommended.',
    '*1/*2A'=
        'This patient is heterozygous for an allele associated with loss of DPYD activity.
        Partially deficient enzymatic activity can be expected.
        In case of treatment with kapecitabin or fluorouracil, 50% of normal start dose is recommended.',
    '*1/*13'=
        'This patient is heterozygous for an allele associated with loss of DPYD activity.
        Partially deficient enzymatic activity can be expected.
        In case of treatment with kapecitabin or fluorouracil, 50% of normal start dose is recommended.',
    'HapB3/HapB3'=
        'This patient is homozygous for an allele associated with deficient DPYD activity.
        Partially deficient enzymatic activity can be expected.
        In case of treatment with kapecitabin or fluorouracil, 50% of normal start dose is recommended.',
    'c.2846A>T/c.2846A>T'=
        'This patient is homozygous for an allele associated with deficient DPYD activity.
        Partially deficient enzymatic activity can be expected.
        In case of treatment with kapecitabin or fluorouracil, 50% of normal start dose is recommended.',
    'c.2846A>T/HapB3'=
        'This patient carries two different alleles associated with deficient DPYD activity.
        Partially deficient enzymatic activity can be expected.
        In case of treatment with kapecitabin or fluorouracil, 50% of normal start dose is recommended.',
    '*2A/HapB3'=
        'This patient carries two alleles associated with deficient and loss of DPYD activity, respectively.
        Severely deficient enzymatic activity can be expected.
        Treatment with kapecitabin or fluorouracil should be avoided.
        If attempted, carefulness and no more than 25% of normal start dose is recommended.',
    '*13/HapB3'=
        'This patient carries two alleles associated with deficient and loss of DPYD activity, respectively.
        Severely deficient enzymatic activity can be expected.
        Treatment with kapecitabin or fluorouracil should be avoided.
        If attempted, carefulness and no more than 25% of normal start dose is recommended.',
    '*2A/c.2846A>T'=
        'This patient carries two alleles associated with deficient and loss of DPYD activity, respectively.
        Severely deficient enzymatic activity can be expected.
        Treatment with kapecitabin or fluorouracil should be avoided.
        If attempted, carefulness and no more than 25% of normal start dose is recommended.',
    '*13/c.2846A>T'=
        'This patient carries two alleles associated with deficient and loss of DPYD activity, respectively.
        Severely deficient enzymatic activity can be expected.
        Treatment with kapecitabin or fluorouracil should be avoided.
        If attempted, carefulness and no more than 25% of normal start dose is recommended.',
    '*2A/*2A'=
        'This patient is homozygous for an allele associated with loss of DPYD activity.
        Minimal enzymatic activity can be expected.
        Treatment with kapecitabin or fluorouracil should be avoided.',
    '*13/*13'=
        'This patient is homozygous for an allele associated with loss of DPYD activity.
        Minimal enzymatic activity can be expected.
        Treatment with kapecitabin or fluorouracil should be avoided.',
    '*13/*2A'=
        'This patient carries two alleles associated with loss of DPYD activity.
        Minimal enzymatic activity can be expected.
        Treatment with kapecitabin or fluorouracil should be avoided.'
)



# Parse alleles -----------------------------------------------------------

{
    infile <- opt$bam_file

    # --- START: Reference Genome Detection ---

    # Check BAM header to identify the reference genome.
    header <- scanBamHeader(infile)
    header_text <- unlist(header[[1]]$text)

    # Define identifiers for hg19/GRCh37 and hg38/GRCh38.
    hg19_identifiers <- c("hg19", "GRCh37", "human_g1k_v37", "Homo_sapiens_assembly19")
    hg38_identifiers <- c("hg38", "GRCh38")

    # Detect which genome build is used in the BAM file.
    is_hg19 <- any(sapply(hg19_identifiers, function(id) any(grepl(id, header_text, ignore.case = TRUE))))
    is_hg38 <- any(sapply(hg38_identifiers, function(id) any(grepl(id, header_text, ignore.case = TRUE))))

    # Select the correct variant coordinates based on the detected genome.
    # The 'pos' column is dynamically created for the merge later.
    if (is_hg19) {
        variants <- copy(variants_table)
        variants[, pos := pos_hg19]
    } else if (is_hg38) {
        variants <- copy(variants_table)
        variants[, pos := pos_hg38]
    } else {
        stop("ERROR: BAM file reference genome is not recognized as hg19 or hg38.")
    }

    # Create the GRanges object for querying the BAM file.
    ranges <- makeGRangesFromDataFrame(unique(variants[, .(chromosome, start = pos, end = pos)]))

    # --- END: Reference Genome Detection ---

    p_params<- PileupParam(max_depth=2e5, min_base_quality= 20, min_mapq= 20,
                           min_nucleotide_depth=1, min_minor_allele_depth= 0,
                           distinguish_strands= FALSE, distinguish_nucleotides= TRUE,
                           ignore_query_Ns= FALSE, include_deletions= TRUE, include_insertions=
                               FALSE)
    b_params <- ScanBamParam(which = ranges)
    bf <- BamFile(infile, yieldSize=2e5)
    data <- as.data.table(pileup(bf, scanBamParam = b_params, PileupParam=p_params))


    # --- Chromosome Naming Fix ---

    # Standardize chromosome names from pileup output (e.g., "chr1" -> "1").
    if(nrow(data) > 0) {
        colnames(data)[1] <- 'chromosome'
        data[, chromosome := str_remove(chromosome, "chr")]
    } else {
        # If pileup returns no data, create an empty table with the correct columns.
        data <- data.table(chromosome=character(), pos=integer(), nucleotide=character(), count=integer())
    }


    # Merge pileup results with variant definitions. 'all=T' keeps all defined
    # variants, ensuring those with zero coverage are reported.
    data <- merge(variants, data[count>2], by=c('chromosome','pos','nucleotide'), all=T)
    data <- data.table(file=infile, data)
}


# Process --------------------------------------------------------

{
    # Calculate total read depth at each unique SNP position.
    data[,depth:=sum(count,na.rm=T),by=c('rsid')]

    # Aggregate counts for each nucleotide at each SNP position.
    data[,count:=sum(count),by=c('rsid','nucleotide')]

    # Strand information is not used, so remove the column.
    data[,strand:=NULL]
    data <- unique(data)

    # Replace NA counts with 0 for sites with no coverage.
    data[is.na(count), count := 0]
    data[,variant_ratio:=count/depth]
    # Handle cases with zero depth to avoid NaN (division by zero).
    data[depth == 0, variant_ratio := 0]

}

# Get genotype ----------------------------------------------------------------

# Define thresholds for calling alleles and flagging issues.
mindepth <- 50 # Minimum depth required to consider a site covered.
mincount <- 5 # Minimum number of reads for an allele to be called.
minratio <- .05 # Minimum allele fraction for an allele to be called.
homratio <- .95 # Allele fraction threshold to call homozygous.
low_warning <- .1 # Lower bound for VAF to warn about potential contamination.
high_warning <- .9 # Upper bound for VAF to warn about potential contamination.


# Determine how many of the 4 target SNPs have sufficient coverage.
data[,variants_covered:=length(unique(rsid[depth>mindepth]))]

# Flag rows corresponding to wild-type (wt) or alternate (alt) alleles
# if they pass the minimum count and ratio thresholds.
data[,is_wt:=0][count>=mincount & variant_ratio>minratio & allele=='wt',is_wt:=1]
data[,is_alt:=0][count>=mincount & variant_ratio>minratio & allele=='alt',is_alt:=1]

# Count the number of SNP sites with confident wt/alt calls.
data[,n_wt:=sum(is_wt), by=file]
data[,n_alt:=sum(is_alt), by=file]


# Determine the number of copies for each allele (0, 1, or 2).
data[,copies:=0]
data[count>=mincount & variant_ratio > minratio, copies:=1]
data[variant_ratio > homratio,copies:=2]

# Sum the total copies of wt and alt alleles across all 4 SNPs.
data[,wt_copies:=sum(copies*is_wt), by=file]
data[,alt_copies:=sum(copies*is_alt), by=file]

# --- Genotype assignment based on a hierarchy of rules ---

# Default genotype if no rules match.
data[,genotype:='Undefined']

# Rule for wild-type (*1/*1).
data[wt_copies==8 & alt_copies==0 & variants_covered==4,genotype:='*1/*1']

# Rule for a single homozygous variant.
data[variants_covered==4 & n_alt==1 & n_wt==3, genotype:=paste0(rep(unique(alias[is_alt==1 & copies==2]),2),collapse = '/')]

# Rule for two different heterozygous variants (compound het).
data[variants_covered==4 & n_alt==2 & n_wt==4, genotype:=paste0(sort(unique(alias[is_alt==1 & copies==1])),collapse = '/')]

# Rule for a single heterozygous variant.
data[variants_covered==4 & n_alt==1 & n_wt==4, genotype:=paste0('*1/',unique(alias[is_alt==1 & copies==1]))]

# Consolidate the genotype call for the sample.
this_genotype <- unique(data$genotype)

# Select the appropriate clinical text based on the final genotype.
this_genotype_text <- 'Analysis did not result in a clear genotype, review evidence table.'

if (length(this_genotype)>1) {
    this_genotype_text <- 'Analysis failed with apparently >1 genotype selected, review evidence table.'
    this_genotype <- 'Review needed'
} else {
    ix <- match(this_genotype,names(genotype_text))
    if (!is.na(ix)) this_genotype_text <- str_squish(genotype_text[[ix]])
}

# Coverage stat -------------------------------------------------------------

# Prepare a summary of sequencing depth for the report.
coverage <- unique(data[,.(rsid,alias,depth)])
coverage[,text:=paste0(rsid,'(',alias,'): ',depth,'x')]

if (all(coverage$depth>mindepth)) {
    text <- 'Relevant DPYD SNPs were sequenced to sufficient coverage:'
} else if (any(coverage$depth>mindepth)) {
    text <- 'Some DPYD SNPs were sequenced to sufficient coverage:'
} else {
    text <- 'DPYD SNPs were not sequenced to sufficient coverage:'
}

coverage_text <- paste(text,paste(coverage$text,collapse = ', '))


# Add warnings for allele fractions in ambiguous ranges (e.g., contamination).
if (any(data[allele=='alt']$variant_ratio > 0 & data[allele=='alt']$variant_ratio < low_warning, na.rm = TRUE)) {
    this_genotype_text <- paste(this_genotype_text,'WARNING: Low-fraction variant observed, review data to verify genotype')
}
if (any(data[allele=='alt']$variant_ratio > high_warning & data[allele=='alt']$variant_ratio < 1, na.rm = TRUE)) {
    this_genotype_text <- paste(this_genotype_text,'WARNING: Low-fraction reference observed, review data to verify genotype')
}

# Add texts to data ----------------------------------------------------------
data[,genotype_text:=this_genotype_text]
data[,coverage_text:=coverage_text]

# Write result ----------------------------------------------------------

out_dir='.'

csv_file <- paste0(out_dir,'/',basename(infile),'.DPYD.csv')
if (!is.null(opt$csv_file)) csv_file <- opt$csv_file

fwrite(data,file = csv_file)


json_file <- paste0(out_dir,'/',basename(infile),'.DPYD.json')
if (!is.null(opt$json_file)) json_file <- opt$json_file

exportJson <- toJSON(list('sample'=infile,'genotype'=this_genotype,'text'=this_genotype_text,'coverage'=coverage_text))
write(exportJson, json_file)

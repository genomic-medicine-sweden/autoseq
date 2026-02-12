#!/usr/bin/env python

import argparse
import pandas as pd
import logging

__version__ = "1.0.0"

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def annotate_cnvs(input_file: str, curation_ann_file: str, cancer_type: str, sample_type: str, output_file: str) -> None:
    """
    Annotates Copy Number Variant (CNV) segmentation data with gene-specific comments
    to facilitate clinical or research curation.

    The function maps genes found in CNV segments against a curated annotation database,
    filtering for the specific cancer project and sample origin (somatic vs. germline).

    Args:
        input_file (str): Path to the TSV file containing CNV segments.
                          Expected to have a 'gene' column (comma-separated list).
        curation_ann_file (str): Path to the CSV containing master gene annotations.
                                 Must include columns: 'gene', 'prefix', 'type', and 'comment'.
        cancer_type (str): The specific cancer abbreviation (e.g., 'BRCA', 'LUAD').
                           Used to derive the project prefix.
        sample_type (str): Origin of the sample. Values 'tumor' or 'cfdna' are treated
                           as 'somatic'; others are treated as 'germline'.
        output_file (str): Path where the final annotated TSV will be saved.

    Returns:
        None: Writes the resulting DataFrame directly to output_file.
    """

    logging.info(f"Starting annotation of CNVs from {input_file} using {curation_ann_file}")
    cnvs_df = pd.read_csv(input_file, sep='\t')

    logging.info(f"Read {len(cnvs_df)} rows from input file")
    ann_df = pd.read_csv(curation_ann_file)

    logging.info(f"cancer type code - {cancer_type}")
    ## filter the annotation file based on the project and type
    ann_df = ann_df[(ann_df.prefix == cancer_type) & (ann_df.type == sample_type)]
    logging.info(f"Filtered annotation file to {len(ann_df)} rows for project {cancer_type} and type {sample_type}")

    # Create a dictionary from the annotation file for quick lookup
    annotation_dict = ann_df.set_index('gene')['comment'].to_dict()

    # Function to annotate a row
    def annotate_row(row) -> pd.Series:
        genes = row['gene'].split(',') if not pd.isnull(row['gene']) else []
        comments = [" ".join([gene, annotation_dict[gene]]) for gene in genes if gene in annotation_dict]
        row['curate'] = ','.join(comments)
        return row

    # Apply the annotation function to each row
    cnvs_df = cnvs_df.apply(annotate_row, axis=1)
    logging.info("Applied annotations to CNV data")

    # Write the updated CNV segmentations to the output file
    cnvs_df.to_csv(output_file, sep='\t', index=False)
    logging.info(f"Written annotated CNV data to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=
        'Annotate CNVs with cancer gene specific information for curation')
    parser.add_argument('-i', '--input', required=True, help="Input CNV segmentations file")
    parser.add_argument('-t', '--sample-type', choices= ['germline', 'somatic'], required=True,
                        help="Sample type: somatic or germline")
    parser.add_argument('-c', '--curation-annotation-file', help="Annotation file for curation")
    parser.add_argument('--cancer-type', default='PANCANCER', help="shorthand code for cancer type; Default: PANCANCER")
    parser.add_argument('-o', '--output', required=True,
                        help="output segmentation file with curate column in the end")
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s {__version__}')
    args = parser.parse_args()

    annotate_cnvs(args.input, args.curation_annotation_file, args.cancer_type, args.sample_type, args.output)

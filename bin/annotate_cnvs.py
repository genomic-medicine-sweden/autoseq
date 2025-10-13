#!/usr/bin/env python

import os
import argparse
import pandas as pd
import logging

__version__ = "1.0.0"

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def annotate_cnvs(input_file, curation_ann_file, output_file):
    """
    Annotate CNVs with cancer gene specific information for curation

    :param input_file: Input CNV segmentations file
    :param curation_ann_file: Annotation file for curation
    :param output_file: Output CNV segmentations file with curate column in the end
    :param type: Sample Type - Normal, Tumor, CFDNA to be added
    :return
    """

    logging.info(f"Starting annotation of CNVs from {input_file} using {curation_ann_file}")
    cnvs_df = pd.read_csv(input_file, sep='\t')

    logging.info(f"Read {len(cnvs_df)} rows from input file")
    ann_df = pd.read_csv(curation_ann_file)

    logging.info(f"Read {len(ann_df)} rows from annotation file")
    project = os.path.basename(input_file).split('-')[0]
    sample_type = os.path.basename(input_file).split('-')[3]
    type = 'somatic' if sample_type in ['T', 'CFDNA'] else 'germline'

    logging.info(f"Project identified as {project}")
    ## filter the annotation file based on the project and type
    ann_df = ann_df[(ann_df.prefix == project) & (ann_df.type == type)]
    logging.info(f"Filtered annotation file to {len(ann_df)} rows for project {project} and type {type}")

    # Create a dictionary from the annotation file for quick lookup
    annotation_dict = ann_df.set_index('gene')['comment'].to_dict()


    # Function to annotate a row
    def annotate_row(row):
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
    parser.add_argument('-c', '--curation-ann', help="Annotation file for curation")
    parser.add_argument('-o', '--output', required=True,
                        help="output segmentation file with curate column in the end")
    parser.add_argument('-v', '--version', action='version', version=f'%(prog)s {__version__}')
    args = parser.parse_args()

    annotate_cnvs(args.input, args.curation_ann, args.output)

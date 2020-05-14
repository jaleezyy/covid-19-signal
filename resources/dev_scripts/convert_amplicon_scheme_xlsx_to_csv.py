#!/usr/bin/env python
#-*- coding: utf-8 -*-

import pandas as pd # needs xlrd to support xlsx
import argparse
import sys, os


def convert_xlsx_to_bed(input_xlsx, output_filepath,
			reference_name):
    """
    Convert xlsx formatted amplicon scheme specific to a BED file that can be
    used by ivar etc to remove primers
    """
    # expects "Primer Name", "Pool", "Start", "End"
    df = pd.read_excel(input_xlsx)

    # check fields exist
    expected_fields = ["Primer Name", "Pool", "Start", "End"]
    for field in expected_fields:
        if field not in df.columns:
            print(f"Required field '{field}' missing in '{input_xlsx}'")
            print("Please add this column and try again")
            sys.exit(1)

    # add stranding
    df.loc[df['Start'] < df['End'], 'strand'] = '+'
    df.loc[df['Start'] > df['End'], 'strand'] = '-'

    # reorder start end so they are in the same orientation
    df['bed_start'] = df[['Start', 'End']].min(axis=1)
    df['bed_end'] = df[['Start', 'End']].max(axis=1)

    # reference name
    df['reference'] = reference_name

    # rename pool and primer name column for consistency
    df = df.rename(columns={'Pool': 'bed_poolname',
                            'Primer Name': 'primer_name'})

    # extract bed columns from dataframe
    bed  = df[['reference', 'bed_start', 'bed_end',
		'primer_name', 'bed_poolname', 'strand']]

    # save to bed file without index or column names
    bed.to_csv(output_filepath, sep='\t', header=False, index=False)

    print(f"BED formatted scheme file saved to {output_filepath} with "
          f"reference: '{reference_name}'")
    return


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Process some integers.')

    parser.add_argument('input_xlsx',
                        help="Path to XLSX primer scheme "
                             "Expects the following fields: "
		             "'Primer Name', 'Pool', 'Start', 'End'")
    parser.add_argument('-r', '--reference_name', default="MN908947.3",
                        required=False,
                        help='Accession of reference genome used')
    parser.add_argument('-o', '--output_path', default="$INPUT_XLSX_PATH.bed",
			required=False,
			help='Filename/path for generated .bed file')

    args = parser.parse_args()

    # similarly if output path not specified, just input with bed suffix
    if args.output_path == "$INPUT_XLSX_PATH.bed":
        args.output_path = os.path.splitext(args.input_xlsx)[0] + '.bed'

    convert_xlsx_to_bed(args.input_xlsx, args.output_path,
			args.reference_name)

#!/usr/bin/env python
import pysam
import sys
import argparse

def filter_reads(contig_name, input_sam_fp, output_bam_fp):

    # use streams if args are None
    if input_sam_fp:
        input_sam = pysam.AlignmentFile(input_sam_fp, 'r')
    else:
        input_sam = pysam.AlignmentFile('-', 'r')

    if output_bam_fp:
        output_bam = pysam.AlignmentFile(output_bam_fp, 'wb',
                                         template=input_sam)
    else:
        output_bam = pysam.AlignmentFile('-', 'wb', template=input_sam)

    # if read isn't mapped or mapped to viral reference contig name
    viral_reads = 0
    human_reads = 0
    unmapped_reads = 0

    # iterate over input from BWA
    for read in input_sam:
        # only look at primary alignments
        if not read.is_supplementary and not read.is_secondary:
            if read.reference_name == contig_name:
                output_bam.write(read)
                viral_reads += 1
            elif read.is_unmapped:
                output_bam.write(read)
                unmapped_reads +=1
            else:
                human_reads += 1

    total_reads = viral_reads + human_reads + unmapped_reads

    print(f"viral read count = {viral_reads} "
            f"({viral_reads/total_reads * 100:.2f}%)", file=sys.stderr)
    print(f"human read count = {human_reads} "
            f"({human_reads/total_reads * 100:.2f}%) ", file=sys.stderr)
    print(f"unmapped read count = {unmapped_reads} "
            f"({unmapped_reads/total_reads * 100:.2f}%)", file=sys.stderr)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Get all reads that are "
                                                  "unmapped and map to a "
                                                  "specific reference "
                                                  "contig")
    parser.add_argument('-i', '--input', required=False, default=False,
                        help="Input SAM formatted file (stdin used if "
                              " not specified)")

    parser.add_argument('-o', '--output', required=False, default=False,
                        help="Output BAM formatted file (stdout used if not "
                             "specified)")

    parser.add_argument('-c', '--contig_name', required=False,
                        default="MN908947.3",
                        help="Contig name to retain e.g. viral")

    args = parser.parse_args()

    filter_reads(args.contig_name, args.input, args.output)

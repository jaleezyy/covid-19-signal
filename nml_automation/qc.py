#!/usr/bin/env python3

from Bio import SeqIO
from functools import reduce
from pybedtools import BedTool
import csv
import subprocess
import re
import pandas as pd
import matplotlib.pyplot as plt
import shlex

"""
This script can incorporate as many QC checks as required
as long as it outputs a csv file containing a final column
headed with 'qc_pass' and rows for each sample indcating
'TRUE' if the overall QC check has passed or 'FALSE' if not.
"""

def make_qc_plot(depth_pos, n_density, samplename, window=200):
    depth_df = pd.DataFrame( { 'position' : [pos[1] for pos in depth_pos], 'depth' : [dep[2] for dep in depth_pos] } )
    depth_df['depth_moving_average'] = depth_df.iloc[:,1].rolling(window=window).mean()

    n_df = pd.DataFrame( { 'position' : [pos[0] for pos in n_density], 'n_density' : [dens[1] for dens in n_density] } )

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    ax1.set_xlabel('Position')

    ax1.set_ylabel('Depth', color = 'g')
    ax1.set_ylim(top=10**5, bottom=1)
    ax1.set_yscale('log')
    ax1.plot(depth_df['depth_moving_average'], color = 'g')

    ax2.set_ylabel('N density', color = 'r')  
    ax2.plot(n_df['n_density'], color = 'r')
    ax2.set_ylim(top=1)

    plt.title(samplename)
    plt.savefig(samplename + '.depth.png')

def read_depth_file(bamfile):
    p = subprocess.Popen(['samtools', 'depth', '-a', '-d', '0', bamfile],
                       stdout=subprocess.PIPE)
    out, err = p.communicate()
    counter = 0

    pos_depth = []
    for ln in out.decode('utf-8').split("\n"):
       if ln:
          pos_depth.append(ln.split("\t"))

    return pos_depth


def get_covered_pos(pos_depth, min_depth):
    counter = 0
    for contig, pos,depth in pos_depth:
        if int(depth) >= min_depth:
            counter = counter + 1
    
    return counter

def get_depth_coverage(pos_depth, ref_length):
    depth_total = 0
    for contig, pos, depth in pos_depth:
        depth_total = depth_total + int(depth)

    return round(depth_total/ref_length)

def get_N_positions(fasta):
    n_pos =  [i for i, letter in enumerate(fasta.seq.lower()) if letter == 'n']

    return n_pos

def get_pct_N_bases(fasta):
    
    count_N = len(get_N_positions(fasta))

    pct_N_bases = count_N / len(fasta.seq) * 100

    return pct_N_bases, count_N

def get_largest_N_gap(fasta):
    n_pos = get_N_positions(fasta)

    n_pos = [0] + n_pos + [len(fasta.seq)]

    n_gaps = [j-i for i, j in zip(n_pos[:-1], n_pos[1:])]

    return sorted(n_gaps)[-1]

def get_ref_length(ref):
    record = SeqIO.read(ref, "fasta")
    return len(record.seq)

def sliding_window_N_density(sequence, window=10):

    sliding_window_n_density = []
    for i in range(0, len(sequence.seq), 1):
        window_mid = i + ( window / 2)
        window_seq = sequence.seq[i:i+window]
        n_count = window_seq.lower().count('n')
        n_density = n_count / window

        sliding_window_n_density.append( [ window_mid, n_density ] )

    return sliding_window_n_density

def get_num_reads(bamfile):

    st_filter = '0x900'
    command = 'samtools view -c -F{} {}'.format(st_filter, bamfile)
    what = shlex.split(command)

    return subprocess.check_output(what).decode().strip()

def get_variants(variants_tsv, variants_list=[], locations=[]):
    with open(variants_tsv) as input_handle:

        for index, line in enumerate(input_handle):

            # Pass the header line so as to not include it!
            if index == 0:
                continue

            row = line.strip('\n').split('\t') # Order is [REGION, POS, REF, ALT, ...]

            variant = '{}{}{}'.format(row[2], row[1], row[3])

            # Checking for duplicate variants that have been an issue
            if variant in variants_list:
                pass
            else:
                variants_list.append(variant)
                locations.append(row[1])
    
    variants = (';'.join(variants_list))

    if variants == '':
        return 'None', None
        
    return variants, locations

def find_primer_mutations(pcr_bed, variants_locations, primer_mutations=[]):
    if variants_locations is None:
        return 'None'
    
    input_bed = BedTool(pcr_bed)

    for gene in input_bed:
        location = range(gene.start, gene.stop + 1) # Plus one to make sure that we get mutations in the final location of the range

        for variant_pos in variants_locations:
            if variant_pos in location:
                primer_mutations.append('Variant position {} overlaps PCR primer {}'.format(variant_pos, gene.name))
    
    if primer_mutations != []:
        statement = '; '.join(primer_mutations)
        return 'Warning: {}'.format(statement)

    return 'None'

def get_lineage(pangolin_csv, sample_name):
    with open(pangolin_csv, 'r') as input_handle:
        reader = csv.reader(input_handle)

        for row in reader: # Row format is ['taxon', 'lineage', 'SH-alrt', 'UFbootstrap', 'lineages_version', 'status', 'note']

            if re.search(sample_name, row[0]):
                return str(row[1])
    
    return 'Unknown'

def parse_ncov_tsv(file_in, sample, negative=False):

    # Try to read file (as negative control may not have data in it)
    try:
        df = pd.read_csv(file_in, sep='\t')

    # If no data, we set up how it should be and then pass it through
    # Could also make is such that runs without negative ctrls just don't have the columns
    except pd.errors.EmptyDataError:
        negative_df = pd.DataFrame(columns=['sample', 'qc', 'genome_covered_bases', 'genome_total_bases', 'genome_covered_fraction', 'amplicons_detected'])
        negative_df.loc[1, 'sample'] = sample
        negative_df.fillna('NA', inplace=True)
    
        return negative_df

    # If these get changed in the input just replace the new file name here
    # First column is the file name
    if negative:
        new_columns = df.columns.values
        new_columns[0] = 'sample'
        df.columns = new_columns
    else:
        # Input is summary_df, drop its lineage column as we pull and create our own
        df.drop(columns=['lineage'], inplace=True)

    file_column = 'sample'

    for index, name in enumerate(df[file_column].tolist()):
        if re.search(sample, name):
            df.loc[index, file_column] = sample
            df.fillna('NA', inplace=True)
            return df.iloc[[index]]
    
    # If sample is not a negative control need to keep columns
    negative_df = pd.DataFrame(columns=new_columns)
    negative_df.loc[1, file_column] = sample
    negative_df.fillna('NA', inplace=True)
    
    return negative_df

def get_samplesheet_info(sample_tsv, sample_name):
    with open(sample_tsv) as input_handle:

        for line in input_handle:
            row = line.strip('\n').split('\t') # Order is [sample, run, barcode, project_id, ct]

            if re.search(sample_name, row[0]):
                return str(row[1]), str(row[2]), str(row[3]), str(row[4]) # [run, barcode, project_id, ct]

    return 'Unknown', 'Unknown', 'Unknown', 'Unknown'
    
def go(args):
    if args.illumina:
        depth = 10
    elif args.nanopore:
        depth = 20

    ## Depth calcs
    ref_length = get_ref_length(args.ref)
    depth_pos = read_depth_file(args.bam)

    depth_covered_bases = get_covered_pos(depth_pos, depth)

    depth_coverage = get_depth_coverage(depth_pos, ref_length)

    pct_covered_bases = depth_covered_bases / ref_length * 100

    ## Number of aligned reads calculaton
    num_reads = get_num_reads(args.bam)

    # Unknown base calcs
    fasta = SeqIO.read(args.fasta, "fasta")

    pct_N_bases   = 0
    largest_N_gap = 0
    qc_pass       = "FALSE"

    if len(fasta.seq) != 0:

        pct_N_bases, count_N = get_pct_N_bases(fasta)
        largest_N_gap = get_largest_N_gap(fasta)

    	# QC PASS / FAIL
        if largest_N_gap >= 10000 or pct_N_bases < 50.0:
                qc_pass = "TRUE"

    ## Added checks ##
    # TSV Variants
    variants, variant_locations = get_variants(args.variants)

    # Find any overlap of variants in the pcr primer regions
    primer_statement = find_primer_mutations(args.pcr_bed, variant_locations)

    # Pangolin Lineages
    lineage = get_lineage(args.pangolin, args.sample)

    # NCOV-Tools Results
    summary_df = parse_ncov_tsv(args.ncov_summary, args.sample)
    negative_df = parse_ncov_tsv(args.ncov_negative, args.sample, negative=True)

    # No samplesheet used for illumina stuff really
    if args.sample_sheet:
        run_name, barcode, project_id, ct = get_samplesheet_info(args.sample_sheet, args.sample)
    
    # IRIDA project is 2 for illumina and 3 for nanopore
    else:
        run_name = 'NA'
        project_id = 2
        ct = 'NA'


    qc_line = {      'sample' : [args.sample],
                 'project_id' : [project_id],
           'num_aligned_reads': [num_reads],
                    'lineage' : [lineage],
                   'variants' : [variants],
'diagnostic_primer_mutations' : [primer_statement],
                'script_name' : ['covid-signal'],
                   'revision' : [args.revision]}

    qc_df = pd.DataFrame.from_dict(qc_line)

    data_frames = [qc_df, summary_df, negative_df]

    # Merge all dataframes together
    out_df = reduce(lambda left,right: pd.merge(left,right,on='sample', how='left'), data_frames)

    # Remove comma's as some of the ncov-tools fields have commas :(
    out_df.replace(',',';', regex=True, inplace=True)

    # Output
    out_df.to_csv(args.outfile, sep=',', index=False)
    N_density = sliding_window_N_density(fasta)
    make_qc_plot(depth_pos, N_density, args.sample)

def main():
    import argparse

    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--nanopore', action='store_true')
    group.add_argument('--illumina', action='store_true')
    parser.add_argument('--outfile', required=True)
    parser.add_argument('--sample', required=True)
    parser.add_argument('--ref', required=True)
    parser.add_argument('--bam', required=True)
    parser.add_argument('--fasta', required=True)
    parser.add_argument('--pangolin', required=True)
    parser.add_argument('--ncov_summary', required=True)
    parser.add_argument('--ncov_negative', required=True)
    parser.add_argument('--revision', required=True)
    parser.add_argument('--variants', required=True)
    parser.add_argument('--pcr_bed', required=True)
    parser.add_argument('--sample_sheet', required=False)

    args = parser.parse_args()
    go(args)

if __name__ == "__main__":
    main()

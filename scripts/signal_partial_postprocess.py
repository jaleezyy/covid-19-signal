#!/bin/env python

import os, sys
import argparse
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
from pathlib import Path
from collections import defaultdict

def check_file(path: str) -> Path:
	"""
	Check an input file exists and is readable
	"""
	path = Path(path)
	if path.exists() and path.is_file():
		return path
	else:
		return False

def parse_lineages(file):
	lineages = pd.read_table(file, sep='\t')
	try:
		lin_df = lineages[['isolate', 
							'pango_lineage', 
							'pangolin_version', 
							'pangoLEARN_version',
							'pango_version',
							'pangolin_qc', 
							'pangolin_note', 
							]]
	except KeyError: # Pangolin v4 modifications
		lin_df = lineages[['isolate', 
							'pango_lineage', 
							'pangolin_version', 
							'version',
							'scorpio_version',
							'constellation_version',
							'pangolin_qc',
							'pangolin_qc_note',
							'pangolin_note', 
							]]
						
	lin_df['isolate'] = lin_df['isolate'].apply(lambda x: x.split("_")[1].split(".")[0] if x.startswith("Consensus") else x)
	return lin_df
	
def parse_stats(stats):
	info = defaultdict(list)
	missing = set()
	for file in stats:
		sample_name = os.path.basename(file).split("_")[0] # intended sample name
		if check_file(file) is False: # if missing or incomplete
			print(f"Missing {file}")
			missing.add(f"{sample_name}")
			continue
		else:
			with open(file) as data:
				file_red = data.readlines()
			if any(s.startswith("Assembly\t") and (match := s) for s in file_red):
				info['isolate'].append(match.split("\t")[1].split(".")[0].strip("\n")) ### QUAST ONLY
			else:
				info['isolate'].append(f"{sample_name}") ### DEFAULT
			### QUAST
			if any(s.startswith("Genome fraction (%)\t") for s in file_red): # and (match := s) for s in file_red):
				#info["Genome Fraction (%)"].append(match.split("\t")[1].strip("\n")) ### QUAST
				if any(t.startswith("Total aligned length\t") and (match_aligned := t) for t in file_red): # define as the total number of aligned bases in the assembly
					covered = int(match_aligned.split("\t")[1].strip("\n"))
				else:
					covered = 0
				if any(r.startswith("Reference length\t") and (match_len := r) for r in file_red):
					length = int(match_len.split("\t")[1].strip("\n"))
				else:
					length = 1
				genome_frac = round((covered/length)*100, 2)
				
				info["Genome Fraction (%)"].append(f"{genome_frac}")
			### STATS
			elif any(s.startswith("[Alignment Statistics]") for s in file_red):
				if any(c.startswith("Covered Bases: ") and (match_cov := c) for c in file_red):
					covered = int(match_cov.split(": ")[1].strip("\n"))
				else:
					covered = 0
				if any(r.startswith("Reference Length: ") and (match_len := r) for r in file_red):
					length = int(match_len.split(": ")[1].strip("\n"))
				else:
					length = 1
				genome_frac = round((covered/length)*100, 2)
				
				info["Genome Fraction (%)"].append(f"{genome_frac}")
			### DEFAULT
			else:
				info["Genome Fraction (%)"].append(" ")

			quast_df = pd.DataFrame(info, columns=['isolate', 'Genome Fraction (%)'])
		
	if len(info) == 0 and len(missing) > 0: # nothing was added after each file, put default values
		for sample_name in missing:
			info['isolate'].append(f"{sample_name}")
			info["Genome Fraction (%)"].append(" ")
			quast_df = pd.DataFrame(info, columns=['isolate', 'Genome Fraction (%)'])

		
	return quast_df
	
def collate_output(lineage, stats, output):
	merged_df = lineage.merge(stats, on='isolate', how='outer')
	merged_df.rename(columns={'isolate': 'SAMPLE ID', 'pango_lineage': "LINEAGE", 'pangolin_qc': 'STATUS', 'pangolin_qc_note': "QC NOTE", 'pangolin_note': "NOTE"}, inplace=True)
	merged_df.to_csv(output, sep='\t', index=False)


if __name__ == '__main__':
	parser = argparse.ArgumentParser('Partial postprocessing script for SIGNAL using input VCFs or consensus FASTAs')
	parser.add_argument("-i", "--input_lineages", required=True, help="Lineage assessment file from SIGNAL")
	parser.add_argument("-q", "--input_sample", nargs='*', required=True, help="TXT file stats, QUAST TSV output, or list of sample IDs")
	parser.add_argument("-o", "--output", default="summary.tsv", help="Output file for collated summary table")
	args = parser.parse_args()
	
	if check_file(args.input_lineages) is False:
		raise argparse.ArgumentTypeError(f"Required lineage assessment file {args.input_lineages} can't be read!")
	lin = parse_lineages(args.input_lineages)
	
	stat = parse_stats(args.input_sample)
	
	collate_output(lin, stat, args.output)

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
		#raise argparse.ArgumentTypeError(f"{path} can't be read")

def identify_stats_file(id):
	"""
	Determine if provided stats file is from QUAST or already provided
	"""

	files = []
	
	for filename in id:
		stats = False
		quast = False
	
		stats_file = os.path.join(snakemake.params['results_dir'], filename, "_stats.txt")
		quast_file = os.path.join(snakemake.params['results_dir'], filename, "_quast_report.tsv")
		
		if check_file(stats_file): stats = True
		if check_file(quast_file): quast = True
	
		if (stats is True) and (quast is True): # select stats file by default
			files.append(stats_file)
			#return stats
		elif (stats is True) and (quast is False):
			files.append(stats_file)
			#return stats
		elif (stats is False) and (quast is True):
			files.append(quast_file)
			#return quast
		else:
			files.append(stats_file)
			#return stats
	
	return files

def parse_lineages(file):
	lineages = pd.read_table(file, sep='\t')
	lin_df = lineages[['isolate', 
						'pango_lineage', 
						'pangolin_conflict', 'pangolin_ambiguity_score', 
						'pangolin_qc', 
						'pangolin_note', 
						'pangolin_version', 
						'pangoLEARN_version',
						'pango_version'
						]]
	lin_df.rename(columns={'pangolin_qc': 'status'}, inplace=True)
	lin_df['isolate'] = lin_df['isolate'].apply(lambda x: x.split("_")[1].split(".")[0] if x.startswith("Consensus") else x)
	return lin_df
	
def parse_stats(stats):
	#quast_df = pd.DataFrame(columns=['isolate', 'Genome Fraction (%)'])
	info = defaultdict(list)
	for file in stats:
		sample_name = os.path.basename(file).split("_")[0] # intended sample name
		if check_file(file) is False: # if missing or incomplete, replace with default values
			print(f"Adding {file}")
			info['isolate'].append(f"{sample_name}")
			info["Genome Fraction (%)"].append(" ")
			quast_df = pd.DataFrame(info, columns=['isolate', 'Genome Fraction (%)'])
		else:
			with open(file) as data:
				file_red = data.readlines()
				#print(file_red)
				# exit()
			if any(s.startswith("Assembly\t") and (match := s) for s in file_red):
				info['isolate'].append(match.split("\t")[1].split(".")[0].strip("\n")) ### QUAST
			else:
				info['isolate'].append(f"{sample_name}") ### DEFAULT
			
			if any(s.startswith("Genome fraction (%)\t") and (match := s) for s in file_red):
				info["Genome Fraction (%)"].append(match.split("\t")[1].strip("\n")) ### QUAST
			elif any(c.startswith("Covered Bases: ") and (match_cov := c) for c in file_red):
				covered = int(match_cov.split(": ")[1].strip("\n")) ### STATS
				if any(r.startswith("Reference Length: ") and (match_len := r) for r in file_red):
					length = int(match_len.split(": ")[1].strip("\n")) ### STATS
				genome_frac = round((covered/length)*100, 2)
				
				info["Genome Fraction (%)"].append(f"{genome_frac}")
			else:
				info["Genome Fraction (%)"].append(" ") ### DEFAULT


		# for ele in file_red:
			# if ele.startswith("Assembly\t"):
				# info['isolate'].append(ele.split("\t")[1].split(".")[0].strip("\n"))
			# if ele.startswith('Genome fraction (%)\t'): 
				# info["Genome Fraction (%)"].append(ele.split("\t")[1].strip("\n"))
			
	
	quast_df = pd.DataFrame(info, columns=['isolate', 'Genome Fraction (%)'])

		
	return quast_df
	
	
	
def collate_output(lineage, stats, output):
	merged_df = lineage.merge(stats, on='isolate', how='outer')
	
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
	
	stats_files = identify_stats_file(args.input_sample)
	stat = parse_stats(stats_files)
	
	collate_output(lin, stat, args.output)

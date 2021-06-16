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
	
def parse_quast(quasts):
	#quast_df = pd.DataFrame(columns=['isolate', 'Genome Fraction (%)'])
	info = defaultdict(list)
	for file in quasts:
		sample_name = os.path.basename(file).split("_")[0] # intended sample name
		if check_file(file) is False: # if missing or incomplete, replace with default values
			print(f"QUAST report {file} can't be read! Skipping!")
			info['isolate'].append(f"{sample_name}")
			info["Genome Fraction (%)"].append(" ")
			quast_df = pd.DataFrame(info, columns=['isolate', 'Genome Fraction (%)'])
			return quast_df
		
		with open(file) as data:
			file_red = data.readlines()
			#print(file_red)
			# exit()
		if any(s.startswith("Assembly\t") and (match := s) for s in file_red):
			info['isolate'].append(match.split("\t")[1].split(".")[0].strip("\n"))
		else:
			info['isolate'].append(f"{sample_name}")
		
		if any(s.startswith("Genome fraction (%)\t") and (match := s) for s in file_red):
			info["Genome Fraction (%)"].append(match.split("\t")[1].strip("\n"))
		else:
			info["Genome Fraction (%)"].append(" ")

		# for ele in file_red:
			# if ele.startswith("Assembly\t"):
				# info['isolate'].append(ele.split("\t")[1].split(".")[0].strip("\n"))
			# if ele.startswith('Genome fraction (%)\t'): 
				# info["Genome Fraction (%)"].append(ele.split("\t")[1].strip("\n"))
			
	
	quast_df = pd.DataFrame(info, columns=['isolate', 'Genome Fraction (%)'])

		
	return quast_df
	
	
	
def collate_output(lineage, quast, output):
	merged_df = lineage.merge(quast, on='isolate', how='outer')
	
	merged_df.to_csv(output, sep='\t', index=False)


if __name__ == '__main__':
	parser = argparse.ArgumentParser('Partial postprocessing script for SIGNAL using input consensus FASTAs')
	parser.add_argument("-i", "--input_lineages", required=True, help="Lineage assessment file from SIGNAL")
	parser.add_argument("-q", "--input_quast", nargs='*', required=True, help="QUAST TSV output")
	parser.add_argument("-o", "--output", default="summary.tsv", help="Output file for collated summary table")
	args = parser.parse_args()
	
	if check_file(args.input_lineages) is False:
		raise argparse.ArgumentTypeError(f"Required lineage assessment file {args.input_lineages} can't be read!")
	lin = parse_lineages(args.input_lineages)
	
	qua = parse_quast(args.input_quast)
	
	collate_output(lin, qua, args.output)

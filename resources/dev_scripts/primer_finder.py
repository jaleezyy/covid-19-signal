#Run in directory with fastqs

import sys
import csv
import os
import glob
import ast 
from Bio import SeqIO
import re
from collections import Counter
import glob, argparse

"""
Search for Rev primers in R1 and forward primers in R2
"""

def counter(fasta_file,read,output):
	for file in sorted(glob.glob("{}.fastq".format(read))):
		sample=(file.replace(".fastq",""))
		primer_seq_list=[]
		for seq_record in SeqIO.parse(file, "fastq"):
			fastq_header = str(seq_record.id)
			fastq_sequence = str(seq_record.seq).upper()
			with open(fasta_file, 'r') as primers:
				for seq_record in SeqIO.parse(primers, "fasta"):
					primer_header = str(seq_record.id)
					primer_sequence = str(seq_record.seq).upper()
					if primer_sequence in fastq_sequence:
						primer_seq_list.append(primer_sequence)

		sequence_counts = Counter(primer_seq_list)
		for s in sequence_counts:
			print("{}\t{}\t{}".format(s,sequence_counts[s],read))
			output.write("{}\t{}\t{}".format(s,sequence_counts[s],read) + "\n")
	

def main(args):
	out = open(args.output, 'w+')
	counter(args.forward_primers, "R2", out)
	counter(args.reverse_primers, "R1", out)
	out.close()

def run():
	parser = argparse.ArgumentParser(description='Count the instances of primers from given forward and reverse primer FASTAs on FASTQs within the current working directory.')
	parser.add_argument('-1','--forward_primers', dest="forward_primers", default=None, required=True, help='fasta forward primers file')
	parser.add_argument('-2','--reverse_primers', dest="reverse_primers", default=None, required=True, help='fasta reverse primers file')
	parser.add_argument('-o','--output', dest="output", default="count.txt", required=False, help='output file name. "count.txt" is the default.')
	args = parser.parse_args()
	main(args)  

if __name__ == '__main__':
	run()



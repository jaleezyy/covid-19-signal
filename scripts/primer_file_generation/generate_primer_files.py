#!/bin/python

import sys
from Bio import SeqIO

try:
	input = sys.argv[1]
except IndexError:
	print("Single argument required leading to a .tsv file")
	exit()

output = str(input).rsplit(".", 1)[0].rsplit("/",1)[0]

with open(output+".fasta", 'w+') as out:
	with open(input) as infile:
		next(infile)
		for line in infile:
			out.write(">" + str(line).split("\t")[0] + "\n" + str(line).split("\t")[2] + "\n")

fw = open(output+"_fw.fasta", 'w+')
rc = open(output+"_rc.fasta", 'w+')

for seq_record in SeqIO.parse(output+".fasta", 'fasta'):
	if "_left" in str(seq_record.description).lower():
		fw.write(">" + str(seq_record.description).rsplit("_",1)[0] + "\n" + str(seq_record.seq).upper() + "\n")
	elif "_right" in str(seq_record.description).lower():
		rc.write(">" + str(seq_record.description).rsplit("_",1)[0] + "_rc\n" + str(seq_record.seq).upper() + "\n")
	else:
		print("Unable to determine direction! Skipping.")
		continue

fw.close()
rc.close()

exit("\nDone.\n")
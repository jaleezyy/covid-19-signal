#!/bin/python

### Separate FC and RC primers
### Can be used with cutadapt to trim

from Bio import SeqIO
import sys

try:
	input_fasta = sys.argv[1]
	name = sys.argv[2]
except IndexError:
	exit("Input FASTA to divide. And basename of output files.")

fw = open(name+"_fw.fasta", 'w+')
rc = open(name+"_rc.fasta", 'w+')

for seq_record in SeqIO.parse(input_fasta, 'fasta'):
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

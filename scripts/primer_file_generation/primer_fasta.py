#!/bin/python

import sys

try:
	input = sys.argv[1]
except IndexError:
	print("Single argument required leading to a .tsv file")
	exit()

output = str(input).rsplit(".", 1)[0].rsplit("/",1)[1] + ".fasta"

with open(output, 'w+') as out:
	with open(input) as infile:
		next(infile)
		for line in infile:
			out.write(">" + str(line).split("\t")[0] + "\n" + str(line).split("\t")[2] + "\n")

exit("\nDone\n")

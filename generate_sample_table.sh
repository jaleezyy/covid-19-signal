#!/bin/bash

set -e # exit if pipeline returns non-zero status
set -o pipefail # return value of last command to exit with non-zero status

database_dir=0
name="sample_table.csv"
existing=0

HELP="""
ASSUMES FASTQ FILES ARE NAMED AS <sample_name>_*_R1_*.fastq(.gz)

Flags:
    -d  :  Path to directory containing fastq(.gz) (Absolute paths preferred for consistency, but can use relative paths)
    -n  :  Name of final sample table (with extension) (default: 'sample_table.csv') - will overwrite if file exists
    -e  :  Name of an existing sample table - will append to the end of the provided table

Select one of '-n' (new sample table) or '-e' (existing sample table). 
If neither provided, a new sample table called 'sample_table.csv' will be created (or overwritten) by default. 
"""

while getopts ":d:n:e:" option; do
    case "${option}" in
        d) database_dir=$OPTARG;;
        n) name=$OPTARG;;
	e) existing=$OPTARG;;
    esac
done

if [ $database_dir = 0 ] ; then
    echo "You must specify a data directory containing fastq(.gz) reads."
    echo "$HELP"
    exit 1
fi

if [ $existing = 0 ] ; then
	echo -e "Creating new sample table called ${name}\n"
 	echo "sample,r1_path,r2_path" > .${name} && echo "sample1,r1_path,r2_path"
else
	name=$(basename $existing)
	if [ -f ""$existing"" ] ; then
		echo -e "Using existing sample table called ${name}\n"
		cat $existing > .${name}
	else
		echo -e "Sample table does not exist. Check that sample table exists or create a new sample table."
		echo "$HELP"
		exit 1
	fi
fi

for file in $database_dir/*.fastq*; do
	sample=$(basename $file | cut -d_ -f 1)
	r1=$(ls $database_dir/$(basename $file | cut -d_ -f1)*R1* | grep /${sample}_)
	r2=$(ls $database_dir/$(basename $file | cut -d_ -f1)*R2* | grep /${sample}_)
	echo ${sample},${r1},${r2} >> .${name} && echo ${sample},${r1},${r2}
done

uniq .${name} > $name
rm .${name}

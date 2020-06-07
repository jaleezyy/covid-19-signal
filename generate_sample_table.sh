#!/bin/bash

set -e # exit if pipeline returns non-zero status
set -o pipefail # return value of last command to exit with non-zero status

database_dir=0
name="sample_table.csv"
existing=0

HELP="""
ASSUMES FASTQ FILES ARE NAMED AS <sample_name>_*_R1_*.fastq(.gz)

Flags:
    -d  :  Absolute path to directory containing fastq(.gz) (can use relative paths, but no guarantee of successful run)
    -n  :  Name of final sample table (with extension) (default "sample_table.csv")
    -e  :  Name of an existing sample table (OPTIONAL) (Will default to creating new sample table, if not provided)
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
	echo "Creating new sample table called ${name}"
 	echo "sample,r1_path,r2_path" > .${name} && echo "sample1,r1_path,r2_path"
else
	name=$(basename $existing)
	echo "Using existing sample table called ${name}"
	cat $existing > .${name}
fi

for file in $database_dir/*.fastq*; do
	sample=$(basename $file | cut -d_ -f 1)
	r1=$(ls $database_dir/$(basename $file | cut -d_ -f1)*R1* | grep ${sample}_)
	r2=$(ls $database_dir/$(basename $file | cut -d_ -f1)*R2* | grep ${sample}_)
	echo ${sample},${r1},${r2} >> .${name} && echo ${sample},${r1},${r2} 
done

uniq .${name} > $name
rm .${name}

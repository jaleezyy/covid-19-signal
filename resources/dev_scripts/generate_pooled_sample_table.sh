#!/bin/bash

#set -e # exit if pipeline returns non-zero status
set -o pipefail # return value of last command to exit with non-zero status
shopt -s extglob

database_dir=0
name="sample_table.csv"
existing=0

HELP="""
ASSUMES FASTQ FILES ARE NAMED AS <sample_name>_*_R1_*.fastq(.gz)

Flags:
    -d  :  Path to directory containing separate directories with pooled fastq(.gz) within them (i.e., <directory>/sample_dir/sample_L001_R1.fastq.gz and <directory>/sample_dir/sample_L001_R2.fastq.gz)
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
 	echo "sample,r1_path,r2_path" > .${name} && echo "sample,r1_path,r2_path"
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

for file in $database_dir/*; do
	if [ -d ""$file"" ]; then # read through directory within directory
		sample=$(basename $(ls $file/*.f?(ast)q* | head -1) | cut -d_ -f 1)
		p1r1=$(ls $file/*L001_R1*.f?(ast)q* | grep /${sample}_)
		p1r2=$(ls $file/*L001_R2*.f?(ast)q* | grep /${sample}_)
		p2r1=$(ls $file/*L002_R1*.f?(ast)q* | grep /${sample}_)
		p2r2=$(ls $file/*L002_R2*.f?(ast)q* | grep /${sample}_)
		echo ${sample},${p1r1},${p1r2} >> .${name} && echo ${sample},${p1r1},${p1r2}
		echo ${sample},${p2r1},${p2r2} >> .${name} && echo ${sample},${p2r1},${p2r2}
	fi
done

samples=()
for file in $database_dir/*.f?(ast)q*; do
	sample=$(basename $file | cut -d_ -f 1)
	if [[ ! " ${samples[@]} " =~ " ${sample} " ]]; then
		p1r1=$(ls $database_dir/$(basename $file | cut -d_ -f1)*L001_R1* | grep /${sample}_)
		p1r2=$(ls $database_dir/$(basename $file | cut -d_ -f1)*L001_R2* | grep /${sample}_)
		p2r1=$(ls $database_dir/$(basename $file | cut -d_ -f1)*L002_R1* | grep /${sample}_)
		p2r2=$(ls $database_dir/$(basename $file | cut -d_ -f1)*L002_R2* | grep /${sample}_)
		samples+=("${sample}")
		echo ${sample},${p1r1},${p1r2} >> .${name} && echo ${sample},${p1r1},${p1r2}
		echo ${sample},${p2r1},${p2r2} >> .${name} && echo ${sample},${p2r1},${p2r2}
	fi
done

uniq .${name} > $name
rm .${name}

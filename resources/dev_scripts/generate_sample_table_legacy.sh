#!/bin/bash

set -e # exit if pipeline returns non-zero status
set -o pipefail # return value of last command to exit with non-zero status
shopt -s extglob


database_dir=0
pooled_dir=0
name="sample_table.csv"
existing=0

HELP="""
ASSUMES FASTQ FILES ARE NAMED AS <sample_name>_*_R1_*.fastq(.gz)

Flags:
    -d  :  Path to directory containing individual sample fastq(.gz) (Absolute paths preferred for consistency, but can use relative paths)
    -p  :  Path to directory containing pooled samples (i.e., 4 FASTQs per sample: L001_R1, L001_R2, L002_R1, L002_R2)
    -n  :  Name of final sample table (with extension) (default: 'sample_table.csv') - will overwrite if file exists
    -e  :  File path to an existing sample table - will append to the end of the provided table

Select one of '-n' (new sample table) or '-e' (existing sample table). 
If neither provided, a new sample table called 'sample_table.csv' will be created (or overwritten) by default. 

Can specify both -d and -p and the generated table will combine both sets of data. It is best to separate the pooled samples from individual samples to avoid confusion. 
"""

while getopts ":d:n:e:p:" option; do
    case "${option}" in
        d) database_dir=$OPTARG;;
        n) name=$OPTARG;;
        e) existing=$OPTARG;;
        p) pooled_dir=$OPTARG;;
    esac
done

if [ $database_dir = 0 ] && [ $pooled_dir = 0 ] ; then
    echo "You must specify a data directory containing fastq(.gz) reads (either -d or -p)."
    echo "$HELP"
    exit 1
fi

if [ $existing = 0 ] ; then
	echo -e "Creating new sample table called ${name}\n"
 	echo "sample,r1_path,r2_path" > ${name} && echo "sample,r1_path,r2_path"
else
	filename=$(basename $existing)
	if [ -f ""$existing"" ] ; then
		echo -e "Using existing sample table called ${filename}\n"
		name=$existing
	else
		echo -e "Sample table does not exist. Check that sample table exists or create a new sample table."
		echo "$HELP"
		exit 1
	fi
fi

if [ ! $database_dir = 0 ]; then
	samples_dir=()
	echo -e "Adding individual samples"
	for file in $database_dir/*.f?(ast)q*; do
		sample=$(basename $file | cut -d_ -f 1)
		if [[ ! " ${samples_dir[@]} " =~ " ${sample} " ]]; then
			r1=$(ls $database_dir/$(basename $file | cut -d_ -f1)*R1* | grep /${sample}_)
			r2=$(ls $database_dir/$(basename $file | cut -d_ -f1)*R2* | grep /${sample}_)
			samples_dir+=("${sample}")
			echo ${sample},${r1},${r2} >> ${name} && echo ${sample},${r1},${r2}
		fi
	done
	echo -e "\n"
fi

if [ ! $pooled_dir = 0 ]; then
	samples_pooled=()
	echo -e "Adding pooled samples"
	for file in $pooled_dir/*.f?(ast)q*; do
		sample=$(basename $file | cut -d_ -f 1)
		if [[ ! " ${samples_pooled[@]} " =~ " ${sample} " ]]; then
			p1r1=$(ls $pooled_dir/$(basename $file | cut -d_ -f1)*_*1_R1* | grep /${sample}_)
			p1r2=$(ls $pooled_dir/$(basename $file | cut -d_ -f1)*_*1_R2* | grep /${sample}_)
			p2r1=$(ls $pooled_dir/$(basename $file | cut -d_ -f1)*_*2_R1* | grep /${sample}_)
			p2r2=$(ls $pooled_dir/$(basename $file | cut -d_ -f1)*_*2_R2* | grep /${sample}_)
			samples_pooled+=("${sample}")
			echo ${sample},${p1r1},${p1r2} >> ${name} && echo ${sample},${p1r1},${p1r2}
			echo ${sample},${p2r1},${p2r2} >> ${name} && echo ${sample},${p2r1},${p2r2}
		fi
	done
fi

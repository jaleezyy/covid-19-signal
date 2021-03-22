#!/bin/bash

set -e # exit if pipeline returns non-zero status
set -o pipefail # return value of last command to exit with non-zero status
shopt -s extglob


database_dir=0
name="sample_table.csv"
existing=0

HELP="""
<<<<<<< HEAD
ASSUMES FASTQ FILES ARE NAMED AS <sample_name>_R{1,2}*.fastq(.gz)
=======
ASSUMES FASTQ FILES ARE NAMED AS <sample_name>_L00#_R{1,2}*.fastq(.gz)
>>>>>>> c77f400b0546bc033b664165898c04c84120f0c0

Flags:
    -d  :  Path to directory containing sample fastq(.gz) files (Absolute paths preferred for consistency, but can use relative paths)
    -n  :  Name or file path for final sample table (with extension) (default: 'sample_table.csv') - will overwrite if file exists
    -e  :  Name or file path for an existing sample table - will append to the end of the provided table

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

echo -e "Adding samples\n"
if [ $existing = 0 ] ; then
	echo -e "Creating new sample table: ${name}\n"
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
	samples_fail=()
	for file in $database_dir/*.f?(ast)q*; do
		sample=$(basename $file | cut -d_ -f 1)
		if [[ ! " ${samples_dir[@]} " =~ " ${sample} " ]] && [[ ! " ${samples_fail[@]} " =~ " ${sample} " ]]; then
<<<<<<< HEAD
			count=$(($(ls $database_dir/${sample}*_R{1,2}*.f?(ast)q* 2>/dev/null | wc -l)/2)) || samples_fail+=("${sample}") # estimate # of files; sample fails if file(s) missing
			for (( i=1; i<=$count; i++ )); do
				r1=$(ls $database_dir/${sample}*_R1* | grep /${sample}_) 
				r2=$(ls $database_dir/${sample}*_R2* | grep /${sample}_)
=======
			count=$(($(ls $database_dir/${sample}*_L00*_R{1,2}*.f?(ast)q* 2>/dev/null | wc -l)/2)) || samples_fail+=("${sample}") # estimate # of files; sample fails if file(s) missing
			for (( i=1; i<=$count; i++ )); do
				r1=$(ls $database_dir/${sample}*_L00${i}_R1* | grep /${sample}_) 
				r2=$(ls $database_dir/${sample}*_L00${i}_R2* | grep /${sample}_)
>>>>>>> c77f400b0546bc033b664165898c04c84120f0c0
				echo ${sample},${r1},${r2} >> ${name} && echo ${sample},${r1},${r2}
			done
			samples_dir+=("${sample}") # sample passed
		fi
	done
	echo -e "\n"
fi

if [ ${#samples_fail[@]} -gt 0 ]; then 
	echo -e "Samples that failed to be added to sample table:"
	for failed in "${samples_fail[@]}"; do
		echo -e $failed
	done
else
	echo -e "Success!"
fi

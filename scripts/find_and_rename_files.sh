#!/usr/bin/env bash

shopt -s extglob

valid='2912419cba616797ed5f2a82c095e6e6'
output=0
raw_data=0
fasta=0
fastq=0
vcf=0
qc=0
test='false'
autoyes='false'

### Functions

check_access () {
	# -r read; -p echo input; -s silent mode; -n no newline
	#while IFS= read -p "$pass_var" -r -s -n 1 letter; do
	#	if [[ $letter == $'\0' ]]; then
	#		break
	#	fi

		# store password
	#	password=password+"$letter"

		# put * in place
	#	pass_var="*"
	#done
	read -p "$pass_var" -s password
	hash=$(echo -n $password | md5sum | cut -d' ' -f1 | tr -d '\n')
	#echo $hash
	if [[ ! $hash == $valid ]]; then
		echo -e "\nAccess denied!"
		exit 1
	else
		echo -e "\nAccess allowed!"
	fi
}

### Parser
HELP="""
Usage:
bash find_and_rename_files.sh -r <raw_illuima_server_data> -f <fasta_dir> -q <fastq_dir> -v <vcf_dir> -s <ncov_qc_summary_dir> [-t] [-y]

This scripts aims to copy over relevant sequencing files from the HHS Illumina server and rename the files accordingly. 

Flags:
	-o  :  Output directory for all sequencing files. If blank, only the draft data will remain under $PWD/draft
	-r  :  (Raw) Sequencing data from HHS Illumina server. Will override -f, -q, -v, -s flags and generate draft directory (see -o)
	-f  :  Directory containing FASTA files
	-q  :  Directory containing FASTQ files
	-v  :  Directory containing VCF files
	-s  :  Directory containing NCOV_TOOLS QC summary files
	-t  :  DEBUG ONLY: Disable all input parameters. Check for access.
	-y  :  Automatically answer "yes" for all checkpoint choices

Credits: Script developed by Jalees A. Nasir, McArthur Lab, @Jaleezyy, 2023
"""

while getopts ":o:r:f:q:v:s:ty" option; do
	case "${option}" in
		o) output=$OPTARG;;
		r) raw_data=$OPTARG;;
		f) fasta=$OPTARG;;
        q) fastq=$OPTARG;;
        v) vcf=$OPTARG;;
		s) qc=$OPTARG;;
        t) test='true';;
		y) autoyes='true';;
	esac
done

### Check if running tests only
if [[ $test == 'true' ]]; then
	output=0
	raw_data=0
	fasta=0
	fastq=0
	vcf=0
	qc=0

	pass_var="Enter password: "
	password=''
	check_access

	echo "dry run"
	read -p "Does the above look correct? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]]  && echo "renaming" || 
	read -p "Do you wish to continue? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || exit 1

	exit 0
fi

### Check parameters (one of either raw_data or directory files)
if [[ $raw_data == 0 ]] && [[ $fasta == 0 || $fastq == 0 || $vcf == 0 || $qc == 0 ]]; then
	echo "$HELP"
	exit 0
fi

if [[ $raw_data != 0 ]]; then
	echo "Generating draft data directory $PWD/draft"
	fasta="$PWD/draft/fasta"
	fastq="$PWD/draft/fastq"
	vcf="$PWD/draft/vcf"
	qc="$PWD/draft/summary_qc"

	if [ -d $fasta ]; then
		rm -r $fasta
		mkdir -p $fasta
	fi
	find $raw_data -name *.consensus.fasta -exec cp {} $fasta \;
	if  [ -d $fastq ]; then
		rm -r $fastq
		mkdir -p $fastq
	fi
	find $raw_data -name *.f?(ast)q?(.gz) -exec cp {} $fastq \;
	if [ -d $vcf ]; then
		rm -r $vcf
		mkdir -p $vcf
	fi
	find $raw_data -name *.ann.vcf -exec cp {} $vcf \;
	if [ -d $qc ]; then
		rm -r $qc
		mkdir -p $qc
	fi
	find $raw_data -name *.tsv -exec cp {} $qc \;

	exit 0
fi

### Password check
pass_var="Enter password: "
password=''
check_access

### Rename files
# FASTA
if [ -d $fasta ]; then
	if [[ $autoyes == 'true' ]]; then
		#rename
	else 
		#rename -n 
		read -p "Does the above look correct? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || #&& rename || 
		read -p "Do you wish to continue? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || exit 1
	fi
fi
# FASTQ
if [ -f $fastq ]; then
	if [[ $autoyes == 'true' ]]; then
		# rename
	else
		# rename -n
		read -p "Does the above look correct? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || #&& rename || 
		read -p "Do you wish to continue? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || exit 1
	fi
fi
# VCF
if [ -d $vcf ]; then
	if [[ $autoyes == 'true' ]]; then
		# rename
	else
		# rename -n
		read -p "Does the above look correct? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || #&& rename || 
		read -p "Do you wish to continue? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || exit 1
	fi
fi
# QC Check
if [ -d $qc ]; then	
	if [[ $autoyes == 'true' ]]; then
		# rename
	else
		read -p "Does the above look correct? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || #&& rename || 
		read -p "Do you wish to continue? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || exit 1
	fi
fi

if [[ $output != 0 ]]; then
	if [ -d $output ]; then
		if [[ $autoyes == 'true' ]]; then
			rm -r $output
			mkdir -p $output
		else
			read -p "Do you wish to overwrite ${output} (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] && 
			rm -r $output && mkdir $output || 
			echo -e "\nPlease update output directory and try again!" && exit 1
		fi
	fi

	for dir in $PWD/draft/*; do 
		cp $dir/* $output
	done
fi

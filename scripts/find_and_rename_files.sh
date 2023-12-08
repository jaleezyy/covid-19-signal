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
maxthreads=$(nproc --all)

### Functions

check_access () {
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
	-y  :  Automatically answer 'yes' for all checkpoint choices

Credits: Script developed by Jalees A. Nasir, McArthur Lab, @Jaleezyy, 2023
Version: 1.1
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

### Password check
#pass_var="Enter password: "
#password=''
#check_access

if [[ $raw_data != 0 ]]; then
	echo "Generating draft data directory $PWD/draft"
	draft="$PWD/draft"
	fasta="$PWD/draft/fasta"
	fastq="$PWD/draft/fastq"
	vcf="$PWD/draft/vcf"
	qc="$PWD/draft/tsv"

	if [ -d $fasta ]; then
		rm -r $fasta
		mkdir -p $fasta
	else
		mkdir -p $fasta 
	fi
	find $raw_data -name *.consensus.fasta -exec cp {} $fasta \;
	if  [ -d $fastq ]; then
		rm -r $fastq
		mkdir -p $fastq
	else
		mkdir -p $fastq
	fi
	find $raw_data -name *.fastq -exec cp {} $fastq \;
	if [ -d $vcf ]; then
		rm -r $vcf
		mkdir -p $vcf
	else
		mkdir -p $vcf
	fi
	find $raw_data -name *.ann.vcf -exec cp {} $vcf \;
	if [ -d $qc ]; then
		rm -r $qc
		mkdir -p $qc
	else
		mkdir -p $qc
	fi
	#find $raw_data -name *.tsv -exec cp {} $qc \;
	find $raw_data -name *_ambiguous_position_report.tsv -exec cp {} $qc \;
	find $raw_data -name *_mixture_report.tsv -exec cp {} $qc \;
	find $raw_data -name *_negative_control_report.tsv -exec cp {} $qc \;
	find $raw_data -name *_summary_qc.tsv -exec cp {} $qc \;
fi

### Rename files
# FASTA
if [ -d $fasta ]; then
	echo "FASTA files:"
	### Remove duplicate samples
	echo "Removing duplicates..."
	for file in $fasta/*; do
		smallest=0
		sample=$(basename $file | cut -d- -f1)
		if [[ $sample == "ExtractPosCtrl" || $sample == "ExtractNegCtrl" ]]; then
			sample=$(basename $file | cut -d- -f1,2)
		fi
		num_files=$(ls $fasta/${sample}-* | wc -l)

		if [ $num_files -eq 2 ]; then
			for f in $fasta/${sample}-*; do
				store=$(basename $f | grep -o "barcode\d*.." | tr -d "barcode")
				if [ $smallest -eq 0 ]; then
					smallest=$store
				elif [ $store -lt $smallest ]; then
					smallest=$store
				fi
			done
		fi
		find $fasta -name ${sample}-barcode${smallest}.* -exec rm {} \;
	done

	if [[ $autoyes == 'true' ]]; then
		cd $fasta
		echo "Renaming samples..."
		rename 's/^1/ON-HRL-23-1/g' *
		rename 's/-barcode\d*\./-v1_/g' *

		echo "Renaming controls..."
		rename 's/^ExtractPosCtrl/Pos/g' *
		rename 's/^ExtractNegCtrl/Neg/g' *
		cd -
		
	else 
		cd $fasta
		echo "Renaming samples..."
		rename -n 's/^1/ON-HRL-23-1/g' *
		read -p "Does the above look correct? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] && echo "Renaming..." && rename 's/^1/ON-HRL-23-1/g' * || 
		read -p "Do you wish to continue? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || exit 1

		rename -n 's/-barcode\d*\./-v1_/g' *
		read -p "Does the above look correct? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] && echo "Renaming..." && rename 's/-barcode\d*\./-v1_/g' * || 
		read -p "Do you wish to continue? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || exit 1

		echo "Renaming controls..."
		rename -n 's/^ExtractPosCtrl/Pos/g' *
		read -p "Does the above look correct? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] && echo "Renaming..." && rename 's/^ExtractPosCtrl/Pos/g' * || 
		read -p "Do you wish to continue? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || exit 1
		
		rename -n 's/^ExtractNegCtrl/Neg/g' *
		read -p "Does the above look correct? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] && echo "Renaming..." && rename 's/^ExtractNegCtrl/Neg/g' * || 
		read -p "Do you wish to continue? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || exit 1
		cd -
	fi
fi
# FASTQ
if [ -d $fastq ]; then
	echo "FASTQ files:"
	### Remove duplicate samples
	echo "Removing duplicates..."
	for file in $fastq/*; do
		smallest=0
		sample=$(basename $file | cut -d- -f1)
		if [[ $sample == "ExtractPosCtrl" || $sample == "ExtractNegCtrl" ]]; then
			sample=$(basename $file | cut -d- -f1,2)
		fi
		num_files=$(ls $fastq/${sample}-* | wc -l)

		if [ $num_files -eq 2 ]; then
			for f in $fastq/${sample}-*; do
				store=$(basename $f | cut -d_ -f1 | grep -o "barcode\d*.." | tr -d "barcode")
				if [ ${smallest} -eq 0 ]; then
					smallest=$store
				elif [ $store -lt $smallest ]; then
					smallest=$store
				fi
			done
		fi
		find $fastq -name ${sample}-barcode${smallest}_barcode${smallest}.* -exec rm {} \;
	done

	if [[ $autoyes == 'true' ]]; then
		cd $fastq
		echo "Renaming samples..."
		rename 's/^1/ON-HRL-23-1/g' *
		rename 's/-barcode\d*_barcode\d*\.fastq/-v1\.fq/g' *

		echo "Renaming controls..."
		rename 's/^ExtractPosCtrl/Pos/g' *
		rename 's/^ExtractNegCtrl/Neg/g' *
		cd -

		echo "Zipping..."
		gzip $fastq/*
	else
		cd $fastq
		echo "Renaming samples..."
		rename -n 's/^1/ON-HRL-23-1/g' *
		read -p "Does the above look correct? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] && echo "Renaming..." && rename 's/^1/ON-HRL-23-1/g' * || 
		read -p "Do you wish to continue? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || exit 1

		rename -n 's/-barcode\d*_barcode\d*\.fastq/-v1\.fq/g' *
		read -p "Does the above look correct? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] && echo "Renaming..." && rename 's/-barcode\d*_barcode\d*\.fastq/-v1\.fq/g' * || 
		read -p "Do you wish to continue? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || exit 1
		
		echo "Renaming controls..."
		rename -n 's/^ExtractPosCtrl/Pos/g' *
		read -p "Does the above look correct? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] && echo "Renaming..." && rename 's/^ExtractPosCtrl/Pos/g' * || 
		read -p "Do you wish to continue? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || exit 1
		
		rename -n 's/^ExtractNegCtrl/Neg/g' *
		read -p "Does the above look correct? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] && echo "Renaming..." && rename 's/^ExtractNegCtrl/Neg/g' * || 
		read -p "Do you wish to continue? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || exit 1
		cd -

		echo "Zipping..."
		#gzip $fastq/*
		parallel -j ${maxthreads} --gnu pigz -p ${maxthreads} ::: $fastq/*
	fi
fi
# VCF
if [ -d $vcf ]; then
	echo "VCF files:"
	### Remove duplicate samples
	echo "Removing duplicates..."
	for file in $vcf/*; do
		smallest=0
		sample=$(basename $file | cut -d- -f1)
		if [[ $sample == "ExtractPosCtrl" || $sample == "ExtractNegCtrl" ]]; then
			sample=$(basename $file | cut -d- -f1,2)
		fi
		num_files=$(ls $vcf/${sample}-* | wc -l)

		if [ $num_files -eq 2 ]; then
			for f in $vcf/${sample}-*; do
				store=$(basename $f | grep -o "barcode\d*.." | tr -d "barcode")
				if [ $smallest -eq 0 ]; then
					smallest=$store
				elif [ $store -lt $smallest ]; then
					smallest=$store
				fi
			done
		fi
		find $vcf -name ${sample}-barcode${smallest}.* -exec rm {} \;
	done

	if [[ $autoyes == 'true' ]]; then
		cd $vcf
		echo "Renaming samples..."
		rename 's/^1/ON-HRL-23-1/g' *
		rename 's/-barcode\d*\.ann/-v1/g' *

		echo "Renaming controls..."
		rename 's/^ExtractPosCtrl/Pos/g' *
		rename 's/^ExtractNegCtrl/Neg/g' *
		cd -
	else
		cd $vcf
		echo "Renaming samples..."
		rename -n 's/^1/ON-HRL-23-1/g' *
		read -p "Does the above look correct? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] && echo "Renaming..." && rename 's/^1/ON-HRL-23-1/g' * || 
		read -p "Do you wish to continue? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || exit 1

		rename -n 's/-barcode\d*\.ann/-v1/g' *
		read -p "Does the above look correct? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] && echo "Renaming..." && rename 's/-barcode\d*\.ann/-v1/g' * || 
		read -p "Do you wish to continue? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || exit 1

		echo "Renaming controls..."
		rename -n 's/^ExtractPosCtrl/Pos/g' *
		read -p "Does the above look correct? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] && echo "Renaming..." && rename 's/^ExtractPosCtrl/Pos/g' * || 
		read -p "Do you wish to continue? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || exit 1
		
		rename -n 's/^ExtractNegCtrl/Neg/g' *
		read -p "Does the above look correct? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] && echo "Renaming..." && rename 's/^ExtractNegCtrl/Neg/g' * || 
		read -p "Do you wish to continue? (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] || exit 1
		cd -
	fi
fi
# QC Check
if [ -d $qc ]; then
	echo "NCOV_QC files:"
	last=0
	for file in $qc/*; do
		sample=$(basename $file | cut -d_ -f1,2)
		if [ $sample != $last ]; then
			if [ ! -f $qc/"${sample}_ambiguous_position_report.tsv" ]; then
				echo -e "Missing ${sample}_ambiguous_position_report.tsv!"
			fi
			if [ ! -f $qc/"${sample}_mixture_report.tsv" ]; then
				echo -e "Missing ${sample}_mixture_report.tsv!"
			fi
			if [ ! -f $qc/"${sample}_negative_control_report.tsv" ]; then
				echo -e "Missing ${sample}_negative_control_report.tsv!"
			fi
			if [ ! -f $qc/"${sample}_summary_qc.tsv" ]; then
				echo -e "Missing ${sample}_summary_qc.tsv!"
			fi
			last=$sample
		fi
	done
fi

if [[ $output != 0 ]]; then
	if [ -d $output ]; then
		if [[ $autoyes == 'true' ]]; then
			rm -r $output
			mkdir -p $output
		else
			read -p "Do you wish to overwrite ${output} (Y/N): " confirm && [[ $confirm == [yY] || $confirm == [yY][eE][sS] ]] && 
			rm -r $output && mkdir -p $output || exit 1
		fi
	else
		mkdir -p $output
	fi

	if [ -d $fasta ]; then
		for f in $fasta; do 
			cp $f/* $output
		done
	fi
	if [ -d $fastq ]; then
		for f in $fastq; do 
			cp $f/* $output
		done
	fi
	if [ -d $vcf ]; then
		for f in $vcf; do 
			cp $f/* $output
		done
	fi
	if [ -d $qc ]; then
		for f in $qc; do 
			cp $f/* $output
		done
	fi

	# clean up
	if [[ $raw_data != 0 ]]; then
		rm -r $draft
	fi
	
	echo "All data found under ${output}"
else
	echo -e "All data found where inputs are. If input was -r, then check to see if ${draft} remains"
fi

exit 0

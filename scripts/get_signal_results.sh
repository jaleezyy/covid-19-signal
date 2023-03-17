#!/bin/env bash

shopt -s extglob

source=0
destination=0
move='false'
copy='false'

HELP="""
Usage:
bash get_signal_results.sh -s <SIGNAL_results_dir> -d <destination_dir> [-m] [-c]

This scripts aims to copy (rsync by default) or move (mv) select output from SIGNAL 'all', 'postprocess', and 'ncov_tools'.

The following files will be transferred over to the specified destination directory (if found):
SIGNAL 'all' & 'postprocess':
-> signal-results/<sample>/<sample>_sample.txt
-> signal-results/<sample>/core/<sample>.consensus.fa
-> signal-results/<sample>/core/<sample>_ivar_variants.tsv
-> signal-results/<sample>/freebayes/<sample>.consensus.fasta
-> signal-results/<sample>/freebayes/<sample>.variants.norm.vcf

SIGNAL 'ncov_tools':
-> ncov_tools-results/qc_annotation/<sample>.ann.vcf
-> ncov-tools-results/qc_reports/<run_name>_ambiguous_position_report.tsv
-> ncov-tools-results/qc_reports/<run_name>_mixture_report.tsv
-> ncov-tools-results/qc_reports/<run_name>_ncov_watch_variants.tsv
-> ncov-tools-results/qc_reports/<run_name>_negative_control_report.tsv
-> ncov-tools-results/qc_reports/<run_name>_summary_qc.tsv

Flags:
	-s  :  SIGNAL results directory
	-d  :  Directory where summary will be outputted
	-m  :  Invoke 'mv' move command instead of 'rsync' copying of results. Optional
	-c  :  Invoke 'cp' copy command instead of 'rsync' copying of results. Optional
"""

while getopts ":s:d:mc" option; do
	case "${option}" in
		s) source=$OPTARG;;
		d) destination=$OPTARG;;
		m) move='true';;
		c) copy='true';;
	esac
done


if [ $source = 0 ] || [ $destination = 0 ] ; then
	echo "You must specify both source and destination locations."
	echo "$HELP"
	exit 1
fi

if [ ! -d $destination ]; then
	echo "Invalid destination directory!"
	exit 1
fi

if [ ! -f $source/summary.html ] && [ ! -f $source/summary.zip ]; then
	echo "Invalid SIGNAL directory! Make sure you've run SIGNAL 'all' and 'postprocess'!"
	exit 1
else
	run_name=$(basename $source)
	final_dir=${destination}/${run_name}
	mkdir -p $final_dir/signal-results
fi
	
if [ ${move} = true ] && [ ${copy} = true ]; then
	echo -e "Only pick one of '-m' or '-c' depending on whether you wish to move or copy files, respectively"
	exit
elif [ ${move} = true ] && [ ${copy} = false ]; then
	cmd='mv'
elif [ ${move} = false ] && [ ${copy} = true ]; then
	cmd='cp'
else
	cmd='rsync -avh'
	# rsync -avh
fi

echo -e "We will use ${cmd} for your files!"

### SIGNAL results_dir
for file in $source/*; do 
	if [ -d $file ]; then # results_dir/sample
		sample=$(basename $file) # sample name, within contain our files
		sample_dest=${final_dir}/'signal-results'/${sample}
		if [[ ! $sample == 'ncov-tools-results' ]]; then
			mkdir -p $sample_dest
		fi
		for d in $file/*; do
			name=$(basename $d)
			if [ -d $d ] && [[ $name == 'core' ]]; then
				mkdir -p $sample_dest/core
				$cmd ${d}/${sample}.consensus.fa $sample_dest/core/${sample}.consensus.fa
				$cmd ${d}/${sample}_ivar_variants.tsv $sample_dest/core/${sample}_ivar_variants.tsv
			elif [ -d $d ] && [[ $name == 'freebayes' ]]; then
				mkdir -p $sample_dest/freebayes
				$cmd ${d}/${sample}.consensus.fasta $sample_dest/freebayes/${sample}.consensus.fasta 
				$cmd ${d}/${sample}.variants.norm.vcf $sample_dest/freebayes/${sample}.variants.norm.vcf
			elif [ -f $d ] && [[ $name =~ '_sample.txt' ]]; then
				$cmd ${d} $sample_dest/$(basename $d)
			else
				continue
			fi
		done
	fi
done

echo "Files from SIGNAL transferred!"

### NCOV-TOOLS
if [ ! -d $source/ncov-tools-results ]; then
	echo "No ncov-tools-results directory found!"
else
	ncov_dest=${final_dir}/ncov-tools-results
	mkdir -p $ncov_dest/qc_{annotation,reports}
	$cmd $source/ncov-tools-results/qc_reports/* $ncov_dest/qc_reports
	$cmd $source/ncov-tools-results/qc_annotation/*.ann.vcf  $ncov_dest/qc_annotation
	
	echo "Files from ncov-tools transferred!"
fi
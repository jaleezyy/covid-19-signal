#!/bin/bash

#set -e # exit if pipeline returns non-zero status
#set -o pipefail # return value of last command to exit with non-zero status

database_dir=0
table=0
name="signal-results"

HELP="""
Generate directory (in the current working directory) with symbolic links to SIGNAL output. 
Flags:
    -d  :  Absolute path to directory containing output from SIGNAL pipeline.
    -n  :  Name or filepath of output directory (default: 'signal-results') - directory must not already exist - should match in config.yaml for ncov-tools.
    -t  :  Name or filepath to sample table used in SIGNAL
"""

while getopts ":d:n:t:" option; do
    case "${option}" in
        d) database_dir=$OPTARG;;
        n) name=$OPTARG;;
		t) table=$OPTARG;;
    esac
done

if [ $database_dir = 0 ] ; then
	echo "You must specify a data directory containing SIGNAL output."
	echo "$HELP"
	exit 1
fi

if [ $table = 0 ] ; then
	echo "You must specify a sample table used in SIGNAL."
	echo "$HELP"
	exit 1
fi

# take sample table and derive sample IDs
samples=()
while read line; do
	id=$(echo $line | cut -d, -f1)
	samples+=("${id}")
done < $table

if [ ""${samples[0]}"" != "sample" ]; then
	echo "Error in sample table. Please double-check and try and again."
	exit 1
fi

# create output directory
if [ -d ${name} ]; then
	echo "Output directory already exists."
	echo "$HELP"
	exit 1
else
	mkdir $name
fi

# examine directories within SIGNAL results, should line up with sample IDs
for dir in ${database_dir}/*; do
	if [ ! -d "${dir}" ]; then 
	# ignore non-directories
		continue
	elif [ "${dir}" == "${database_dir}/.snakemake" ] ; then 
	# ignore snakemake file
		continue
	elif [[ " ${samples[@]} " =~ " $(basename $dir) " ]] && [ -d "${dir}/core" ]; then 
	# directory should match sample and contain a 'core' directory within it
		sample=$(basename ${dir})
		echo "${dir}/core/${sample}.consensus.fa"
		ln -s ${dir}/core/${sample}.consensus.fa ${name}/${sample}.consensus.fa 2>/dev/null || echo "Not found"
		echo "${dir}/core/${sample}_viral_reference.mapping.primertrimmed.sorted.bam"
		ln -s ${dir}/core/${sample}_viral_reference.mapping.primertrimmed.sorted.bam ${name}/${sample}.sorted.bam 2>/dev/null || echo "Not found"
		echo "${dir}/core/${sample}_ivar_variants.tsv"
		ln -s ${dir}/core/${sample}_ivar_variants.tsv ${name}/${sample}.variants.tsv 2>/dev/null || echo "Not found"
	fi
done

# echo -e "\nIgnore any errors from cleanup.\n"

# rm ${name}/summary*
# rm ${name}/config.yaml*
# rm ${name}/*.csv{.consensus.fa,.sorted.bam,.variants.tsv}
# rm ${name}/*.tsv{.consensus.fa,.sorted.bam,.variants.tsv}
# rm ${name}/bioproject*
# rm ${name}/results*
# rm ${name}/ncov-tools-results.{consensus.fa,sorted.bam,variants.tsv}

echo -e "\nDon't forget to update the config.yaml file prior to running ncov-tools."

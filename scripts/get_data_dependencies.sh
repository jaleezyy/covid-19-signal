#!/usr/bin/env bash

#set -e # exit if pipeline returns non-zero status
#set -o pipefail # return value of last command to exit with non-zero status
source ~/.bashrc

if [ -z "$CONDA_EXE" ] ; then
  echo "This script requires conda but it looks like you don't have conda activated" >&2
  exit 1
fi

database_dir=0
accession="MN908947.3"

HELP="""
Flags:
	-d  :  Directory to configure database within (~10GB)
	-a  :  Accession to use as viral reference (default=MN908947.3)
"""

while getopts ":d:a:" option; do
	case "${option}" in
		d) database_dir=$OPTARG;;
		a) accession=$OPTARG;;
	esac
done

if [ $database_dir = 0 ] ; then
	echo "You must specify a data directory to install data dependencies."
	echo "$HELP"
	exit 1
fi

echo -e "Warning: \n - final databases require ~10GB of storage\n - building databases temporarily requires a peak of ~35GB of storage and ~4GB of memory \n - script takes up to ~1.5 hours (system depending)"

# make database dir and get abspath to it
if [ ! -d $database_dir ]; then mkdir -p $database_dir; fi
database_dir=$(realpath $database_dir)

# use curl to grab "simple data dependencies"
if [ ! -f $database_dir/'NexteraPE-PE.fa' ]; then
	curl -s "https://raw.githubusercontent.com/timflutre/trimmomatic/3694641a92d4dd9311267fed85b05c7a11141e7c/adapters/NexteraPE-PE.fa" > $database_dir/NexteraPE-PE.fa
fi
if [ ! -f $database_dir/${accession}.gbk ]; then
	curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${accession}&rettype=gb&retmode=txt" > $database_dir/$accession.gbk
fi
if [ ! -f $database_dir/${accession}.gff3 ]; then
	curl -s "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=${accession}" > $database_dir/$accession.gff3
fi
if [ ! -f $database_dir/${accession}.fasta ]; then
	curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${accession}&rettype=fasta&retmode=txt" > $database_dir/$accession.fasta
fi

# install and activate env for kraken/bwa to build their databases/index
CONDA_BASE=$($CONDA_EXE info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda create -n data_dependencies -c conda-forge -c bioconda -y kraken2=2.1.1 bwa
conda activate data_dependencies

# get the GRCh38 human genome
# as per https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use
if [ ! -f $database_dir/"GRC38_no_alt_analysis_set.fna" ]; then
	curl -s "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz" > $database_dir/GRC38_no_alt_analysis_set.fna.gz
	gunzip $database_dir/GRC38_no_alt_analysis_set.fna.gz
fi

# create composite reference of human and virus for competitive bwt mapping 
# based host removal
if [ ! -f $database_dir/'composite_human_viral_reference.fna' ]; then
cat $database_dir/GRC38_no_alt_analysis_set.fna $database_dir/$accession.fasta > $database_dir/composite_human_viral_reference.fna
fi
for file in $database_dir/composite_human_viral_reference.fna.{amb,ann,bwt,pac,sa}; do
	if [ ! -f $file ]; then
		bwa index $database_dir/composite_human_viral_reference.fna
		break
	else
		continue
	fi
done

# get kraken2 viral db
mkdir -p $database_dir/Kraken2/db
for file in $database_dir/Kraken2/db/{hash,opts,taxo}.k2d; do
	if [ ! -f $file ]; then
		curl -s "https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20210517.tar.gz" > $database_dir/Kraken2/db/k2_viral_20210517.tar.gz
		cd $database_dir/Kraken2/db
		tar xvf k2_viral_20210517.tar.gz
		break
	else
		continue
	fi
done

# create blank fasta for 'phylo_include_seqs'
touch $database_dir/blank.fasta

#!/bin/bash

set -e # exit if pipeline returns non-zero status
set -o pipefail # return value of last command to exit with non-zero status

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

echo -e "Warning: \n - final databases require ~10GB of storage\n - building databases temporarily requires a peak of ~35GB of storage and ~2GB of memory \n - script takes up to ~1.5 hours (system depending)"

# make database dir and get abspath to it
mkdir -p $database_dir
database_dir=$(realpath $database_dir)

# use curl to grab "simple data dependencies"
curl -s "https://raw.githubusercontent.com/timflutre/trimmomatic/3694641a92d4dd9311267fed85b05c7a11141e7c/adapters/NexteraPE-PE.fa" > $database_dir/NexteraPE-PE.fa
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${accession}&rettype=gb&retmode=txt" > $database_dir/$accession.gbk
curl -s "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=${accession}" > $database_dir/$accession.gff3
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${accession}&rettype=fasta&retmode=txt" > $database_dir/$accession.fasta

# install and activate env for lmat/kraken to build their databases
conda create -n data_dependencies -c conda-forge -c bioconda -y kraken2=2.0.8
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate data_dependencies

# get kraken2, and clean db after building
kraken2-build --download-taxonomy --db $database_dir/Kraken2/db --threads 10 --use-ftp
kraken2-build --download-library viral --db $database_dir/Kraken2/db --threads 10 --use-ftp
kraken2-build --build --threads 10 --db $database_dir/Kraken2/db
kraken2-build --clean --threads 10 --db $database_dir/Kraken2/db

# get the GRCh38 BWA index human genome
curl -s "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.bwa_index.tar.gz" > $database_dir/GRCh38_bwa.tar.gz
tar xvf $database_dir/GRCh38_bwa.tar.gz -C $database_dir

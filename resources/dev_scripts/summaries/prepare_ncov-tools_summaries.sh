### usage ./prepare_ncov-tools_summaries.sh /workspace/nasirja/re-infection-round2/pipeline

accession=MN908947.3
amplicon_bed=nCov-2019.V3.regions.bed
cores=20

## git clone ncov-tools
git clone https://github.com/jts/ncov-tools.git

cd ncov-tools

## install environment using
conda env create -f environment.yml

## acticate environment
conda acticate ncov-qc

## install samtools
conda install -c conda-forge -c bioconda hisat2=2.1.0
conda install -c conda-forge -c bioconda samtools=1.7

mkdir all_samples
mkdir resources
mkdir logs
cp etc/${amplicon_bed} resources


curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=${accession}&rettype=fasta&retmode=txt" > resources/${accession}.fasta

## write config file
echo """
# path to the top-level directory containing the analysis results
data_root: all_samples

# path to the file containing the amplicon regions (not the primer sites, the actual amplicons)
amplicon_bed: resources/${amplicon_bed}

# path to the nCov reference genome
reference_genome: resources/${accession}.fasta

# the naming convention for the bam files
# this can use the variables {data_root} (as above) and {sample}
# As per the example above, this will expand to run_200430/sampleA.sorted.bam for sampleA
bam_pattern: '{data_root}/{sample}.sorted.bam'

# the naming convention for the consensus sequences
consensus_pattern: '{data_root}/{sample}.consensus.fasta'

# when building a tree of the consensus genomes you can optionally include other sequences
# in the tree by providing a fasta file here
#tree_include_consensus: some_genomes_from_gisaid.fasta

# some plots can by annotated with external metadata. this
# file contains the metadata in simple tab-delimited format
# one column must be 'sample'
#metadata: metadata.tsv

""" > config.yaml

## write run script
echo """
# Build the sequencing QC plots (coverage, allele frequencies)
snakemake -p -s qc/Snakefile all_qc_sequencing --keep-incomplete --keep-going --cores ${cores} > logs/sequencing.log 2>&1

# Build the analysis QC plots (tree with annotated mutations)
snakemake -p -s qc/Snakefile all_qc_analysis --keep-incomplete --keep-going --cores ${cores} > logs/analysis.log 2>&1

# Build the plots that require a metadata file with a 'ct' column
#snakemake -s qc/Snakefile all_qc_by_ct --keep-incomplete --keep-going --cores ${cores} > logs/qc_by_ct.log 2>&1

""" > run.sh

chmod +x run.sh

## build index for reference 
hisat2-build resources/MN908947.3.fasta resources/genome

## copy consensus fastas and host_removed r1 and r2

base_directory_for_all_samples=$1

for f in $base_directory_for_all_samples/*
do
  sample=$(basename $f)
  if [ -f ${base_directory_for_all_samples}/${sample}/consensus/virus.consensus.fa ]; 
  then
  echo "${sample}"
  mkdir ${sample}
  ## copy consensus fastas and rename
  cp ${base_directory_for_all_samples}/${sample}/consensus/virus.consensus.fa  ${sample}/${sample}.consensus.fa
  cat ${sample}/${sample}.consensus.fa | perl -ane "if(/\>/){print \">${sample}\n\"; } else { print; }" > ${sample}/${sample}.consensus.fasta
  ## copy host removed r1 and r2
  cp ${base_directory_for_all_samples}/${sample}/host_removed/R*.fastq.gz ${sample}/

  ## align reads to reference
	hisat2 --threads ${cores} \
	-x resources/genome \
	-1 ${sample}/R1.fastq.gz \
	-2 ${sample}/R2.fastq.gz \
	--summary-file ${sample}/hisat2_summary.txt \
	-S ${sample}/output.sam \

  ## Convert SAM to BAM, with sort
	samtools view --threads ${cores} -b \
	${sample}/output.sam \
	| samtools sort \
	> ${sample}/${sample}.sorted.bam

fi

done

## copy ${sample}.sorted.bam and ${sample}.consensus.fasta in all_samples
cp */*.sorted.bam */*.consensus.fasta all_samples

## run the run.sh script to create summaries
./run.sh


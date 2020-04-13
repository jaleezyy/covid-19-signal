# example: ./remove_host_sequences.sh Wuhan1-LN
SAMPLE=$1

# build db
hisat2-build -p 25 MN908947_3.fasta sars-cov-2

# align reads to reference
# bowtie2 mapping against host sequence database, keep both mapped and unmapped reads (paired-end reads)
# please note that ${SAMPLE}-R1_paired.fastq.gz and ${SAMPLE}-R2_paired.fastq.gz should be from after merging, trimming, and qc
hisat2 -x sars-cov-2 --threads 20 -1 ${SAMPLE}-R1_paired.fastq.gz -2 ${SAMPLE}-R2_paired.fastq.gz  -S  ${SAMPLE}_mapped_and_unmapped.sam

# convert file .sam to .bam
samtools view --threads 20 -bS ${SAMPLE}_mapped_and_unmapped.sam > ${SAMPLE}_mapped_and_unmapped.bam

# SAMtools SAM-flag filter: get mapped pairs (both ends mapped)
samtools view --threads 20 -b -f 3 -F 4 ${SAMPLE}_mapped_and_unmapped.bam > ${SAMPLE}_both_ends_mapped.bam

# sort bam file by read name (-n) to have paired reads next to each other as required by bedtools
samtools sort --threads 20 -T tmp -n ${SAMPLE}_both_ends_mapped.bam -o ${SAMPLE}_both_ends_mapped_sorted.bam

# write reads for paired reads i.e R1 and R2
bedtools bamtofastq -i ${SAMPLE}_both_ends_mapped_sorted.bam -fq ${SAMPLE}_host_removed_R1.fastq -fq2 ${SAMPLE}_host_removed_R2.fastq

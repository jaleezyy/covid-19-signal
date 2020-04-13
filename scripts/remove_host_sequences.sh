# example: ./remove_host_sequences.sh Wuhan1-LN
SAMPLE=$1
threads=20
# build db
hisat2-build -p ${threads} MN908947_3.fasta sars-cov-2

# align reads to reference
# bowtie2 mapping against host sequence database, keep both mapped and unmapped reads (paired-end reads)
# please note that ${SAMPLE}-R1_paired.fastq.gz and ${SAMPLE}-R2_paired.fastq.gz should be from after merging, trimming, and qc
hisat2 -x sars-cov-2 --threads ${threads} -1 ${SAMPLE}-R1_paired.fastq.gz -2 ${SAMPLE}-R2_paired.fastq.gz  -S  ${SAMPLE}_mapped_and_unmapped.sam

# convert file .sam to .bam
samtools view --threads ${threads} -bS ${SAMPLE}_mapped_and_unmapped.sam > ${SAMPLE}_mapped_and_unmapped.bam

# SAMtools SAM-flag filter: get mapped pairs (both ends mapped)
samtools view --threads ${threads} -b -f 3 -F 4 ${SAMPLE}_mapped_and_unmapped.bam > ${SAMPLE}_both_ends_mapped.bam

# sort bam file by read name (-n) to have paired reads next to each other as required by bedtools
samtools sort --threads ${threads} -T tmp -n ${SAMPLE}_both_ends_mapped.bam -o ${SAMPLE}_both_ends_mapped_sorted.bam

# write reads for paired reads i.e R1 and R2
bedtools bamtofastq -i ${SAMPLE}_both_ends_mapped_sorted.bam -fq ${SAMPLE}_host_removed_R1.fastq -fq2 ${SAMPLE}_host_removed_R2.fastq

# assembly
unicycler \
-1 ${SAMPLE}_host_removed_R1.fastq \
-2 ${SAMPLE}_host_removed_R2.fastq \
--threads ${threads} \
--mode conservative \
--verbosity 2 \
-o ${SAMPLE}

# compress fastq
gzip *_host_removed_R1.fastq
gzip *_host_removed_R2.fastq

# rename  assemblies and graph
cp ${SAMPLE}/assembly.fasta ${SAMPLE}-unicycler-round2.fasta
cp ${SAMPLE}/assembly.gfa ${SAMPLE}-unicycler-round2.gfa

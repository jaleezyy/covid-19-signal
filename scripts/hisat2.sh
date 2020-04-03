
# build reference
hisat2-build assembly.fasta genome
# align
hisat2 -x genome --threads 20 -1 *R1_paired.fastq.gz -2 *R2_paired.fastq.gz --summary-file -S output.sam
# sam to bam
samtools view --threads 20 -b  output.sam > output.bam


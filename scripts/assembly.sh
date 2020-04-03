# DECOMPRESS
gunzip *.fastq.gz

# CONCAT
cat MT-swab-Iran-Liverpool-pool1-baited_S15_L001_R1_001.fastq MT-swab-Iran-Liverpool-pool2-baited_S19_L001_R1_001.fastq > Iran1-LA-R1.fastq
cat MT-swab-Iran-Liverpool-pool1-baited_S15_L001_R2_001.fastq MT-swab-Iran-Liverpool-pool2-baited_S19_L001_R2_001.fastq > Iran1-LA-R2.fastq

# SORT
cat Iran1-LA-R1.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > Iran1-LA-R1_sorted.fastq
cat Iran1-LA-R2.fastq | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > Iran1-LA-R2_sorted.fastq

# COMPRESS
gzip *.fastq

# QC
fastqc *_sorted.fastq.gz
mkdir fastqc_on_concat_reads
mv *.html *.zip fastqc_on_concat_reads

# REMOVE PRIMERS
cutadapt \
-a file:/home/raphenar/sars-cov-2/wuhan_primers_28.01.20_trim_RC.fa \
-A file:/home/raphenar/sars-cov-2/wuhan_primers_28.01.20_trim_FW.fa \
-o Iran1-LA-R1_paired-primers-removed.fastq \
-p Iran1-LA-R2_paired-primers-removed.fastq \
Iran1-LA-R1_sorted.fastq.gz  Iran1-LA-R2_sorted.fastq.gz \
-j 25 > cutadapt.log 2>&1

# compress
gzip  *-primers-removed.fastq

# TRIM
java -mx20G -jar /workspace/tsangkk2/ORF/scripts/Trimmomatic/Trimmomatic-0.36/trimmomatic-0.36.jar PE \
-threads 20 \
Iran1-LA-R1_paired-primers-removed.fastq.gz  Iran1-LA-R2_paired-primers-removed.fastq.gz \
Iran1-LA-R1_paired.fastq Iran1-LA-R1_unpaired.fastq Iran1-LA-R2_paired.fastq Iran1-LA-R2_unpaired.fastq \
ILLUMINACLIP:/workspace/tsangkk2/ORF/scripts/Trimmomatic/Trimmomatic-0.36/adapters/NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:4:20 > trim.log 2>&1

# compress
gzip *.fastq

# QC
fastqc *paired.fastq.gz
mkdir fastqc_on_trimmed_reads
mv *.html *.zip fastqc_on_trimmed_reads

# ASSEMBLY
unicycler \
-1 Iran1-LA-R1_paired.fastq.gz \
-2 Iran1-LA-R2_paired.fastq.gz \
--threads 20 \
--mode conservative \
--verbosity 2 \
-o Iran1-LA


## Pipeline steps and running example:

**Note:** Slightly out of date, will udpate soon!

For each step in the pipeline, we give a high-level description,
followed by the shell commands executed (obtained from `snakemake -p`),
in a running example with the following input files:
```
        fastq_inputs/MT-swab-Iran-Liverpool-pool1_S3_L001_R1_001.fastq.gz
        fastq_inputs/MT-swab-Iran-Liverpool-pool1_S3_L001_R2_001.fastq.gz
        fastq_inputs/MT-swab-Iran-Liverpool-pool2_S7_L001_R1_001.fastq.gz
        fastq_inputs/MT-swab-Iran-Liverpool-pool2_S7_L001_R2_001.fastq.gz
```
1. Concatenate and sort
```
        zcat fastq_inputs/MT-swab-Iran-Liverpool-pool1_S3_L001_R1_001.fastq.gz \
             fastq_inputs/MT-swab-Iran-Liverpool-pool2_S7_L001_R1_001.fastq.gz \
        | paste - - - - \
        | sort -k1,1 -t " " \
        | tr "\t" "\n" \
        | gzip \
        > fastq_sorted/Iran1_R1.fastq.gz

        zcat fastq_inputs/MT-swab-Iran-Liverpool-pool1_S3_L001_R2_001.fastq.gz \
             fastq_inputs/MT-swab-Iran-Liverpool-pool2_S7_L001_R2_001.fastq.gz \
        | paste - - - - \
        | sort -k1,1 -t " " \
        | tr "\t" "\n" \
        | gzip \
        > fastq_sorted/Iran1_R2.fastq.gz
```
2. FASTQC analysis of sorted files
```
        fastqc fastq_sorted/Iran1_R1.fastq.gz
        fastqc fastq_sorted/Iran1_R2.fastq.gz
```
3. Remove amplification primer sequences using cutadapt
```
        cutadapt -j 25 \
          -a file:/home/kmsmith/data/wuhan_primers_28.01.20_trim_RC.fa \
          -A file:/home/kmsmith/data/wuhan_primers_28.01.20_trim_FW.fa \
          -o fastq_primers_removed/Iran1_R1.fastq.gz \    # R1 output file
          -p fastq_primers_removed/Iran1_R2.fastq.gz \    # R2 output file
           fastq_sorted/Iran1_R1.fastq.gz \                # R1 input file
           fastq_sorted/Iran1_R2.fastq.gz                  # R2 input file
```
- **Question:** The Wuhan primer files are currently hardcoded here. Is this okay or should we do something different?

4. Run Trimmomatic, to remove Illumina adapter sequences, and trim/remove low quality sequences
```
        trimmomatic PE -threads 20 \
        fastq_primers_removed/Iran1_R1.fastq.gz \
        fastq_primers_removed/Iran1_R2.fastq.gz \
        fastq_trimmed/Iran1_R1_paired.fastq.gz \
        fastq_trimmed/Iran1_R1_unpaired.fastq.gz \
        fastq_trimmed/Iran1_R2_paired.fastq.gz \
        fastq_trimmed/Iran1_R2_unpaired.fastq.gz \
        ILLUMINACLIP:/home/kmsmith/data/NexteraPE-PE.fa:2:30:10 \
        SLIDINGWINDOW:4:20
```
- **Question:** the NexteraPE file is currently hardcoded here. Is this okay or should we do something different?

5. FASTQC analysis of trimmed files
```
        fastqc fastq_trimmed/Iran1_R1_unpaired.fastq.gz
        fastqc fastq_trimmed/Iran1_R1_paired.fastq.gz
        fastqc fastq_trimmed/Iran1_R2_unpaired.fastq.gz
        fastqc fastq_trimmed/Iran1_R2_paired.fastq.gz
```
6. Use breseq to analyze mutation list relative to MN908947.3
```
        breseq --reference /home/kmsmith/data/MN908947_3.gbk \
          --num-processors 80 --polymorphism-prediction --brief-html-output \
          --output Iran1_breseq \
          fastq_trimmed/Iran1_R1_paired.fastq.gz \
          fastq_trimmed/Iran1_R2_paired.fastq.gz
```
- **Question:** MN908947.3 is currently hardcoded here. Is this okay or should we do something different?

7. Use kraken2 to derive percentage of reads derived from SARS-Cov-2 RNA
```
        mkdir -p Iran1_kraken2
        cd Iran1_kraken2
	
        kraken2 --db /home/kmsmith/data/Kraken2/db \
          --threads 20 --quick --unclassified-out unclassified-sequences# --classified-out classified-sequences# \
          --output kraken2.out \
          --paired --gzip-compressed \
          ../fastq_trimmed/Iran1_R1_paired.fastq.gz \
          ../fastq_trimmed/Iran1_R2_paired.fastq.gz \
          --report report
```
- **Note/question:** I'm using the Kraken2 DB from `galaxylab:/home/raphenar/sars-cov-2/round1/human/Iran1-LA/kraken/db`.  Is this okay or should we use something different?

8. Use unicycler to assemble genome
```
        unicycler \
          -1 fastq_trimmed/Iran1_R1_paired.fastq.gz \
          -2 fastq_trimmed/Iran1_R2_paired.fastq.gz \
          --threads 20 --mode conservative --verbosity 2 \
          -o Iran1_assembled
```
9. QUAST relative to MN908947.3
```
        quast Iran1_assembled/assembly.fasta \
            -r /home/kmsmith/data/MN908947_3.fasta \
            -g /home/kmsmith/data/MN908947_3.gff3 \
            --output-dir Iran1_quast \
            --threads 20
```
- **Question:** MN908947.3 is currently hardcoded here. Is this okay or should we do something different?

10. 250bp tiled LMAT
```
        fatile Iran1_assembled/assembly.fasta 250 > Iran1_lmat/assembly.tiled.fasta

        # First few lines are docker boilerplate
        docker run \
          -v /home/kmsmith/data/LMAT-1.2.6/data:/data \
          -v /home/kmsmith/data/LMAT-1.2.6/runtime_inputs:/runtime_inputs \
          -v /home/kmsmith/git/covid-19-sequencing/pipeline/Iran1/Iran1_lmat:/pipeline \
          finlaymaguire/lmat:1.2.6 \
        bash /bin/run_rl.sh \   # "Real" LMAT command line starts here
          --db_file=/data/kML+Human.v4-14.20.g10.db \
          --query_file=/pipeline/assembly.tiled.fasta \
          --odir=/pipeline \
          --overwrite --verbose --threads=80

       cd Iran1_lmat
       parseLMAT > parseLMAT_output.txt
```
- **Question:** The LMAT DB `kML+Human.v4-14.20.g10.db` is currently hardcoded here. Is this okay or should we do something different?

11. HISAT2
```
      mkdir -p Iran1_hisat2
      hisat2-build Iran1_assembled/assembly.fasta Iran1_hisat2/genome
      
      hisat2 --threads 20 \
        -x Iran1_hisat2/genome \
        -1 fastq_trimmed/Iran1_R1_paired.fastq.gz \
        -2 fastq_trimmed/Iran1_R2_paired.fastq.gz \
        --summary-file Iran1_hisat2/summary.txt \
        -S Iran1_hisat2/output.sam
	
      samtools view --threads 20 \
         -b Iran1_hisat2/output.sam
	 > Iran1_hisat2/output.bam
```

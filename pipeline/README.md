
## Pipeline walkthrough

- For convenience in validating the pipeline, we list all shell commands
    executed by the pipeline below (cut-and-pasted from `snakemake -np`).

- This is for a running example with the following input files:
```
        fastq_inputs/MT-swab-Iran-Liverpool-pool1_S3_L001_R1_001.fastq.gz
        fastq_inputs/MT-swab-Iran-Liverpool-pool1_S3_L001_R2_001.fastq.gz
        fastq_inputs/MT-swab-Iran-Liverpool-pool2_S7_L001_R1_001.fastq.gz
        fastq_inputs/MT-swab-Iran-Liverpool-pool2_S7_L001_R2_001.fastq.gz
```

- Numbering of pipeline stages follows the following diagram.

![Workflow Version 5](../Workflow_Version_5.png)

### 1. Concatenate and sort R1 and R2 reads
```
        zcat sample1/fastq_inputs/MT-swab-Iran-Liverpool-pool1_S3_L001_R1_001.fastq.gz \
	     sample1/fastq_inputs/MT-swab-Iran-Liverpool-pool2_S7_L001_R1_001.fastq.gz \
	| paste - - - - \
	| sort -k1,1 -t " " \
	| tr "\t" "\n" \
	| gzip \
	> sample1/fastq_sorted/R1.fastq.gz

        zcat sample1/fastq_inputs/MT-swab-Iran-Liverpool-pool1_S3_L001_R2_001.fastq.gz \
	     sample1/fastq_inputs/MT-swab-Iran-Liverpool-pool2_S7_L001_R2_001.fastq.gz \
	| paste - - - - \
	| sort -k1,1 -t " " \
	| tr "\t" "\n" \
	| gzip \
	> sample1/fastq_sorted/R2.fastq.gz
```

### 2. Cutadapt primers + Trimmomatic
```
        cutadapt -j 16 \
          -a file:/home/kmsmith/data/wuhan_primers_28.01.20_trim_RC.fa
	  -A file:/home/kmsmith/data/wuhan_primers_28.01.20_trim_FW.fa \
	  -o sample1/fastq_primers_removed/R1.fastq.gz \   # R1 output file
	  -p sample1/fastq_primers_removed/R2.fastq.gz \   # R2 output file
	  sample1/fastq_sorted/R1.fastq.gz \               # R1 input file
	  sample1/fastq_sorted/R2.fastq.gz \               # R2 input file
	  >sample1/fastq_primers_removed/cutadapt.log

 	trimmomatic PE -threads 2 \
	  sample1/fastq_primers_removed/R1.fastq.gz \
	  sample1/fastq_primers_removed/R2.fastq.gz \
	  sample1/fastq_trimmed/R1_paired.fastq.gz \
	  sample1/fastq_trimmed/R1_unpaired.fastq.gz \
	  sample1/fastq_trimmed/R2_paired.fastq.gz \
	  sample1/fastq_trimmed/R2_unpaired.fastq.gz \
	  ILLUMINACLIP:/home/kmsmith/data/NexteraPE-PE.fa:2:30:10 \
	  SLIDINGWINDOW:4:20 \
	  2>sample1/fastq_trimmed/trim.log
```

- **Note:** Using Wuhan primer files here.

- **Note:** Using NexteraPE file here.

- **Question:** I'm probably being paranoid, but are we sure we're not swapping the `-a` and `-A` args to cutadapt?
     The cutadapt [documentation](https://cutadapt.readthedocs.io/en/stable/guide.html#trimming-paired-end-reads)
     suggests that we should be using "forward" primers for the `-a` argument, and "reverse" primers for `-A`.
     However, I had a hard time figuring out cutadapt's exact conventions (e.g. whether primers are read 5'-3'
     or vice versa) from its documentation, and I also wasn't sure what conventions were being used in our Wuhan files.

### 3. FastQC analysis of untrimmed reads
```
        fastqc sample1/fastq_sorted/R1.fastq.gz 2>sample1/fastq_sorted/R1_fastqc.log
        fastqc sample1/fastq_sorted/R2.fastq.gz 2>sample1/fastq_sorted/R2_fastqc.log
```

### 4. FastQC analysis of trimmed reads
```
        fastqc sample1/fastq_trimmed/R1_paired.fastq.gz 2>sample1/fastq_trimmed/R1_paired_fastqc.log
        fastqc sample1/fastq_trimmed/R2_paired.fastq.gz 2>sample1/fastq_trimmed/R2_paired_fastqc.log

        fastqc sample1/fastq_trimmed/R1_unpaired.fastq.gz 2>sample1/fastq_trimmed/R1_unpaired_fastqc.log
        fastqc sample1/fastq_trimmed/R2_unpaired.fastq.gz 2>sample1/fastq_trimmed/R2_unpaired_fastqc.log
```

**Note (minor):** currently running fastqc on both paired and unpaired reads, but planning on dropping unpaired reads,
 since these aren't actually used in the downstream pipeline.

### 5. Removing host DNA
```
        hisat2-build /home/kmsmith/data/MN908947_3.fasta \    # input reference genome
	  sample1/host_removed/sars-cov-2 \                   # output .ht2 prefix
	  >sample1/host_removed/hisat2-build.log 2>&1

        hisat2 --threads 2 \
	  -x sample1/host_removed/sars-cov-2 \                # input .ht2 prefix (from hisat-build)
	  -1 sample1/fastq_trimmed/R1_paired.fastq.gz \       # file with #1 mates (output from trimming, step 2)
	  -2 sample1/fastq_trimmed/R2_paired.fastq.gz \       # file with #2 mates (output from trimming, step 2)
	  --summary-file sample1/fastq_hist_removed/hisat2_summary.txt \
	  -S sample1/host_removed/mapped_and_unmapped.sam \   # SAM output file
	  2>sample1/host_removed/hisat2.log

        # Convert SAM to BAM
        samtools view -bS sample1/host_removed/mapped_and_unmapped.sam \
	  > sample1/host_removed/mapped_and_unmapped.bam

        # Select paired reads such that both ends map against reference genome
        samtools view -b -f 3 -F 4 sample1/host_removed/mapped_and_unmapped.bam \
	  > sample1/host_removed/both_ends_mapped.bam

        # Sort mapped reads by location in reference genome
	# (Needed for 'samtools mpileup' in step 10.)
        samtools sort -T tmp sample1/host_removed/both_ends_mapped.bam \
	  -o sample1/host_removed/both_ends_mapped_lsorted.bam

        # Sort mapped reads by name (needed for fastq conversion in next few lines)
	samtools sort -T tmp -n sample1/host_removed/both_ends_mapped.bam \
	  -o sample1/host_removed/both_ends_mapped_nsorted.bam

        # Convert to fastq
        bedtools bamtofastq -i sample1/host_removed/both_ends_mapped_nsorted.bam \
	  -fq sample1/host_removed/R1.fastq \
	  -fq2 sample1/host_removed/R2.fastq

        gzip sample1/host_removed/R1.fastq
        gzip sample1/host_removed/R2.fastq
```

- **Note:** using MN908947.3 reference genome here.

- **Note (minor):** should combine hisat2 and "sam to bam" into pipeline

### 6. Kraken2 analysis of trimmed reads (before removing host DNA), to derive percentage of reads derived from SARS-Cov-2 RNA
```
        cd sample1/kraken2
	
	kraken2 --db /home/kmsmith/data/Kraken2/db \
	  --threads 1 --quick --unclassified-out unclassified-sequences# --classified-out classified-sequences# \
	  --output kraken2.out \
	  --paired --gzip-compressed \
	  ../../sample1/fastq_trimmed/R1_paired.fastq.gz \
	  ../../sample1/fastq_trimmed/R2_paired.fastq.gz \
	  --report report \
	  2>../../sample1/kraken2/kraken2.log
```
- **Note:** Using the Kraken2 DB from `galaxylab:/home/raphenar/sars-cov-2/round1/human/Iran1-LA/kraken/db`.

### 7. Hisat2 confirmation of removal of human data, submit to NCBI BioSample

- **Placeholder:** nothing currently implemented here.

- The first step is to map against a human genome using hisat2,
  so we need a human genome `.ht2`, right?  (Or a human genome `.fasta`
  to run hisat2-build on.)  Can someone point me to this?
  

### 8-9. Mutation list relative to MN908948.3, and relative to clinical diagnostic primers
```
        breseq --reference /home/kmsmith/data/MN908947_primer_annotated_prot_clinical.gb \
	  --num-processors 1 --polymorphism-prediction --brief-html-output \
	  --output sample1/breseq \
	  sample1/host_removed/R1.fastq.gz \
	  sample1/host_removed/R2.fastq.gz \
	  >sample1/breseq/breseq.log 2>&1
```

- **Note:** for steps 8-9 we do a single breseq run, with genome `MN908947_primer_annotated_prot_clinical.gb`
  made by Jalees, which contains both the reference (MN908947.3) and clinical diagnostic primers.

### 10. Variant detection and consensus assembly
```
        samtools mpileup -A -d 0 \
            --reference /home/kmsmith/data/MN908947_3.fasta \
            -B -Q 0 \
            sample1/host_removed/both_ends_mapped_lsorted.bam \
        | ivar variants \
            -r /home/kmsmith/data/MN908947_3.fasta \
            -m 10 -p sample1/ivar_variants/ivar_variants -q 20 -t 0.25

        samtools mpileup -A -d 100000 -Q0 \
            sample1/host_removed/both_ends_mapped_lsorted.bam \
        | ivar consensus -t 0.75 -m 10 -n N \
            -p sample1/consensus/virus.consensus \
            2>sample1/consensus/ivar.log
```

- **Note:** using MN908947.3 reference genome here.

- **Question:** In the 'consensus' rule, the command-line help for `ivar consensus` recommends
 passing the `-aa` flag to `samtools mpileup`, but the ncov2019-artic workflow doesn't use this
 flag ([source](https://github.com/connor-lab/ncov2019-artic-nf/blob/master/modules/illumina.nf#L136)).
 Which is better?

### 11. Quast relative to MN908947.3
```
        quast sample1/consensus/virus.consensus.fa \
	    -r /home/kmsmith/data/MN908947_3.fasta \
	    -g /home/kmsmith/data/MN908947_3.gff3 \
	    --output-dir sample1/quast \
	    --threads 1 \
	>sample1/quast/quast.log
```
### 12. 250bp tiled LMAT
```
        perl ./fatile sample1/consensus/virus.consensus.fa 250 \
	    > sample1/lmat/assembly.tiled.fasta

        LMAT_DIR=/home/kmsmith/data/LMAT-1.2.6/runtime_inputs \   # LMAT needs this env variable
	bash $(which run_rl.sh)  \                  # LMAT needs to be manually run through bash
	    --db_file=/home/kmsmith/data/LMAT-1.2.6/data/kML+Human.v4-14.20.g10.db \
	    --query_file=sample1/lmat/assembly.tiled.fasta \
	    --odir=sample1/lmat \
	    --overwrite --verbose --threads=1

       cd sample1/lmat \
       perl ../../parseLMAT > parseLMAT_output.txt
```

- **Note:** Using LMAT DB `kML+Human.v4-14.20.g10.db`.

- **Note:** Now running LMAT via Fin's LMAT conda recipe, rather than a docker container.

### 13. Average depth of coverage against assembly (hisat2/bedtools)
```
       hisat2-build sample1/consensus/virus.consensus.fa \  # input reference genome (from 'ivar consensus', step 10)
           sample1/coverage/genome \                        # output .ht2 prefix
       >sample1/coverage/hisat2-build.log 2>&1

       hisat2 --threads 2 \
           -x sample1/coverage/genome \           # .ht2 prefix (from hisat-build in previous command)
	   -1 sample1/host_removed/R1.fastq.gz \  # file with #1 mates (from host removal, step 5)
	   -2 sample1/host_removed/R2.fastq.gz \  # file with #2 mates (from host removal, step 5)
	   --summary-file sample1/coverage/hisat2_summary.txt \
	   -S sample1/coverage/output.sam \       # SAM output file
	   2>sample1/coverage/hisat2.log

       # Convert SAM to BAM, with sort
       samtools view --threads 12 -b \
          sample1/coverage/output.sam \
       | samtools sort \
       > sample1/coverage/output.bam

       bedtools genomecov -d -ibam \
          sample1/coverage/output.bam \
       >sample1/coverage/depth.txt
```
- **Note:** Add postprocessing of depth file (eg compute min/max/median)? **Replaced ngsCAT with bedtools to generate per-base coverage.** 

- **Note:** On test data, I find that the first ~200 base pairs have zero coverage. Is that normal?
    Should we trim initial zeros?

### 14. QC control, submit to GSIAID / Nextstrain

- **Placeholder:** nothing currently implemented here

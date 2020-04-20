## Quick start:

A quick test run of the pipeline can be done as follows:

  - Get Fin's LMAT 1.2.6 docker container (`docker pull finlaymaguire/lmat:1.2.6`)

  - Create the directory `$HOME/data`, which should contain a copy of (or be a symlink to)
    `galaxylab:/home/kmsmith/data` (warning: 31 GB!)
    
  - From the `pipeline` subdirectory of this git repository, do:
    ```
    # Create pipeline
    ./c19_make_pipeline.py -o Iran1 $HOME/data/MT-swab-Iran-Liverpool*.fastq.gz

    # Run pipeline (cacheing conda envs in $HOME/.snakemake)
    cd Iran1/   # directory created by 'c19_make_pipeline.py'
    snakemake -p --cores=16 --use-conda --conda-prefix=$HOME/.snakemake all
    ```

## Software components:

  - fastqc 0.11.8 ([docs](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
  - cutadapt 1.18 ([docs](https://cutadapt.readthedocs.io/en/stable/)
  - trimmomatic 0.36 ([docs](http://www.usadellab.org/cms/?page=trimmomatic))
  - kraken2 2.0.7-beta ([docs](https://ccb.jhu.edu/software/kraken2/))
  - unicycler 0.4.8 ([github](https://github.com/rrwick/Unicycler))
  - quast 5.0.2 ([docs](http://quast.sourceforge.net/quast))
  - hisat 2.1.0 ([docs](http://daehwankimlab.github.io/hisat2/))
  - lmat, "sourceforge" version ([sourceforge](https://sourceforge.net/projects/lmat/))
  - samtools 1.7 ([docs](http://www.htslib.org/))
  - bedtools 2.26.0 ([docs](https://bedtools.readthedocs.io/en/latest/))

**Question:** These are the versions installed on galaxylab (with the exception of kraken2, see below).
 Are these the versions we want to use in all cases?

**Note:** I wanted to use kraken2 version 2.0.8-beta, but got segfaults with the bioconda version,
 so I used 2.0.7-beta instead.  I'll reinvestigate this later.
 
**Note:** We're using the "sourceforge" version of LMAT, even though predates the current github
 version by at least 4 years.  I plan to reinvestigate this too.

## Workflow diagram:

![Workflow Version 5](../Workflow_Version_5.png)
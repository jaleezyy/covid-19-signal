# SARS-CoV-2 Illumina GeNome Assembly Line (SIGNAL) 

This is a complete standardized workflow the assembly and subsequent analysis for short-read viral sequencing.
This core workflow is compatible with the [illumina artic nf pipeline](https://github.com/connor-lab/ncov2019-artic-nf) and produces the same consensus and variants using `ivar` (1.3) [Grubaugh, 2019](https://doi.org/10.1186/s13059-018-1618-7).
However, it performs far more extensive quality control and visualisation of results including an interactive HTML summary of run results.

Briefly, raw reads undergo qc using `fastqc` [Andrews](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) before removal of host-related reads by competitive mapping against a composite host and human reference with `BWA-MEM` (0.7.5) [Li, 2013](https://arxiv.org/abs/1303.3997), `samtools` (1.9) [Li, 2009](https://academic.oup.com/bioinformatics/article/25/16/2078/204688), and a [custom script](scripts/filter_non_human_reads.py).
This is to ensure raw as data as possible can be deposited in central databases.
After this, reads undergo adapter trimming and further qc with `trim-galore` (0.6.5) [Martin](https://doi.org/10.14806/ej.17.1.200).
Residual truseq sequencing adapters are then removed through another [custom script](scripts/filter_residual_adapters.py).
Reads are then mapped to the viral reference with `BWA-MEM`, and amplicon primer sequences trimmed using `ivar` (1.3) [Grubaugh, 2019](https://doi.org/10.1186/s13059-018-1618-7).
Fastqc is then used to perform a QC check on the reads that map to the viral reference.
After this, `ivar` is used to generate a consensus genome and variants are called using both `ivar variants` and `breseq` (0.35) [Deatherage, 2014](https://link.springer.com/protocol/10.1007/978-1-4939-0554-6_12).
Coverage statistics are calculated using bedtools before a final QC via `quast` and a `kraken2` taxonomic classification of mapped reads.
Finally, data from all samples are collated via a [post-processing script](scripts/signal_postprocess.py) into an interactive summary for exploration of results and quality control.
Optionally, users can run [ncov-tools](https://github.com/jts/ncov-tools/) to generate additional quality control and summary plots and statistics.

If you use this software please [cite](https://doi.org/10.3390/v12080895):

    Nasir, Jalees A., Robert A. Kozak, Patryk Aftanas, Amogelang R. Raphenya, Kendrick M. Smith, Finlay Maguire, Hassaan Maan et al. "A Comparison of Whole Genome Sequencing of SARS-CoV-2 Using Amplicon-Based Sequencing, Random Hexamers, and Bait Capture." Viruses 12, no. 8 (2020): 895.
    https://doi.org/10.3390/v12080895

## Setup/Execution

0. Clone the git repository (`--recursive` only needed to run`ncov-tools` postprocessing)
    
        git clone --recursive https://github.com/jaleezyy/covid-19-signal

1. Install `conda` and `snakemake` (version >5) e.g.

        wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
        bash Miniconda3-latest-Linux-x86_64.sh # follow instructions
        source $(conda info --base)/etc/profile.d/conda.sh
        conda create -n signal -c conda-forge -c bioconda -c defaults snakemake pandas
        conda activate signal 

There are some issues with `conda` failing to install newer versions of snakemake
so alternatively install `mamba` and use that (snakemake has beta support for it within the workflow)
    
        conda install -c conda-forge mamba
        mamba create -c conda-forge -c bioconda -n signal snakemake conda
        conda activate signal

Additional software dependencies are managed directly by `snakemake` using conda environment files:

  - trim-galore 0.6.5 ([docs](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
  - kraken2 2.0.7-beta ([docs](https://ccb.jhu.edu/software/kraken2/))
  - quast 5.0.2 ([docs](http://quast.sourceforge.net/quast))
  - bwa 0.7.17 ([docs](http://bio-bwa.sourceforge.net/))
  - samtools 1.7/1.9 ([docs](http://www.htslib.org/))
  - bedtools 2.26.0 ([docs](https://bedtools.readthedocs.io/en/latest/))
  - breseq 0.35.0 ([docs](https://barricklab.org/twiki/bin/view/Lab/ToolsBacterialGenomeResequencing))
  - ivar 1.3 ([docs](https://github.com/andersen-lab/ivar))
  - ncov-tools postprocessing scripts require additional dependencies (see [file](ncov-tools/workflow/envs/environment.yml)).

2. Download necessary database files

The pipeline requires:
 
 - Amplicon primer scheme sequences
 - SARS-CoV2 reference fasta
 - SARS-CoV2 reference gbk 
 - SARS-CoV2 reference gff3
 - kraken2 viral database
 - Human GRCh38 reference fasta (for composite human-viral BWA index)

        bash scripts/get_data_dependencies.sh -d data -a MN908947.3

3. Configure your `config.yaml` file

Either using the convenience python script or 
through modifying the `example_config.yaml` to suit your system

4. Specify your samples in CSV format (e.g. `sample_table.csv`)

See the example table `example_sample_table.csv` for an idea of how to organise this table. You can attempt to use `generate_sample_table.sh` to circumvent manual creation of the table.

5. Execute pipeline (optionally explicitly specify `--cores`)

      `snakemake -kp --configfile config.yaml --cores=NCORES --use-conda --conda-prefix=$PWD/.snakemake/conda all`
   
   If the `--conda-prefix` is not set as this then all envs will be reinstalled for each
   time you change the `results_dir` in the `config.yaml`.

6. Postprocessing analyses:

      `snakemake -p --configfile config.yaml --cores=NCORES --use-conda --conda-prefix=$PWD/.snakemake/conda postprocess`

After postprocessing finishes, you'll see the following summary files:

```
  - summary.html                top-level summary, with links to per-sample summaries
  - {sample_name}/sample.html   per-sample summaries, with links for more detailed info
  - {sample_name}/sample.txt    per-sample summaries, in text-file format instead of HTML
  - summary.zip                 zip archive containing all of the above summary files.
```
Note that the pipeline postprocessing ('snakemake postprocess') is separated from
the rest of the pipeline ('snakemake all').  This is because in a multi-sample run,
it's likely that at least one pipeline stage will fail.  The postprocessing script
should handle failed pipeline stages gracefully, by substituting placeholder values
when expected pipeline output files are absent.  However, this confuses snakemake's
dependency tracking, so there seems to be no good alternative to separating piepline
processing and postprocessing into 'all' and 'postprocess' targets.

Related: because pipeline stages can fail, we recommend running 'snakemake all'
with the -k flag ("Go on with independent jobs if a job fails").

Additionally, you can run @jts' [ncov-tools](https://github.com/jts/ncov-tools)
to generate phylogenies and alternative summaries.

    snakemake --use-conda --cores 10 ncov_tools

Signal manages installing the dependencies and will invoke `ncov-tools` if 
it has been cloned as a sub-module and a fasta containing sequences to include in 
the tree has been specified using `phylo_include_seqs:` in the main signal `config.yaml`.

Outputs will be written as specified within the `ncov-tools` folder and documentation.

### Docker

Alternatively, the pipeline can be deployed using Docker (see `resources/Dockerfile_pipeline` for specification).
To pull from dockerhub:

        docker pull finlaymaguire/signal

Download data dependencies into a data directory that already contains your reads (`data` is this example but whatever name you wish to use):

        mkdir -p data && docker run -v $PWD/data:/data finlaymaguire/signal:1.0.0 bash scripts/get_data_dependencies.sh -d /data

Generate your `config.yaml`and `sample_table.csv` (with paths to the readsets underneath `/data`) and place them into the data directory:

        cp config.yaml sample_table.csv $PWD/data

*WARNING* `result_dir` in `config.yaml` must be within `/data` e.g. `/data/results` to automatically be copied to your host system.  Otherwise they will be automatically deleted when the container finishes running (unless docker is run interactively).


Then execute the pipeline:

        docker run -v $PWD/data:/data finlaymaguire/signal conda run -n snakemake snakemake --configfile /data/config.yaml --use-conda --conda-prefix /covid-19-signal/.snakemake/conda --cores 8 all

## Summaries:

  - `postprocessing` and `ncov_tools` as described above generate many summaries including interactive html reports.`

  - Generate summaries of BreSeq among many samples, [see](resources/dev_scripts/summaries/README.md)


## Pipeline details:

For a step-by-step walkthrough of the pipeline, see [pipeline/README.md](PIPELINE.md).

A diagram of the workflow is shown below.

![Workflow Version 8](./resources/Workflow_Version_8.png)

## Possible Artefacts

- @jts: Host derived poly-A reads that sneak through the composite host removal stage can align to the viral poly-A tail giving it enough coverage to be called in the consensus.  Having this poly-A tail in the consensus can mess up the later analyses that require MSA.  If a sample is causing issues, check for host-derived poly-A reads. 


The conda environments defined here were created as follows (following suggestion by Fin in a slack message dated 2020-04-08):
```
conda create -n trim_qc -c conda-forge -c bioconda fastqc=0.11.8 cutadapt=1.18 trimmomatic=0.36 kraken2=2.0.7
conda create -n assembly -c conda-forge -c bioconda unicycler=0.4.8
conda create -n assembly_qc -c conda-forge -c bioconda quast=5.0.2
conda create -n snp_mapping -c conda-forge -c bioconda hisat2=2.1.0 breseq=0.35.0 samtools=1.7 bedtools=2.26.0
conda create -n lmat -c fmaguire lmat
```
Conda environments were exported to yaml with:
```
conda env export --name ENVNAME >ENVNAME.yaml
```
(Note: I also edited the yaml file to remove the `name:` at the top and the `prefix:` at the bottom, following the
[example in the snakemake documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management).)

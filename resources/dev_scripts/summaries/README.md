## Compare multiple samples breseqs

## Requirements

- breseq 0.29.0  ([docs](https://barricklab.org/twiki/pub/Lab/ToolsBacterialGenomeResequencing/documentation/))
- gdtools 0.29.0 ([docs](https://barricklab.org/twiki/pub/Lab/ToolsBacterialGenomeResequencing/documentation/gd_usage.html))

These can be installed using command: `conda install -c bioconda breseq`. Best method is to create a conda environment.

## Run breseq_compare.sh

Before running `breseq_compare.sh` provide paths to reference genome in genbank format (gb), and GenomeDiff or `.gd` files for multiple breseq runs. 

## QC plots for coronavirus sequencing data using ncov-tools

 - reference: https://github.com/jts/ncov-tools

 - run script `prepare_ncov-tools_summaries.sh` to generate the plots

 - provide absolute path to results from SIGNAL pipeline example:

 ```
 ./prepare_ncov-tools_summaries.sh /path/to/SIGNAL/results
 ```
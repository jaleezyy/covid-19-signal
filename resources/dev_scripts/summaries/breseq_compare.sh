## conda install -c bioconda breseq

## reference genomes
reference=MN908947_primer_annotated_prot_clinical.gb
## path to .gd files for all breseq runs
genomes=genomes/*.gd
## output name
output=compare

## html output
gdtools COMPARE \
--collapse \
--repeat-header 0 \
--format html \
--output ${output}.html \
--reference ${reference} ${genomes}

## tab output
gdtools COMPARE \
--collapse \
--repeat-header 0 \
--output ${output}.tsv \
--format tsv \
--reference ${reference} ${genomes}

## phylip output
gdtools COMPARE \
--collapse \
--phylogeny-aware \
--format phylip \
-o ${output}.phylip \
--reference ${reference} ${genomes}


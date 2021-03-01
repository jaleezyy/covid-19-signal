#!/bin/bash

# Need scheme to grab the correct files and config
scheme=$1

# Need the working directory to properly do stuff
work_dir=$2

# Activate env
eval "$(conda shell.bash hook)"

# Default name is ncov-qc from https://github.com/jts/ncov-tools/blob/master/workflow/envs/environment.yml
# May have to switch to a container or something similar later to make it easier
conda activate ncov-qc

# Path to scheme bed files
if [ "$scheme" == "articV3" ]
then
    echo "amplicon_bed: $work_dir/resources/primer_schemes/artic_v3/ncov-qc_V3.scheme.bed" >> ./ncov-tools/config.yaml
    echo "primer_bed: $work_dir/resources/primer_schemes/artic_v3/nCoV-2019.bed" >> ./ncov-tools/config.yaml

elif [ "$scheme" == "freed" ]
then
    echo "amplicon_bed: $work_dir/resources/primer_schemes/freed/ncov-qc_freed.scheme.bed" >> ./ncov-tools/config.yaml
    echo "primer_bed: $work_dir/resources/primer_schemes/freed/nCoV-2019.bed" >> ./ncov-tools/config.yaml

elif [ "$scheme" == "resende" ]
then
    echo "amplicon_bed: $work_dir/resources/primer_schemes/2kb_resende/ncov-qc_resende.scheme.bed" >> ./ncov-tools/config.yaml
    echo "primer_bed: $work_dir/resources/primer_schemes/2kb_resende/nCoV-2019.bed" >> ./ncov-tools/config.yaml

elif [ "$scheme" == "V2resende" ]
then
    echo "amplicon_bed: $work_dir/resources/primer_schemes/2kb_resende_v2/nCoV-2019.bed" >> ./ncov-tools/config.yaml
    echo "primer_bed: $work_dir/resources/primer_schemes/2kb_resende_v2/ncov-qc_resende.scheme.bed" >> ./ncov-tools/config.yaml

else
    echo "Please specify 'articV3', 'freed', 'resende', or 'V2resende' as a scheme"
    exit;
fi

# If we have a matching negative control, we modify the config to make sure its gotten
# If we don't find any, then no negative controls are added
if $(ls signal_results/ | grep -q -i "negative\|ntc\|water\|blank\|neg")
then
    ls -1 signal_results/ > name_list.txt
    negative_list=$(grep -i -e ntc -e negative -e water -e blank -e neg name_list.txt | cut -f 1 | sed 's/^/"/g' | sed 's/$/"/g' | tr "\n" ',' | sed 's/^/[/' | sed 's/$/]/')
    echo "negative_control_samples: ${negative_list}" >> ./ncov-tools/config.yaml
fi

# Configs gotten and set, now we need the files!
# Files are: *.consensus.fa, *_ivar_variants.tsv, and *_viral_reference.mapping.primertrimmed.sorted.bam
# File names are consistent so we can easily do this!
mkdir -p ncov-tools/files/
find ./signal_results/ -type f -name *.consensus.fa -exec cp {} ncov-tools/files/ \;
find ./signal_results/ -type f -name *.sorted.bam* -exec cp {} ncov-tools/files/ \;
find ./signal_results/ -type f -name *variants.tsv -exec cp {} ncov-tools/files/ \;

# Rename files to give consistent structure that is easier to work with
cd ./ncov-tools/files
$work_dir/scripts/prename "s/_viral_reference././" *.sorted.bam
$work_dir/scripts/prename "s/_ivar_/./" *_variants.tsv
# Remove extra header info
sed -i 's|.consensus_threshold_0.75_quality_20||' *.consensus.fa

# Back to ncov-tools root directory to run it
cd ../
snakemake -kp -s workflow/Snakefile all --cores 3
snakemake -s workflow/Snakefile --cores 1 build_snpeff_db
snakemake -s workflow/Snakefile --cores 2 all_qc_annotation

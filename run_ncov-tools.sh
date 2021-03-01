#!/bin/bash

# Need scheme to grab the correct files and config
scheme=$1

# Need the working directory to properly do stuff
work_dir=$2

# Name to watch on sq
run_name=$3

# Activate env
eval "$(conda shell.bash hook)"

conda activate ncov-tools

# Path to scheme bed files
if [ "$scheme" == "articV3" ]
then
    echo "amplicon_bed: $work_dir/resources/primer_schemes/artic_v3/ncov-qc_V3.scheme.bed" >> ./ncov-tools/config.yaml
    echo "primer_bed: $work_dir/resources/primer_schemes/artic_v3/nCoV-2019.bed" >> ./ncov-tools/config.yaml

elif [ "$scheme" == "freed" ]
then
    echo "amplicon_bed: $work_dir/resources/primer_schemes/Freed.distinct_regions.bed" >> ./ncov-tools/config.yaml
    echo "primer_bed: $work_dir/resources/primer_schemes/Freed_signal.bed" >> ./ncov-tools/config.yaml

elif [ "$scheme" == "resende" ]
then
    echo "amplicon_bed: $work_dir/resources/primer_schemes/Resende.distinct_regions.bed" >> ./ncov-tools/config.yaml
    echo "primer_bed: $work_dir/resources/primer_schemes/signal_resende_fixed.bed" >> ./ncov-tools/config.yaml

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
find ./signal_results/ -type f -name *.consensus.fa -exec cp {} ncov-tools/files/ \;
find ./signal_results/ -type f -name *.sorted.bam* -exec cp {} ncov-tools/files/ \;
find ./signal_results/ -type f -name *variants.tsv -exec cp {} ncov-tools/files/ \;

# Now hopefully we don't need to re-name and if we do it will be below...
# Needed to rename :(
cd ./ncov-tools/files
#<placeholder>prename "s/_viral_reference././" *.sorted.bam
#<placeholder>prename "s/_ivar_/./" *_variants.tsv
sed -i 's|.consensus_threshold_0.75_quality_20||' *.consensus.fa

# And then finally run it!!!!!
# back to ncov-tools root directory
cd ../
snakemake -kp -s workflow/Snakefile all --cores 12
snakemake -s workflow/Snakefile --cores 1 build_snpeff_db
snakemake -s workflow/Snakefile --cores 2 all_qc_annotation

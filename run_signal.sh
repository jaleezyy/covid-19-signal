#!/bin/bash

################
### SETTINGS ###
################

# DEFAULTS #
############
set -e # exit if pipeline returns non-zero status
set -o pipefail # return value of last command to exit with non-zero status
SCRIPTPATH="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )" # Get script location which will be the signal base dir

# Function to check if element is in an array
containsElement () {
    local e match="$1"
    shift
    for e; do [[ "$e" == "$match" ]] && return 0; done
    return 1
}

# Values to Check #
schemeArray=('articV3' 'freed' 'resende' 'V2resende')
envArray=('signal' 'ncov-qc' 'snpdist_signal' 'qc_summary_signal')

# Base Parameters #
FASTQ_PAIRS=""
PRIMER_SCHEME=""
CORES="3"
RUNNAME="nml"
### END DEFAULTS ###

# INPUTS #
##########
while getopts ":hf:p:c:n:" opt; do
  case ${opt} in
    h )
      echo "Usage:"
      echo "    bash run_signal.sh -h                      Display this help message."
      echo "    bash $SCRIPTPATH/run_signal.sh -f PATH_TO_PAIRED_FASTQ_DIR -p PRIMER_SCHEME
    
Flags:
    -f      :  Path to paired fastq file directory
    -p      :  Specify input data primer scheme
                Available Primer Schemes: articV3, freed, resende, V2resende
    -c      :  (OPTIONAL) Number of Cores to use in Signal. Default is 3
    -n      :  (OPTIONAL) Run name for final outputs. Default is 'nml'
    "
        exit 0
        ;;
    f)
        FASTQ_PAIRS=${OPTARG%/}
            if [ -d "$FASTQ_PAIRS" ]; then
                echo "Directory '$FASTQ_PAIRS' exists"
            else
                echo "ERROR: Directory '$FASTQ_PAIRS' does not exist"
                exit 1
            fi       
        ;;
    p)
        PRIMER_SCHEME=$OPTARG
            if containsElement "$PRIMER_SCHEME" "${schemeArray[@]}"; then
                echo "Using primer scheme $PRIMER_SCHEME"
            else
                echo "ERROR: $PRIMER_SCHEME unavailable"
                echo "Primer schemes available are ${schemeArray[@]}"
                exit 1
            fi
        ;;
    t)
        CORES=$OPTARG
            if [[ $CORES == +([0-9]) ]]; then
                echo "Using $CORES"
            else
                echo "ERROR: Cores input (-c) not an integer"
                exit 1
            fi
        ;;
    n)
        RUNNAME=$OPTARG
        ;;
   \? )
     echo "Invalid Option: -$OPTARG" 1>&2
     exit 1
     ;;
  esac
done
### END INPUTS ###


################################
### OVERALL AUTOMATION SETUP ###
################################

# PATHES #
##########

# Set name for output and full path to reads #
pushd $FASTQ_PAIRS
run_name=${PWD##*/}
fastq_dir_path=$PWD
popd

# If we cannot find the data directory, we will fail out and let user know it isn't there
if [ ! -d "$SCRIPTPATH/data/" ]; then
    echo "data/ directory not found in $SCRIPTPATH"
    echo "Please run bash $SCRIPTPATH/scripts/get_data_dependencies.sh -d data -a MN908947.3 to make this folder or move an already made data folder here"
    exit 1
fi
### END PATHES ###

# CONDA #
#########
eval "$(conda shell.bash hook)"

# Conda Env Check #
for ENV in ${envArray[@]}; do
    if [ $(conda env list | awk '{print $1}' | grep "^$ENV"'$') ]; then
        echo "$ENV found"
    else
        echo "Conda env '$ENV' doesn't exist."
        # If its the ncov-tools env, need to use mamba
        if [ $ENV = "ncov-qc" ]; then
            conda install -y mamba -n base -c conda-forge
            mamba env create -f $SCRIPTPATH/conda_envs/$ENV.yml
        else
            conda env create --file $SCRIPTPATH/conda_envs/$ENV.yml
        fi
    fi
done
### END CONDA ###

# FINAL OUTPUT LOCATION #
#########################
timestamp=`date +%b%d_%H%M`
root="signal_${run_name}_${timestamp}"
mkdir -p $root
cd $root

# Get output root full path for later parts to the root directory #
root_path=$PWD
### END OUTPUT LOCATION ###

# ALL SNAKEMAKE CONFIGURATIONS #
################################
SIGNAL_PROFILE="$SCRIPTPATH/resources/profile"
cp "$SCRIPTPATH/resources/profile/parameters.yaml" .
SIGNAL_CONFIG="parameters.yaml"

# Get always current ncov-tools with clone and copy resources as we write to the config
echo "Setting up files in $root_path"
git clone --depth 1 https://github.com/jts/ncov-tools
# Copy the ncov-tools config as we will write data to it later
cp -r $SCRIPTPATH/resources/ncov-tools_files/* ./ncov-tools/

# Echo in correct data locations to parameters config
# Done as for the primer pairs to work we need the full path
if [ "$PRIMER_SCHEME" == "articV3" ]; then
    # SIGNAL Parameters
    echo "scheme_bed: 'resources/primer_schemes/artic_v3/nCoV-2019.primer.bed'" >> $SIGNAL_CONFIG
    echo "amplicon_loc_bed: 'resources/primer_schemes/artic_v3/ncov-qc_V3.scheme.bed'" >> $SIGNAL_CONFIG
    echo "primer_pairs_tsv: '-f $SCRIPTPATH/resources/primer_pairs/articV3_primer_pairs.tsv'" >> $SIGNAL_CONFIG

    # NCOV-TOOLS Parameters
    echo "amplicon_bed: $SCRIPTPATH/resources/primer_schemes/artic_v3/ncov-qc_V3.scheme.bed" >> ./ncov-tools/config.yaml
    echo "primer_bed: $SCRIPTPATH/resources/primer_schemes/artic_v3/nCoV-2019.bed" >> ./ncov-tools/config.yaml

elif [ "$PRIMER_SCHEME" == "freed" ]; then
    # SIGNAL Parameters
    echo "scheme_bed: 'resources/primer_schemes/freed/nCoV-2019.primer.bed'" >> $SIGNAL_CONFIG
    echo "amplicon_loc_bed: 'resources/primer_schemes/freed/ncov-qc_freed.scheme.bed'" >> $SIGNAL_CONFIG
    echo "primer_pairs_tsv: '-f $SCRIPTPATH/resources/primer_pairs/freed_primer_pairs.tsv'" >> $SIGNAL_CONFIG

    # NCOV-TOOLS Parameters
    echo "amplicon_bed: $SCRIPTPATH/resources/primer_schemes/freed/ncov-qc_freed.scheme.bed" >> ./ncov-tools/config.yaml
    echo "primer_bed: $SCRIPTPATH/resources/primer_schemes/freed/nCoV-2019.bed" >> ./ncov-tools/config.yaml

elif [ "$PRIMER_SCHEME" == "resende" ]; then
    # SIGNAL Parameters
    echo "scheme_bed: 'resources/primer_schemes/2kb_resende/nCoV-2019.primer.bed'" >> $SIGNAL_CONFIG
    echo "amplicon_loc_bed: 'resources/primer_schemes/2kb_resende/ncov-qc_V3.scheme.bed'" >> $SIGNAL_CONFIG
    echo "primer_pairs_tsv: '-f $SCRIPTPATH/resources/primer_pairs/resende_primer_pairs.tsv'" >> $SIGNAL_CONFIG

    # NCOV-TOOLS Parameters
    echo "amplicon_bed: $SCRIPTPATH/resources/primer_schemes/2kb_resende/ncov-qc_resende.scheme.bed" >> ./ncov-tools/config.yaml
    echo "primer_bed: $SCRIPTPATH/resources/primer_schemes/2kb_resende/nCoV-2019.bed" >> ./ncov-tools/config.yaml

elif [ "$PRIMER_SCHEME" == "V2resende" ]; then
    # SIGNAL Parameters
    echo "scheme_bed: 'resources/primer_schemes/2kb_resende_v2/nCoV-2019.primer.bed'" >> $SIGNAL_CONFIG
    echo "amplicon_loc_bed: 'resources/primer_schemes/2kb_resende_v2/ncov-qc_V3.scheme.bed'" >> $SIGNAL_CONFIG
    echo "primer_pairs_tsv: '-f $SCRIPTPATH/resources/primer_pairs/resende_primer_pairs.tsv'" >> $SIGNAL_CONFIG

    # NCOV-TOOLS Parameters
    echo "amplicon_bed: $SCRIPTPATH/resources/primer_schemes/2kb_resende_v2/nCoV-2019.bed" >> ./ncov-tools/config.yaml
    echo "primer_bed: $SCRIPTPATH/resources/primer_schemes/2kb_resende_v2/ncov-qc_resende.scheme.bed" >> ./ncov-tools/config.yaml
fi
echo "run_name: $RUNNAME" >> ./ncov-tools/config.yaml
### END SNAKEMAKE CONFIGS ###


# GET WANTED SIGNAL #
ln -s $SCRIPTPATH/* .


##################
### RUN SIGNAL ###
##################
# Activate env
# Echo out relevant info to double check it looks ok
conda activate signal
echo "Running ${run_name} with scheme ${PRIMER_SCHEME}"
echo "config: $SIGNAL_CONFIG"
echo "profile: $SIGNAL_PROFILE"

# Run signal and postprocessing
bash generate_sample_table.sh -d $fastq_dir_path
snakemake -s Snakefile --configfile $SIGNAL_CONFIG --profile $SIGNAL_PROFILE --cores=$CORES --conda-prefix=$SCRIPTPATH/.snakemake/conda all
snakemake -s Snakefile --configfile $SIGNAL_CONFIG --profile $SIGNAL_PROFILE --cores=$CORES --conda-prefix=$SCRIPTPATH/.snakemake/conda postprocess


##################
### NCOV-TOOLS ###
##################
# Run NCOV-Tools on separate script at the moment as it won't work on same one. IDK why the same code won't work here
bash nml_automation/run_ncov-tools.sh $CORES $SCRIPTPATH


##################
### FINAL DATA ###
##################

# Set Up #
##########
cd $root_path
mkdir -p ./summary_csvs

# Setting up relative paths from the root dir to needed data
ncov_qc="./ncov-tools/qc_reports/${RUNNAME}_summary_qc.tsv"
ncov_neg="./ncov-tools/qc_reports/${RUNNAME}_negative_control_report.tsv"
ref="./data/MN908947.3.fasta"
pangolin="./ncov-tools/lineages/${RUNNAME}_lineage_report.csv"

# Check to see if we have a negative control. If so do nothing. If not make a "fake" one to not error out
if [ -f "$ncov_neg" ];
then
    echo "Negative contol data found"
else
    touch ./ncov-tools/qc_reports/${RUNNAME}_negative_control_report.tsv
fi
### END Set Up ###

# QC Summary #
##############
conda activate qc_summary_signal

for RESULT in $(ls -d signal_results/*/)
do
    # Getting and keeping track of the name
    found_name="$(basename $RESULT)"
    relative_sample_path="$RESULT/core/"
    snpeff="./ncov-tools/qc_annotation/${found_name}_aa_table.tsv"
    echo "Processing $found_name"

    # Run the python summary
    python ./nml_automation/qc.py --illumina \
        --outfile ./summary_csvs/${found_name}.qc.csv \
        --sample ${found_name} \
        --ref ${ref} \
        --bam ${relative_sample_path}${found_name}_viral_reference.mapping.primertrimmed.sorted.bam \
        --fasta ${relative_sample_path}${found_name}.consensus.fa \
        --tsv_variants ${relative_sample_path}${found_name}_ivar_variants.tsv \
        --pangolin ${pangolin} \
        --ncov_summary ${ncov_qc} \
        --ncov_negative ${ncov_neg} \
        --revision v1.2.1 \
        --script_name covid-19-signal \
        --sequencing_technology illumina \
        --scheme $PRIMER_SCHEME \
        --snpeff_tsv $snpeff \
        --pcr_bed ./resources/pcr_primers.bed
done

# Concatenate all summaries
csvtk concat ./summary_csvs/*.qc.csv > ./summary_csvs/combined.csv
# Check negative controls
python ./nml_automation/negative_control_fixes.py --qc_csv ./summary_csvs/combined.csv --output_prefix $RUNNAME

# SNP DISTS #
#############
conda deactivate
conda activate snpdist_signal

if [ -f "./ncov-tools/qc_analysis/${RUNNAME}_aligned.fasta" ];
then
    snp-dists ./ncov-tools/qc_analysis/${RUNNAME}_aligned.fasta > matrix.tsv
else
    echo "No ./ncov-tools/qc_analysis/${RUNNAME}_aligned.fasta, skipping the snp-dist check"
fi


### END RUNNING ###

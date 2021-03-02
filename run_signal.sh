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

# Allowed scheme values #
schemeArray=('articV3' 'freed' 'resende' 'V2resende')

# Base Parameters #
FASTQ_PAIRS=""
PRIMER_SCHEME=""
SIGNAL_DIR="$SCRIPTPATH"
THREADS="3"
### END DEFAULTS ###

# INPUTS #
##########
while getopts ":hf:p:s:t:" opt; do
  case ${opt} in
    h )
      echo "Usage:"
      echo "    bash run_signal.sh -h                      Display this help message."
      echo "    bash $SCRIPTPATH/run_signal.sh -f PATH_TO_PAIRED_FASTQ_DIR -p PRIMER_SCHEME -s PATH_TO_SIGNAL_DIRECTORY
    
Flags:
    -f      :  Path to paired fastq file directory
    -p      :  Specify input data primer scheme
                Available Primer Schemes: articV3, freed, resende, V2resende
    -s      :  (OPTIONAL) Path to wanted base Signal Directory. Default is $SCRIPTPATH
    -t      :  (OPTIONAL) Number of Threads to use in Signal. Default is 3
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
    s) 
        SIGNAL_DIR=${OPTARG%/}
            if [ -d "$SIGNAL_DIR" ]; then
                echo "Directory '$SIGNAL_DIR' exists"
            else
                echo "ERROR: Directory '$SIGNAL_DIR' does not exist"
                exit 1
            fi
        ;;
    t)
        THREADS=$OPTARG
            if [[ $THREADS == +([0-9]) ]]; then
                echo "Using $THREADS"
            else
                echo "ERROR: Threads input (-t) not an integer"
                exit 1
            fi
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
pushd $SIGNAL_DIR
SIGNAL_DIR=$PWD
popd

# Set name for output and full path to reads #
pushd $FASTQ_PAIRS
run_name=${PWD##*/}
fastq_dir_path=$PWD
popd

# If we cannot find the data directory, we will fail out and let user know it isn't there
if [ ! -d "$SIGNAL_DIR/data/" ]; then
    echo "data/ directory not found in $SIGNAL_DIR"
    echo "Please run bash $SIGNAL_DIR/scripts/get_data_dependencies.sh -d data -a MN908947.3 to make this folder"
    exit 1
fi
### END PATHES ###

# CONDA #
#########
eval "$(conda shell.bash hook)"

# Env check
# Write check for if we have correct conda envs and if not, create them
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
### END SNAKEMAKE CONFIGS ###


# GET WANTED SIGNAL #
ln -s $SIGNAL_DIR/* .


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
snakemake -s Snakefile --configfile $SIGNAL_CONFIG --profile $SIGNAL_PROFILE --cores=$THREADS --conda-prefix=$SIGNAL_DIR/.snakemake/conda all
snakemake -s Snakefile --configfile $SIGNAL_CONFIG --profile $SIGNAL_PROFILE --cores=$THREADS --conda-prefix=$SIGNAL_DIR/.snakemake/conda postprocess


##################
### NCOV-TOOLS ###
##################
# Run NCOV-Tools on separate script at the moment as it won't work on same one? IDK why
bash nml_automation/run_ncov-tools.sh $PRIMER_SCHEME $THREADS $SCRIPTPATH

# Post NCOV-TOOLS run processing and uploads of data
# We are back in the root dir. Need the ncov-tools data along with the signal data for final output table
# individual data to the ./summary_csvs directory made by script
# bash $SIGNAL_DIR/nml_automation/final_cleanup.sh

### END RUNNING ###

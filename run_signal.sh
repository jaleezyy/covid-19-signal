#!/bin/bash

### Settings ###
set -e # exit if pipeline returns non-zero status
set -o pipefail # return value of last command to exit with non-zero status

FASTQ_PAIRS=""
PRIMER_SCHEME=""
SIGNAL_DIR=""


while getopts ":hf:p:s:n:" opt; do
  case ${opt} in
    h )
      echo "Usage:"
      echo "    run_prov_signal.sh -h                      Display this help message."
      echo "    run_prov_signal.sh -f <directory> -p <name> -s <directory> -n <directory>      USAGE: run_prov_signal.sh -f PATH_TO_PAIRED_FASTQ_DIR -p PRIMER_SCHEME -s PATH_TO_SIGNAL_DIRECTORY
    
Flags:
    -f      :  Path to paired fastq file directory
    -p      :  Specify input data primer scheme
                Available Primer Schemes: articV3, freed, resende, V2resende
    -s      :  Path to signal directory
    "
        exit 0
        ;;
    f)
        FASTQ_PAIRS=${OPTARG%/}
            if [ -d "$FASTQ_PAIRS" ]; then
                echo "Directory '$FASTQ_PAIRS' exists"
            else
                echo "Directory '$FASTQ_PAIRS' does not exists"
            fi       
        ;;
    p)
        PRIMER_SCHEME=$OPTARG
        ;;
    s) 
        SIGNAL_DIR=${OPTARG%/}
        if [ -d "$SIGNAL_DIR" ]; then
                echo "Directory '$SIGNAL_DIR' exists"
            else
                echo "Directory '$SIGNAL_DIR' does not exists"
            fi
        ;;
   \? )
     echo "Invalid Option: -$OPTARG" 1>&2
     exit 1
     ;;
  esac
done


pushd $SIGNAL_DIR
SIGNAL_DIR=$PWD
popd

SIGNAL_PROFILE=""
SIGNAL_CONFIG=""
if [ "$PRIMER_SCHEME" == "articV3" ];
then
    SIGNAL_PROFILE=$SIGNAL_DIR/resources/signal_profiles/articV3
    SIGNAL_CONFIG=$SIGNAL_DIR/resources/signal_profiles/articV3/parameters.yaml
elif [ "$PRIMER_SCHEME" == "freed" ];
then
    SIGNAL_PROFILE=$SIGNAL_DIR/resources/signal_profiles/freed
    SIGNAL_CONFIG=$SIGNAL_DIR/resources/signal_profiles/freed/parameters.yaml
elif [ "$PRIMER_SCHEME" == "resende" ];
then
    SIGNAL_PROFILE=$SIGNAL_DIR/resources/signal_profiles/resende
    SIGNAL_CONFIG=$SIGNAL_DIR/signal_profiles/resende/parameters.yaml
elif [ "$PRIMER_SCHEME" == "V2resende" ];
then
    SIGNAL_PROFILE=$SIGNAL_DIR/resources/signal_profiles/V2resende
    SIGNAL_CONFIG=$SIGNAL_DIR/resources/signal_profiles/V2resende/parameters.yaml
else
    echo "ERROR: $PRIMER_SCHEME is an invalid Primer Scheme"
    echo "$HELP"
    exit 1
fi

### Set Up Paths and File Structure ###
#######################################
# Define our path to signal directory #
if [ "$SIGNAL_DIR" = "" ];
then
	echo "$SIGNAL_DIR : ERROR: signal directory does not exist"
	echo "$HELP"
	exit 1
fi

# Set name for output and full path to reads #
pushd $FASTQ_PAIRS
run_name=${PWD##*/}
fastq_dir_path=$PWD
popd


# Setup Output Folder structure #
timestamp=`date +%b%d_%H%M`

# Load Conda
eval "$(conda shell.bash hook)"

root="signal_${run_name}_${timestamp}"
mkdir -p $root
cd $root

# Get output root full path for later parts to the root directory #
root_path=$PWD

### GET INPUT FILES ###
#######################
# Get always current ncov-tools with clone and copy resources as we write to the config
echo "Setting up files in $PWD"
git clone --depth 1 https://github.com/jts/ncov-tools
# Copy the ncov-tools config as we will write data to it later
cp -r $SIGNAL_DIR/resources/ncov-tools_files/* ./ncov-tools/

# Softlinking signal data and our automation scripts to run them but also to save space
ln -s $SIGNAL_DIR/* .

# Generate the sample_list.csv
bash generate_sample_table.sh -d $fastq_dir_path
### END FILE SETUP ###

### RUN SIGNAL AND QC ###
#########################
# Activate env
# Echo out relevant info to double check it looks ok
conda activate signal
echo "Running ${run_name} with scheme ${PRIMER_SCHEME}"
echo "config: $SIGNAL_CONFIG"
echo "profile: $SIGNAL_PROFILE"

# Run signal and postprocessing
snakemake -s Snakefile --configfile $SIGNAL_CONFIG --profile $SIGNAL_PROFILE --cores=3 --conda-prefix=$SIGNAL_DIR/.snakemake/conda all
snakemake -s Snakefile --configfile $SIGNAL_CONFIG --profile $SIGNAL_PROFILE --cores=3 --conda-prefix=$SIGNAL_DIR/.snakemake/conda postprocess

# Setup and run NCOV-TOOLS
bash nml_automation/run_ncov-tools.sh $PRIMER_SCHEME $root_path

# Post NCOV-TOOLS run processing and uploads of data
# We are back in the root dir. Need the ncov-tools data along with the signal data for final output table
# individual data to the ./summary_csvs directory made by script
# bash $SIGNAL_DIR/nml_automation/final_cleanup.sh

### END RUNNING ###

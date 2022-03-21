#!/bin/bash

################
### SETTINGS ###
################

# DEFAULTS #
############
SCRIPTPATH="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )" # Get script location which will be the signal base dir
mkdir -p $SCRIPTPATH/.snakemake

# Function to check if element is in an array
containsElement () {
    local e match="$1"
    shift
    for e; do [[ "$e" == "$match" ]] && return 0; done
    return 1
}

# Base Parameters #
FASTQ_PAIRS=0
PRIMER_SCHEME=0
CORES="3"
RUNNAME="nml"
METADATA_TSV=0
PDF=false

# SUBSTITUTE PARAMS #
BED=0
AMPLICON=0
CUSTOM=false

# Conda Parameters #
BASE_ENV_PATH="$SCRIPTPATH/.snakemake/conda"
SIGNAL_ENV="signal"
USER_SIGNAL=false
NCOV_ENV="ncov-qc-pangolin3"
USER_NCOV=false

# Values to Check #
schemeArray=('articV3' 'freed' 'resende' 'V2resende' 'articV4')

HELP="
USAGE:
    bash $SCRIPTPATH/run_signal.sh -d PATH_TO_PAIRED_FASTQ_DIR -p PRIMER_SCHEME <OPTIONAL FLAGS>
    bash $SCRIPTPATH/run_signal.sh --update
    
Flags:
    NEEDED:
    -d  --directory      :  Path to paired fastq file directory
    -p  --primer-scheme  :  Specify input data primer scheme
                Available Primer Schemes: articV3, articV4, freed, resende, V2resende, custom by passing '--bed' and '--amplicon'

    SUBSTITUTE (Can be used instead of a primer scheme but must be used together):
    --bed       :  Path to custom primer bed file to be used instead of default schemes
    --amplicon  :  Path to custom amplicon bed file to be used instead of default schemes

    OPTIONAL:
    -c  --cores           :  Number of Cores to use in Signal. Default is 3
    -n  --run-name        :  Run name for final ncov-tools outputs. Default is 'nml'
    -m  --metadata        :  Add metadata to the run. Must be in TSV format with atleast a column called 'sample'
    --pdf                 :  If you have pdflatex installed runs ncov-tools pdf output
    --signal-env          :  Name of signal conda env. Default is '$BASE_ENV_PATH/signal'
    --ncov-tools-env      :  Name of ncov-tools env. Default is '$BASE_ENV_PATH/ncov-qc-pangolin3'
                **NOTE** It is highly recommended to let the script generate the environments as it will
                          only occur once and you won't have to pass the path each time

    OTHER:
    --update  :  Passing --update will update ncov-tools pip dependencies, pangolin and pangoLEARN along with this repo and then exit
"
### END DEFAULTS ###

# INPUTS #
##########

# Check for Args #
if [ $# -eq 0 ]; then
    echo "$HELP"
    exit 0
fi

# Set Arguments #
while [ "$1" = "--directory" -o "$1" = "-d" -o "$1" = "--primer-scheme" -o "$1" = "-p" -o "$1" = "--cores" -o "$1" = "-c" -o "$1" = "--bed" -o "$1" = "--amplicon" -o "$1" = "--run-name" -o "$1" = "-n" -o "$1" = "-m" -o "$1" = "--metadata" -o "$1" = "--pdf" -o "$1" = "--signal-env" -o "$1" = "--ncov-tools-env" -o "$1" = "--update" ];
do
    if [ "$1" = "--directory" -o "$1" = "-d" ]; then
        shift
        FASTQ_PAIRS=$1
        shift
    elif [ "$1" = "--primer-scheme" -o "$1" = "-p" ]; then
        shift
        PRIMER_SCHEME=$1
        shift
    elif [ "$1" = "--cores" -o "$1" = "-c" ]; then
        shift
        CORES=$1
        shift
    elif [ "$1" = "--bed" ]; then
        shift
        BED=$1
        shift
    elif [ "$1" = "--amplicon" ]; then
        shift
        AMPLICON=$1
        shift
    elif [ "$1" = "--run-name" -o "$1" = "-n" ]; then
        shift
        RUNNAME=$1
        shift
    elif [ "$1" = "--metadata" -o "$1" = "-m" ]; then
        shift
        METADATA_TSV=$1
        shift
    elif [ "$1" = "--pdf" ]; then
        PDF=true
        shift
    elif [ "$1" = "--signal-env" ]; then
        shift
        SIGNAL_ENV=$1
        USER_SIGNAL=true
        shift
    elif [ "$1" = "--ncov-tools-env" ]; then
        shift
        NCOV_ENV=$1
        USER_NCOV=true
        shift
    elif [ "$1" = "--update" ]; then
        shift
        # Scripts
        cd $SCRIPTPATH
        git pull

        # ncov-tools Environment (not managed by snakemake unfortunately)
        eval "$(conda shell.bash hook)"
        printf "\n Updating ncov-tools environment at $BASE_ENV_PATH/$NCOV_ENV \n\n"
        conda activate $BASE_ENV_PATH/$NCOV_ENV
        pip install git+https://github.com/cov-lineages/pango-designation.git --upgrade
        pangolin --update
        pip install ncov-parser --upgrade
        pip install git+https://github.com/jts/ncov-watch.git --upgrade
        exit
    else
        shift
    fi
done

# CHECK INPUTS #
################

if [ -d "$FASTQ_PAIRS" ]; then
    echo "Directory '$FASTQ_PAIRS' exists"
else
    if [ $FASTQ_PAIRS = 0 ]; then
        echo "ERROR: Please input a paired fastq directory with '-d'"
        echo "$HELP"
        exit 1
    else
        echo "ERROR: Directory '$FASTQ_PAIRS' does not exist"
        echo "Please input a valid paired fastq directory"
        exit 1
    fi
fi

if containsElement "$PRIMER_SCHEME" "${schemeArray[@]}"; then
    echo "Using primer scheme $PRIMER_SCHEME"
elif [[ $BED != 0 ]] && [[ $AMPLICON != 0 ]]; then
    echo "Using custom amplicons"
    PATHBED=$(realpath $BED)
    PATHAMPLICON=$(realpath $AMPLICON)
    CUSTOM=true
else
    # Check for custom scheme
    if [[ $BED != 0 ]] && [[ $AMPLICON = 0 ]]; then
        echo "ERROR: Passed '--bed' without also passing '--amplicon'."
        echo "       Re-run the command with the missing argument to continue"
        exit 1
    elif [[ $BED = 0 ]] && [[ $AMPLICON != 0 ]]; then
        echo "ERROR: Passed '--amplicon' without also passing '--bed'"
        echo "       Re-run the command with the missing argument to continue"
        exit 1
    fi

    # Check if no scheme given or if scheme name was wrong
    if [ $PRIMER_SCHEME = 0 ]; then
        echo "ERROR: Please specify a primer scheme"
        echo "Primer schemes available are ${schemeArray[@]}"
        exit 1
    else
        echo "ERROR: $PRIMER_SCHEME unavailable"
        echo "Primer schemes available are ${schemeArray[@]}"
        exit 1
    fi
fi

if [[ $CORES == +([0-9]) ]]; then
    echo "Using $CORES cores for analysis"
else
    echo "ERROR: Cores input (-c) not an integer"
    exit 1
fi

if [ $METADATA_TSV = 0 ]; then
    :
elif [ -f $METADATA_TSV ]; then
    echo "$METADATA_TSV file found, using it"
    FULL_PATH_METADATA=$(realpath $METADATA_TSV)
else
    echo "ERROR: Metadata input $METADATA_TSV was not found"
    exit 1
fi

# Envs needed
envArray=($SIGNAL_ENV $NCOV_ENV 'snpdist_signal' 'qc_summary_signal' 'type_variants_signal')
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

# Make mamba env or activate it
if [[ $(conda env list | awk '{print $1}' | grep "^mamba-signal"'$') ]]; then
    conda activate mamba-signal
else
    echo "Generating 'mamba-signal' environment to install needed dependencies"
    conda create --yes -n mamba-signal -c conda-forge -c defaults mamba
    conda activate mamba-signal
fi

# Conda Env Check #
for ENV in ${envArray[@]}; do
    # If user given, don't add BASE_ENV_PATH to it
    if [[ "$USER_SIGNAL" = true && "$ENV" = "$SIGNAL_ENV" ]] || [[ "$USER_NCOV" = true && "$ENV" = "$NCOV_ENV" ]]; then
        FULL_ENV="$ENV"
    else
        FULL_ENV="$BASE_ENV_PATH/$ENV"
    fi

    # Check if env exists in user envs. NOTE if it is a env not listed in `conda env list` it will error out
    if [[ $(conda env list | awk '{print $1}' | grep "^$FULL_ENV"'$') ]]; then
        echo "$FULL_ENV found"
    else
        echo "Conda env '$FULL_ENV' wasn't found in the env list. Attempting to make the environment in $SCRIPTPATH/.snakemake"
        mamba env create -f=$SCRIPTPATH/conda_envs/$ENV.yaml -p $FULL_ENV
    fi
done
# Deactivate mamba, not needed anymore once envs made
conda deactivate
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

# Metadata #
if [ $METADATA_TSV = 0 ]; then
    sed -i -e 's/^metadata/#metadata/' ./ncov-tools/config.yaml
    cleanup_metadata=""
else
    # Check if metadata has correct ncov-tools columns, if not we only append it to final output csv
    if $(head -n 1 $FULL_PATH_METADATA | grep -q ct) && $(head -n 1 $FULL_PATH_METADATA | grep -q date); then
        echo "Metadata contains correct ncov-tools headers"
    else
        echo "Metadata is missing ncov-tools headers, appending it only to final QC output"
        sed -i -e 's/^metadata/#metadata/' ./ncov-tools/config.yaml
    fi
    cp $FULL_PATH_METADATA ./ncov-tools/metadata.tsv
    cleanup_metadata="--sample_sheet $FULL_PATH_METADATA"
fi

# Echo in correct data locations to parameters config
# Done as for the primer pairs to work we need the full path
if [ "$CUSTOM" = true ]; then
    # SIGNAL Parameters
    echo "scheme_bed: '$PATHBED'" >> $SIGNAL_CONFIG
    echo "amplicon_loc_bed: '$PATHAMPLICON'" >> $SIGNAL_CONFIG
    echo "primer_pairs_tsv: ''" >> $SIGNAL_CONFIG

    # NCOV-TOOLS Parameters
    echo "amplicon_bed: $PATHAMPLICON" >> ./ncov-tools/config.yaml
    echo "primer_bed: $PATHBED" >> ./ncov-tools/config.yaml

elif [ "$PRIMER_SCHEME" == "articV3" ]; then
    # SIGNAL Parameters
    echo "scheme_bed: 'resources/primer_schemes/artic_v3/nCoV-2019.primer.bed'" >> $SIGNAL_CONFIG
    echo "amplicon_loc_bed: 'resources/primer_schemes/artic_v3/ncov-qc_V3.scheme.bed'" >> $SIGNAL_CONFIG
    echo "primer_pairs_tsv: '-f $SCRIPTPATH/resources/primer_pairs/articV3_primer_pairs.tsv'" >> $SIGNAL_CONFIG

    # NCOV-TOOLS Parameters
    echo "amplicon_bed: $SCRIPTPATH/resources/primer_schemes/artic_v3/ncov-qc_V3.scheme.bed" >> ./ncov-tools/config.yaml
    echo "primer_bed: $SCRIPTPATH/resources/primer_schemes/artic_v3/nCoV-2019.bed" >> ./ncov-tools/config.yaml

elif [ "$PRIMER_SCHEME" == "articV4" ]; then
    # SIGNAL Parameters
    echo "scheme_bed: 'resources/primer_schemes/artic_v4/nCoV-2019.primer.bed'" >> $SIGNAL_CONFIG
    echo "amplicon_loc_bed: 'resources/primer_schemes/artic_v4/ncov-qc_V4.scheme.bed'" >> $SIGNAL_CONFIG
    echo "primer_pairs_tsv: '-f $SCRIPTPATH/resources/primer_pairs/articV4_primer_pairs.tsv'" >> $SIGNAL_CONFIG

    # NCOV-TOOLS Parameters
    echo "amplicon_bed: $SCRIPTPATH/resources/primer_schemes/artic_v4/ncov-qc_V4.scheme.bed" >> ./ncov-tools/config.yaml
    echo "primer_bed: $SCRIPTPATH/resources/primer_schemes/artic_v4/nCoV-2019.bed" >> ./ncov-tools/config.yaml

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
    echo "amplicon_loc_bed: 'resources/primer_schemes/2kb_resende/ncov-qc_resende.scheme.bed'" >> $SIGNAL_CONFIG
    echo "primer_pairs_tsv: '-f $SCRIPTPATH/resources/primer_pairs/resende_primer_pairs.tsv'" >> $SIGNAL_CONFIG

    # NCOV-TOOLS Parameters
    echo "amplicon_bed: $SCRIPTPATH/resources/primer_schemes/2kb_resende/ncov-qc_resende.scheme.bed" >> ./ncov-tools/config.yaml
    echo "primer_bed: $SCRIPTPATH/resources/primer_schemes/2kb_resende/nCoV-2019.bed" >> ./ncov-tools/config.yaml

elif [ "$PRIMER_SCHEME" == "V2resende" ]; then
    # SIGNAL Parameters
    echo "scheme_bed: 'resources/primer_schemes/2kb_resende_v2/nCoV-2019.primer.bed'" >> $SIGNAL_CONFIG
    echo "amplicon_loc_bed: 'resources/primer_schemes/2kb_resende_v2/ncov-qc_resende.scheme.bed'" >> $SIGNAL_CONFIG
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
if [ "$USER_SIGNAL" = true ]; then
    conda activate $SIGNAL_ENV
else
    conda activate $BASE_ENV_PATH/$SIGNAL_ENV
fi

# Echo out relevant info to double check it looks ok
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
if [ "$USER_NCOV" = true ]; then
    FULL_NCOV_ENV="$NCOV_ENV"
else
    FULL_NCOV_ENV="$BASE_ENV_PATH/$NCOV_ENV"
fi
bash nml_automation/run_ncov-tools.sh $CORES $SCRIPTPATH $FULL_NCOV_ENV $PDF


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
conda activate $BASE_ENV_PATH/qc_summary_signal

for RESULT in $(ls -d signal_results/*/)
do
    # Getting and keeping track of the name
    found_name="$(basename $RESULT)"
    relative_sample_path="$RESULT/core/"
    snpeff="./ncov-tools/qc_annotation/${found_name}_aa_table.tsv"
    echo "Processing $found_name"

    # gzip variants vcf from freebayes to allow it to go into summary pipeline
    gzip $RESULT/freebayes/${found_name}.variants.norm.vcf

    # Run the python summary
    python ./nml_automation/qc.py --illumina \
        --outfile ./summary_csvs/${found_name}.qc.csv \
        --sample ${found_name} \
        --ref ${ref} \
        --bam ${relative_sample_path}${found_name}_viral_reference.mapping.primertrimmed.sorted.bam \
        --fasta ${relative_sample_path}${found_name}.consensus.fa \
        --vcf $RESULT/freebayes/${found_name}.variants.norm.vcf.gz \
        --pangolin ${pangolin} \
        --ncov_summary ${ncov_qc} \
        --ncov_negative ${ncov_neg} \
        --revision v1.5.0 \
        --script_name covid-19-signal \
        --sequencing_technology illumina \
        --scheme $PRIMER_SCHEME \
        --snpeff_tsv $snpeff \
        --pcr_bed ./resources/pcr_primers.bed \
        $cleanup_metadata
done

# Concatenate all summaries
csvtk concat ./summary_csvs/*.qc.csv > ./summary_csvs/combined.csv
# Check negative controls
python ./nml_automation/negative_control_fixes.py --qc_csv ./summary_csvs/combined.csv --output_prefix $RUNNAME

# SNP DISTS #
#############
conda deactivate
conda activate $BASE_ENV_PATH/snpdist_signal

if [ -f "./ncov-tools/qc_analysis/${RUNNAME}_aligned.fasta" ];
then
    snp-dists ./ncov-tools/qc_analysis/${RUNNAME}_aligned.fasta > matrix.tsv
else
    echo "No ./ncov-tools/qc_analysis/${RUNNAME}_aligned.fasta, skipping the snp-dist check"
fi

# TYPE VARIANTS #
#################
conda deactivate
conda activate $BASE_ENV_PATH/type_variants_signal

cat ./ncov-tools/files/*.fasta > all_seq.fasta
mafft --auto --keeplength --addfragments all_seq.fasta ./ncov-tools/nCoV-2019.reference.fasta > aligned_mafft.fasta

for CSV in $SCRIPTPATH/resources/variant_profile_csvs/*.csv
do
    CSV_NAME="$(basename $CSV)"

    type_variants.py --fasta-in aligned_mafft.fasta --variants-config $CSV --reference ./ncov-tools/nCoV-2019.reference.fasta --variants-out $CSV_NAME --append-genotypes
done
### END RUNNING ###

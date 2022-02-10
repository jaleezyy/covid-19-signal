#!/usr/bin/env bash

# Script runs under the ncov-qc conda environment through snakemake integration
# SIGNAL auto-generates filesystem so directory and configuration files are already set

set -e # exit if pipeline returns non-zero status
set -o pipefail # return value of last command to exit with non-zero status
shopt -s extglob

CORES=1
SIGNAL=0

HELP="""
This script attempts to execute ncov-tools for postprocessing of SIGNAL sequencing data.
SIGNAL ncov_tools already prepares the filesystem but we cannot directly execute a snakemake workflow within another workflow.
This script acts as a liason between python snakemake and terminal bash, executing snakemake outside the background python workflow.
Results files will be moved to SIGNAL results to consolidate all needed files. 


First argument provided for this script is '-c' for the number of cores. Default value is 1.
Second argument is '-s' for SIGNAL results directory (namely the name) - much of the script operates on relative pathing

"""

while getopts ":c:s:" option; do
	case "${option}" in
		c) CORES=$OPTARG;; # number of cores
		s) SIGNAL=$OPTARG;; # SIGNAL results directory
	esac
done

if [ $1 = 'help' ]; then
	echo "$HELP"
	exit 1
fi

if [ $SIGNAL = 0 ] ; then
    echo "You must specify the name of the directory holding SIGNAL results."
    echo "$HELP"
    exit 1
fi

# Start point for executing from ncov-tools.py is SIGNAL results directory
RESULTS=$PWD

# determine where ncov-tools is located *add check*
# change directory to ncov-tools
cd ../ncov-tools

# run ncov-tools
snakemake -s workflow/Snakefile --cores ${CORES} all

# move ncovresults to SIGNAL results directory
mv ${SIGNAL}'_ncovresults' ${RESULTS}/ncov-tools-results

# return success
exit 0
# PDF summary to be set as future parameter ('all' is required)

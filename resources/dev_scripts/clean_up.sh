#!/bin/bash

set -e # exit if pipeline returns non-zero status
set -o pipefail # return value of last command to exit with non-zero status

sample_table=0
output_dir="output"

HELP="""
Clean up SIGNAL output directories and files in to one consolidated directory. Uses sample table to derive sample directories.

Flags:
    -s  :  Path to sample table
    -n  :  Name of final output directory (default 'output')
"""

while getopts ":s:n:" option; do
    case "${option}" in
        s) sample_table=$OPTARG;;
        n) output_dir=$OPTARG;;
    esac
done

if [ $sample_table = 0 ] ; then
    echo "You must specify the sample table used for SIGNAL analysis."
    echo "$HELP"
    exit 1
fi

if [ -d $output_dir ]; then
	echo "Directory ${output_dir} already exists."
	echo "$HELP"
	exit 1
else
	echo -e "Making output directory called ${output_dir}\n"
	mkdir $output_dir
fi

for file in $(sed '1d' ${sample_table} | cut -d, -f1); do
	if [ -d $file ] ; then
		echo $file
		mv $file $output_dir
	fi
done

if [ -f "summary.html" ] ; then
	echo "summary.html"
	mv summary.html $output_dir
fi

if [ -f "summary.zip" ] ; then
	echo "summary.zip"
	mv summary.zip $output_dir
fi

if [ -f "summary1.png" ] ; then
	echo "summary1.png"
	mv summary1.png $output_dir
fi

if [ -f "summary2.png" ] ; then
	echo "summary2.png" 
	mv summary2.png $output_dir
fi

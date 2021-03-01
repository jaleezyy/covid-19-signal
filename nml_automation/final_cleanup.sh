#!/bin/bash

mkdir -p ./summary_csvs

eval "$(conda shell.bash hook)"

conda activate <placeholder>

# Setting up relative paths from the root dir to needed data
ncov_qc="./ncov-tools/qc_reports/nml_summary_qc.tsv"
ncov_neg="./ncov-tools/qc_reports/nml_negative_control_report.tsv"
ref="./data/MN908947.3.fasta"
pangolin="./ncov-tools/lineages/nml_lineage_report.csv"

# Other inputs
# Revision can be found on the signal github page for what the code is, **manual update is needed below**
# https://github.com/jaleezyy/covid-19-signal
rev="Release-v1.1.0"
pcr_primer_bed="./data/pcr_primers.bed"

# Check to see if we have a negative control. If so do nothing. If not make a "fake" one
# Need a fake one for the qc.py to not error out
if [ -f "$ncov_neg" ];
then
    echo "Negative contol data found" >> ./summary_csvs/process.log
else
    touch ./ncov-tools/qc_reports/nml_negative_control_report.tsv
fi

# Using the folder names (as each sample is in its sample-name named folder) we can do the qc checks
# Not making too many instances of this as it doesn't take too long per and I don't want to have a harder time
# linking up all of the processes after
# If we are being bottlenecked here, can make each its own sbatch and put the IDs to a file for the next step to run off of!
for i in $(ls -d signal_results/*/)
do
    # Getting and keeping track of the name
    name="$(basename $i)"
    relative_sample_path="$i/core/"
    echo "Processing $name" >> ./summary_csvs/process.log

    # Run the python summary
    python ./nml_automation/qc.py --illumina --outfile ./summary_csvs/${name}.qc.csv --sample ${name} --ref ${ref} --bam ${relative_sample_path}${name}_viral_reference.mapping.primertrimmed.sorted.bam --fasta ${relative_sample_path}${name}.consensus.fa --variants ${relative_sample_path}${name}_ivar_variants.tsv --pangolin ${pangolin} --ncov_summary ${ncov_qc} --ncov_negative ${ncov_neg} --revision ${rev} --pcr_bed ${pcr_primer_bed}

done

echo "Done all data processing" >> ./summary_csvs/process.log

# Small annoying png files made (as using the connor lab script)
# Thus move them to somewhere no one will look :p
mkdir -p qc_depth
mv *.png qc_depth
echo "Moved QC depth plots" >> ./summary_csvs/process.log

# Combine all outputs
csvtk concat ./summary_csvs/*.qc.csv > ./summary_csvs/combined.csv

# Finally, correct for the negative control(s)!
python ./nml_automation/negative_control_fixes.py --qc_csv ./summary_csvs/combined.csv --output_prefix nml

echo "Done negative control fixes" >> ./summary_csvs/process.log
echo "\nDone!!\n" >> ./summary_csvs/process.log

# snp-dist
#conda deactivate
#conda activate <placeholder>

if [ -f "./ncov-tools/qc_analysis/nml_aligned.fasta" ];
then
    snp-dists ./ncov-tools/qc_analysis/nml_aligned.fasta > matrix.tsv
else
    echo 'No ./ncov-tools/qc_analysis/nml_aligned.fasta, skipping the snp-dist check'
fi

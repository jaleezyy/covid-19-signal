#!/bin/bash

set -e # exit if pipeline returns non-zero status
set -o pipefail # return value of last command to exit with non-zero status

DATABASE=0
OUTPUT=0
FORWARD=0
REVERSE=0
THREADS=1
#THREADS=$(nproc) # set maximum number of processors

HELP="""

Flags:

	-d	:	FASTA file to generate database to be used for mapping with HiSAT2
	-1	:	Forward reads (paired)
	-2	:	Reverse reads (paired)
	-o	:	Output directory (& filename) - directory will have '_align' appended
	-t	:	Number of threads to use for HiSAT2. Will default to 1.

"""

while getopts ":d:o:h:1:2:t:" option; do
	case "${option}" in
#		f) SOURCE=$OPTARG;; # FASTQ
    d) DATABASE=$OPTARG;;
    o) OUTPUT=$OPTARG;;
    1) FORWARD=$OPTARG;;
    2) REVERSE=$OPTARG;;
    t) THREADS=$OPTARG
		
esac
done

if [ $DATABASE = 0 ] ; then
	echo "You must specify a FASTA file to generate the database file."
  echo "$HELP"
	exit 1
fi


if [ $FORWARD = 0 ] || [ $REVERSE = 0 ] ; then
  echo "You must provide both forward and reverse reads (i.e. R1 and R2)."
  echo "$HELP"
  exit 1
fi

if [ $OUTPUT = 0 ] ; then
  echo "Please provide a name for the final output files/directories."
  echo "$HELP"
  exit 1
fi 

# Generate HiSAT2 database
mkdir ${OUTPUT}_tmp
hisat2-build $DATABASE ${OUTPUT}_tmp/ref_align

# Make final output directory
mkdir ${OUTPUT}_align

# HiSAT2 alignment
echo "Running HiSAT2 alignment."
hisat2 -x ${OUTPUT}_tmp/ref_align -1 $FORWARD -2 $REVERSE --threads $THREADS --summary-file ${OUTPUT}_align/${OUTPUT}_alignment_summary.log -S ${OUTPUT}_align/align.sam 
echo "Alignment complete (log found as ${OUTPUT}_alignment_summary.log). Processing now."
echo -e "\nSummary Statistics:\n" >> ${OUTPUT}_align/${OUTPUT}_alignment_summary.log

# SAM -> BAM process
samtools view -b ${OUTPUT}_align/align.sam | samtools sort > ${OUTPUT}_align/${OUTPUT}.bam
samtools index ${OUTPUT}_align/${OUTPUT}.bam

# Calculate average depth of coverage, median, minimum, and maximum
bedtools genomecov -d -ibam ${OUTPUT}_align/${OUTPUT}.bam | awk '{sum+=$3} END {print "Average Depth of Coverage to '$(basename ${DATABASE})' =",sum/NR"x"}' | tee -a ${OUTPUT}_align/${OUTPUT}_alignment_summary.log
bedtools genomecov -d -ibam ${OUTPUT}_align/${OUTPUT}.bam | awk '{print $3}' | sort -n | awk '{a[NR]=$0} END {print "Median Coverage =",((NR%2==1)?a[int(NR/2)+1]:(a[NR/2]+a[NR/2+1])/2)"x"}' | tee -a ${OUTPUT}_align/${OUTPUT}_alignment_summary.log
bedtools genomecov -d -ibam ${OUTPUT}_align/${OUTPUT}.bam | awk '{print $3}' | echo 'Minimum Coverage = '$(jq -s min)'x' | tee -a ${OUTPUT}_align/${OUTPUT}_alignment_summary.log
bedtools genomecov -d -ibam ${OUTPUT}_align/${OUTPUT}.bam | awk '{print $3}' | echo 'Maximum Coverage = '$(jq -s max)'x' | tee -a ${OUTPUT}_align/${OUTPUT}_alignment_summary.log

echo "Cleaning up."

rm ${OUTPUT}_align/align.sam
rm -rf ${OUTPUT}_tmp

echo "Done."

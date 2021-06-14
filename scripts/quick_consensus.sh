#!/bin/bash

set -euo pipefail

function usage() {
	echo "Usage: quick_consensus.sh [--help|-h] [--reference MN908947.3.fna] [--] vcf_files..." 1>&2
    exit 1
}

ref=''
while [[ $# -gt 0 ]]; do
    case $1 in
        --reference) shift; ref="$1" ;;
        --help|-h) usage ;;
        --) shift; break ;;
        -*) echo "Unknown flag '$1'" 1>&2; usage ;;
        *) break ;;
    esac
    shift
done

if [ -z "$ref" ] ; then
    echo "You must provide a reference genome"
    usage
    exit 1
fi

for vcf in "$@"; do 
    sample_name=$(echo "$(basename $vcf)" | cut -d '.' -f1 | cut -d '_' -f1)
    sample_dir="$(dirname $vcf)"
    bgzip --force --force "$vcf"
    bcftools index "$vcf".gz 
    bcftools consensus --haplotype R --fasta-ref "$ref" "$vcf".gz | sed "s/>.*/>$sample_name/"  > "${sample_dir}/${sample_name}".vcf_consensus.fna
done

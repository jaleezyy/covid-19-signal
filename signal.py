#!/bin/bash

import argparse
import subprocess
import re
from pathlib import Path

def check_directory(path: str) -> Path:
    """
    Check an input directory exists and is readable
    """
    path = Path(path)
    if path.exists() and path.is_dir():
        return path
    else:
        raise argparse.ArgumentTypeError(f"{path} can't be read")

def check_files_and_resolve_paths(project_directory):
    """
    Check read files and resolve paths
    """
    sample_data = []
    for r1 in project_directory.rglob("*_L001_R1_001.fastq.gz"):
        r1 = r1.resolve()
        r2 = Path(str(r1).replace("_L001_R1_001.fastq.gz", "_L001_R2_001.fastq.gz"))
        if r2.exists():
            r2 = r2.resolve()
        else:
            raise ValueError(f"R2 cannot be found {r2}")
        sample_name = "_".join(r1.name.replace("_L001_R1_001.fastq.gz", '').split("_")[:-1])
        sample_data.append([str(sample_name), str(r1), str(r2)])

    print(f"{len(sample_data)} samples detected")
    return sample_data


def write_sample_table(sample_data, output_table):
    """
    Generate sample table from parsed folder read data
    """
    header = ",".join(['sample', 'r1_path', 'r2_path']) + '\n'
    with open(output_table, 'w') as out_fh:
        out_fh.write(header)
        for sample in sample_data:
            out_fh.write(",".join(sample) + '\n')

def write_config_file(run_name, config_file):

	config = f"""amplicon_loc_bed: resources/primer_schemes/artic_v4/SARS-CoV-2.scheme.bed
breseq_reference: data/MN908947.3.gbk
composite_reference: data/composite_human_viral_reference.fna
kraken2_db: data/Kraken2/db
min_len: 20
min_qual: 20
mpileup_depth: 100000
negative_control_prefix: ['']
phylo_include_seqs: data/extra_data.fasta
primer_pairs_tsv: ''
result_dir: {run_name}_results_dir
run_breseq: False
run_freebayes: True
samples: {run_name}_sample_table.csv
scheme_bed: resources/primer_schemes/artic_v4/SARS-CoV-2.primer.bed
var_freq_threshold: 0.75
var_min_coverage_depth: 10
var_min_freq_threshold: 0.25
var_min_variant_quality: 20
viral_reference_contig_name: MN908947.3
viral_reference_feature_coords: data/MN908947.3.gff3
viral_reference_genome: data/MN908947.3.fasta"""

	with open(config_file, 'w') as fh:
		fh.write(config)


if __name__ == '__main__':

    parser = argparse.ArgumentParser("Generate sample sample and check reads")
    parser.add_argument('-d', '--directory', required=True, type=check_directory,
                        help="Path to directory containing reads")
    args = parser.parse_args()

    run_name = args.directory.name

    sample_data = check_files_and_resolve_paths(args.directory)

    write_sample_table(sample_data, run_name + "_sample_table.csv")#args.output)

    config_file = run_name + "_config.yaml"
    write_config_file(run_name, config_file)

    subprocess.run(f"snakemake --conda-frontend mamba --configfile {config_file} --cores=12 --use-conda --conda-prefix=$PWD/.snakemake/conda all -kp", shell=True, check=True)

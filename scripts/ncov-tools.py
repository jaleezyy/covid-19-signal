#!/usr/bin/env python

import os
import shutil
import subprocess
import fileinput

def set_up():
    print("Writing config for ncov to ncov-tools/config.yaml")

    data_root = os.path.abspath('ncov-tools/data')
    if os.path.exists(data_root):
        shutil.rmtree(data_root)
    os.mkdir(data_root)

    config = {'data_root': f"'{data_root}'",
              'amplicon_bed': f"'{os.path.abspath(snakemake.config['scheme_bed'])}'", #grab from signal snakemake config
              'reference_genome': f"'{os.path.abspath(snakemake.config['viral_reference_genome'])}'", #grab from signal snakemake config
              'bam_pattern': "'{data_root}/{sample}.bam'", # symlink files following this
              'consensus_pattern': "'{data_root}/{sample}.consensus.fasta'", # symlink files following this
              'tree_include_consensus': f"'{os.path.abspath(snakemake.config['phylo_include_seqs'])}'",
              #'metadata': snakemake.config['samples'], # causes error get example
              'assigned_lineages': 'true'}

    with open('ncov-tools/config.yaml', 'w') as fh:
        for key, value in config.items():
            fh.write(f"{key}: {value}\n")

    print("Linking files to ncov")
    for bam in snakemake.input['bams']:
        sample = bam.split('/')[0]
        ln_path = f"{data_root}/{sample}.bam"
        if not os.path.exists(ln_path):
            os.link(bam, ln_path)

    for consensus in snakemake.input['consensus']:
        sample = consensus.split('/')[0]
        ln_path = f"{data_root}/{sample}.consensus.fasta"
        if not os.path.exists(ln_path):
            os.link(consensus, ln_path)
        for line in fileinput.input(ln_path, inplace=True):
            if line.startswith(">"):
                new_header = str(">"+sample)
                new_line = line.replace(line, new_header)
                print(new_line, end='\n')
            else:
                print(line, end='\n')

    os.chdir('ncov-tools')

def run_all_qc_sequencing():
    subprocess.run(f"snakemake -s qc/Snakefile --cores {snakemake.threads} all_qc_sequencing",
                    shell=True)

def run_all_qc_analysis():
    os.system(f"snakemake -s qc/Snakefile --cores {snakemake.threads} all_qc_analysis")

def run_all_qc_summary():
    os.system(f"snakemake -s qc/Snakefile --cores {snakemake.threads} all_qc_summary")

if __name__ == '__main__':

    set_up()
    run_all_qc_sequencing()
    run_all_qc_analysis()
    # seems to include hard-coded issues
    #run_all_qc_summary()

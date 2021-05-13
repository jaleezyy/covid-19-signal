#!/usr/bin/env python

import argparse
import subprocess
from pathlib import Path
import time
import pandas as pd
import shutil


def check_file(path: str) -> Path:
    """
    Check an input file exists and is readable
    """
    path = Path(path)
    if path.exists() and path.is_file():
        return path
    else:
        raise argparse.ArgumentTypeError(f"{path} can't be read")


def update_pangolin():
    """
    Ensure pangolin is updated to the latest release
    """
    subprocess.check_output(["pangolin", "--update"])


def update_nextclade():
    """
    Ensure nextclade is updated to the latest release
    """
    subprocess.check_output(["npm", "install", "-g", "@neherlab/nextclade"])


def run_nextclade(input_genomes, threads):
    """
    Execute nextclade and collect assignments
    """
    output_file = Path(f"nextclade_temp_{time.time()}.csv")
    subprocess.check_output(f"nextclade -i {input_genomes} -j {threads} "
                            f"-c {str(output_file)}".split(),
                            stderr=subprocess.DEVNULL)
    if not output_file.exists():
        raise FileNotFoundError(f"{str(output_file)} not created, check "
                                 "nextclade install")
    nextclade_df = pd.read_csv(str(output_file), sep=";")

    # get version information
    nextclade_version = subprocess.run(f"nextclade --version".split(),
                                       stdout=subprocess.PIPE)
    nextclade_version = nextclade_version.stdout.decode('utf-8').strip()
    nextclade_df['nextclade_version'] = f"nextclade {nextclade_version}"

    # tidy up dataframe
    nextclade_df = nextclade_df.rename(columns={'seqName': 'isolate',
                                          'clade': 'nextstrain_clade',
                                          'qc.overallStatus': 'nextclade_qc',
                                          'errors': 'nextclade_errors'})
    nextclade_df = nextclade_df.drop([qc_col for qc_col in \
                                        nextclade_df.columns \
                                        if qc_col.startswith('qc.')], axis=1)

    output_file.unlink()

    return nextclade_df


def run_pangolin(input_genomes, threads):
    """
    Execute pangolin and collect assignments
    """
    output_dir = Path(f"pangolin_tmp_{time.time()}")
    subprocess.check_output(f"pangolin {input_genomes} -t {threads} "
                            f"-o {str(output_dir)}".split(),
                            stderr=subprocess.DEVNULL)

    output_path = output_dir / "lineage_report.csv"
    if not output_path.exists():
        raise FileNotFoundError(f"{str(output_path)} not created, check "
                                 "pangolin install")

    pangolin_df = pd.read_csv(str(output_path), sep=',')

    # get version information for pangolin, pangolearn model, and lineages
    pangolearn_version = subprocess.run("pangolin --pangoLEARN-version".split(),
                                        stdout=subprocess.PIPE)
    pangolearn_version = pangolearn_version.stdout.decode('utf-8').strip()
    pangolin_df['pangoLEARN_version'] = pangolearn_version

    pangolin_version = subprocess.run("pangolin --version".split(),
                                      stdout=subprocess.PIPE)
    pangolin_version = pangolin_version.stdout.decode('utf-8').strip()
    pangolin_df['pangolin_version'] = pangolin_version

    # tidy up the dataframe
    if 'probability' in pangolin_df:
        pangolin_df = pangolin_df.rename(columns={'taxon': 'isolate',
                                              'lineage': 'pangolin_lineage',
                                              'status': 'pangolin_qc',
                                              'note': 'pangolin_note',
                                              'probability': 'pangolin_lineage_score'})
    elif 'conflict' in pangolin_df:
        pangolin_df = pangolin_df.rename(columns={'taxon': 'isolate',
                                              'lineage': 'pangolin_lineage',
                                              'status': 'pangolin_qc',
                                              'note': 'pangolin_note',
                                              'conflict': 'pangolin_lineage_score'})


    # remove temp output
    shutil.rmtree(output_dir)

    return pangolin_df


def collate_output(nextclade, pangolin, output):
    """
    Merge and tidy lineage assignments
    """
    merged_df = pangolin.merge(nextclade, on='isolate', how='outer')

    merged_df = merged_df[['isolate', 'pangolin_lineage',
                           'pangolin_lineage_score', 'pangolin_note',
                           'pangolin_qc', 'nextstrain_clade',
                           'nextclade_qc', 'nextclade_errors',
                           'totalGaps', 'totalInsertions', 'totalMissing',
                           'totalMutations', 'totalNonACGTNs',
                           'totalPcrPrimerChanges',
                           'substitutions', 'deletions', 'insertions',
                           'missing', 'nonACGTNs',
                           'pcrPrimerChanges', 'aaSubstitutions',
                           'totalAminoacidSubstitutions',
                           'aaDeletions', 'totalAminoacidDeletions',
                           'alignmentStart', 'alignmentEnd', 'alignmentScore',
                           'pangolin_version',
                           'pangoLEARN_version', 'nextclade_version']]
    merged_df.to_csv(output, sep='\t', index=False)


if __name__ == '__main__':

    parser = argparse.ArgumentParser('Assign pangolin and nextstrain lienages '
                                     'for a set of SARS-CoV-2 genomes')
    parser.add_argument("-i", "--input_genomes", type=check_file, required=True,
                        help="Concatenated fasta containing consensus genomes")
    parser.add_argument("-t", "--threads", default=8, type=int,
                        help="Number of threads for pangolin/nextclade")
    parser.add_argument("-o", "--output", default="lineage_assignments.tsv",
                        help="Output file for collated assignment table")
    args = parser.parse_args()

    update_pangolin()
    pangolin = run_pangolin(args.input_genomes, args.threads)

    update_nextclade()
    nextclade = run_nextclade(args.input_genomes, args.threads)

    collate_output(nextclade, pangolin, args.output)

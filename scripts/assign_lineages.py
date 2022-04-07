#!/usr/bin/env python

import argparse
import subprocess
from pathlib import Path
import time
import pandas as pd
import shutil
import os, sys


def check_file(path: str) -> Path:
    """
    Check an input file exists and is readable
    """
    path = Path(path)
    if path.exists() and path.is_file():
        return path
    else:
        raise argparse.ArgumentTypeError(f"{path} can't be read")

def update_latest_pangolin():
    """
    Ensure pangolin is updated to the latest release
    """
    try:
        subprocess.check_output(["pangolin", "--update"])
    except subprocess.CalledProcessError:
        print("Something went wrong updating Pangolin! No changes were made! Attempting to proceed...")

def update_pangolin(vers):
    """
    Update pangolin to a specific version
    """
    script_dir = os.path.dirname(sys.argv[0])
    script = os.path.join(script_dir, "pangolin_specific_version_update.py")
    subprocess.run([script, '--versions_file', vers])

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
    
   # check final versions for pangolin
    subprocess.check_output(["pangolin", "--all-versions"])
    
    output_dir = Path(f"pangolin_tmp_{time.time()}")
    subprocess.check_output(f"pangolin {input_genomes} -t {threads} "
                            f"-o {str(output_dir)}".split(),
                            stderr=subprocess.DEVNULL)

    output_path = output_dir / "lineage_report.csv"
    if not output_path.exists():
        raise FileNotFoundError(f"{str(output_path)} not created, check "
                                 "pangolin install")

    pangolin_df = pd.read_csv(str(output_path), sep=',')

    # tidy up the dataframe
    pangolin_df = pangolin_df.rename(columns={'taxon': 'isolate',
                                              'lineage': 'pango_lineage',
                                              'status': 'pangolin_qc',
                                              'note': 'pangolin_note',
                                              'conflict': 'pangolin_conflict',
                                              'ambiguity_score': 'pangolin_ambiguity_score'})

    # remove temp output
    shutil.rmtree(output_dir)

    return pangolin_df


def collate_output(nextclade, pangolin, output):
    """
    Merge and tidy lineage assignments
    """
    merged_df = pangolin.merge(nextclade, on='isolate', how='outer')

    try:
        merged_df = merged_df[['isolate', 'pango_lineage',
                               'pangolin_conflict', 'pangolin_ambiguity_score',
                               'pangolin_note', 'scorpio_call', 'scorpio_support',
                               'scorpio_conflict',
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
                               'pangolin_version', 'pango_version',
                               'pangoLEARN_version', 'pango_version', 'nextclade_version']]
    except KeyError: # adjust for pangolin v4
        merged_df = merged_df[['isolate', 'pango_lineage',
                               'pangolin_conflict', 'pangolin_ambiguity_score',
                               'pangolin_note', 'scorpio_call', 'scorpio_support',
                               'scorpio_conflict',
                               'nextstrain_clade',
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
                               'version','pangolin_version', 'nextclade_version']]
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
    parser.add_argument("-p", "--pangolin_ver", type=check_file, required=False, default=None,
                        help="Input file containing version information for PANGOLIN tools")
    parser.add_argument("--skip", action="store_true", help="Skip updates to pangolin and nextclade")
    args = parser.parse_args()

    if args.skip is False:
        if args.pangolin_ver is None: 
            update_latest_pangolin()
        else:
            update_pangolin(args.pangolin_ver)
        update_nextclade()

    pangolin = run_pangolin(args.input_genomes, args.threads)
    nextclade = run_nextclade(args.input_genomes, args.threads)

    collate_output(nextclade, pangolin, args.output)

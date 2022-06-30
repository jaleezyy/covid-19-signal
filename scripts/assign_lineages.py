#!/usr/bin/env python

import argparse
import subprocess
from pathlib import Path
import time
import pandas as pd
import shutil
import os, sys
from datetime import datetime


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
	if vers is None:
		return None
	script_dir = os.path.dirname(sys.argv[0])
	script = os.path.join(script_dir, "pangolin_specific_version_update.py")
	subprocess.run([script, '--versions_file', vers])

def update_nextclade():
	"""
	DEPRECIATED!
	Ensure nextclade is updated to the latest release
	"""
	subprocess.check_output(["npm", "install", "-g", "@neherlab/nextclade"])

def update_nextclade_dataset(vers, skip):
	"""
	Ensure nextclade dataset is updated to the latest dataset, placed within scripts.
	Reference accession will be set by params.accession (viral_reference_contig_name).
	"""
	output_dir = os.path.join(os.path.dirname(sys.argv[0]), 'nextclade')
	nextclade_version = subprocess.run(f"nextclade --version".split(),
						stdout=subprocess.PIPE).stdout.decode('utf-8').strip().lower()
	if nextclade_version.startswith("nextclade"):
		nextclade_version = nextclade_version.split()[1]

	if skip or (vers is None):
		return output_dir, nextclade_version
	if not os.path.exists(output_dir):
		os.mkdir(output_dir)
	with open(vers) as fh:
		line = fh.readlines()
		assert len(line) == 3 # should only be nextclade, nextclade-data, and recomb
		software_ver = str(line[0].split(":", 1)[1]).strip().lstrip('v')
		requested_ver = str(line[1].split(":", 1)[1]).strip()
		recomb = eval(str(line[2].split(":")[1]).strip())

	# check current version of nextclade, if failed, we stick with the latest
	# search conda for latest version, this will be the default
	# PackagesNotFoundError for invalid version
	print("\n## Existing nextclade install:")
	print("Nextclade: " + nextclade_version + "\n")
	print("## Changing installed versions as needed:")
	try:
		if software_ver != "None": # specific version requested, check if available
			try:
				softrequest = subprocess.check_output(f"conda search -c bioconda -f nextclade={software_ver}".split()).decode('utf-8').strip().split()[-3]
				# check if already installed
				if softrequest == nextclade_version:
					print(f"Nextclade {softrequest} already installed! Skipping update!")
				else:
					print(f"Changing Nextclade from {nextclade_version} to {softrequest}!")
					subprocess.run(f"conda install -q -y -c bioconda nextclade={softrequest}", shell=True, check=True)
			except subprocess.CalledProcessError:
				print("Cannot find version requested, will ensure latest version!")
				softrequest = subprocess.check_output(f"conda search -c bioconda -f nextclade".split()).decode('utf-8').strip().split()[-3]
				# check if already installed
				if softrequest == nextclade_version:
					print(f"Nextclade {softrequest} already installed! Skipping update!")
				else:
					subprocess.run(f"conda install -q -y -c bioconda nextclade={softrequest}", shell=True, check=True)
		else:
			print(f"Installing latest version of Nextclade!")
			softrequest = subprocess.check_output(f"conda search -c bioconda -f nextclade".split()).decode('utf-8').strip().split()[-3]
			if softrequest == nextclade_version:
				print(f"Nextclade {softrequest} already installed! Skipping update!")
			else:
				subprocess.run(f"conda install -q -y -c bioconda nextclade={softrequest}", shell=True, check=True)
	except subprocess.CalledProcessError:
		print(f"Something went wrong updating Nextclade! Skipping update!")

	# check nextclade_ver, if None, assign today's date    
	try:
		if requested_ver != "None":
			# assume yyyy-mm-dd (so only yyyy and mm expected to stay ccnsistent)
			submitted = requested_ver.split("-")
			submitted_date = [s.strip() for s in submitted]
			assert len(submitted_date) == 3
			year = str(submitted_date[0])
			month = str(submitted_date[1])
			if (len(submitted_date[2].split(" ")) == 2) or (len(submitted_date[2].split("T")) == 2): # date and time part of provided tag
				if submitted_date[2].count("T") == 1: # only applies if starting input was in quotations itself in the config file
					day = str(submitted_date[2]).split("T")[0].strip()
					timestamp = str(submitted_date[2]).split("T")[1].split("+", 1)[0].split(":")
				else: 
					day = str(submitted_date[2]).split(" ")[0].strip()
					timestamp = str(submitted_date[2]).split(" ")[1].split("+", 1)[0].split(":")
				tags = [timestamp[0], timestamp[1], timestamp[2].strip("Z")]
			else: # only a date provided, assume timestamp
				day = str(submitted_date[2])
				tags = ["12", "00", "00"]
			requested = str("%s-%s-%sT%s:%s:%sZ" %(year, month, day, tags[0], tags[1], tags[2]))
		else:
			requested = None
	except (AssertionError, TypeError, ValueError): # some other input that isn't in yyyy-mm-dd date format
		print(f"\nProvided Nextclade dataset tag invalid! Downloading latest...")
		requested = None

	if recomb:
		dataset = 'sars-cov-2'
	else:
		dataset = 'sars-cov-2-no-recomb'

	# If specific tag requested, attempt to install, otherwise install latest
	accession = 'MN908947'
	if requested is not None:
		try:
			print(f"\nDownloading Nextclade {dataset} dataset tagged {requested} for reference {accession}!")
			subprocess.run(f"nextclade dataset get "
						f"--name '{dataset}' "
						f"--reference '{accession}' "
						f"--tag {requested} "
						f"--output-dir '{output_dir}'", shell=True, check=True)
		except subprocess.CalledProcessError:
			print(f"\nDatabase not found! Please check whether {requested} tag exists! Downloading latest Nextclade {dataset} dataset for reference {accession}...")
			subprocess.run(f"nextclade dataset get "
						f"--name '{dataset}' "
						f"--reference '{accession}' "
						f"--output-dir '{output_dir}'", shell=True, check=True)
	else:
		print(f"\nDownloading latest Nextclade {dataset} dataset for reference {accession}!")
		subprocess.run(f"nextclade dataset get "
					f"--name '{dataset}' "
					f"--reference '{accession}' "
					f"--output-dir '{output_dir}'", shell=True, check=True)

	# Obtain final version information for output
	nextclade_version = subprocess.run(f"nextclade --version".split(), stdout=subprocess.PIPE).stdout.decode('utf-8').strip().lower()
	if nextclade_version.startswith("nextclade"):
		nextclade_version = nextclade_version.split()[1]
	if requested is None:
		today = datetime.today().strftime('%Y-%m-%d')
		requested = f"Latest as of {today}"
	with open('final_nextclade_versions.txt', 'w+') as out:
		print("\n## Nextclade and datasets now:")
		print("Nextclade: " + nextclade_version)
		print("Reference: %s" %(accession))
		print("Dataset: %s" %(dataset))
		print("Dataset version: %s" %(requested))
		# Output to file
		print("## Nextclade and datasets now:", file=out)
		print("Nextclade: " + nextclade_version, file=out)
		print("Reference: %s" %(accession), file=out)
		print("Dataset: %s" %(dataset), file=out)
		print("Dataset version: %s" %(requested), file=out)

	return output_dir, nextclade_version


def run_nextclade(input_genomes, dataset, threads, version):
	"""
	Execute nextclade and collect assignments
	"""
	output_dir = Path(f"tmp_nextclade")
	output_file = Path(f"{output_dir}/nextclade_temp_{time.time()}.csv")
	# run altered commands for Nextclade v1 and v2+
	if version.startswith("1"):
		subprocess.check_output(f"nextclade -i {input_genomes} -j {threads} --input-dataset {dataset} -c {str(output_file)} --output-dir {output_dir}".split(), stderr=subprocess.DEVNULL)
	elif version.startswith("2"):
		subprocess.check_output(f"nextclade run -j {threads} --input-dataset {dataset} -c {str(output_file)} {input_genomes}".split(), stderr=subprocess.DEVNULL)
	else:
		print("Unknown version of nextclade!")
	if not output_file.exists():
		raise FileNotFoundError(f"{str(output_file)} not created, check "
								 "nextclade install")
	nextclade_df = pd.read_csv(str(output_file), sep=";")

	# get version information
	#nextclade_version = subprocess.run(f"nextclade --version".split(), stdout=subprocess.PIPE)
	#nextclade_version = nextclade_version.stdout.decode('utf-8').strip()
	nextclade_df['nextclade_version'] = f"nextclade {version}"

	# tidy up dataframe
	nextclade_df = nextclade_df.rename(columns={'seqName': 'isolate',
										'clade': 'nextstrain_clade',
										'qc.overallStatus': 'nextclade_qc',
										'errors': 'nextclade_errors'})
	nextclade_df = nextclade_df.drop([qc_col for qc_col in nextclade_df.columns if qc_col.startswith('qc.')], axis=1)

	# remove temp output
	output_file.unlink()
	shutil.rmtree(output_dir)

	return nextclade_df


def run_pangolin(input_genomes, threads, mode):
	"""
	Execute pangolin and collect assignments
	"""
   # check analysis mode
	if mode == "fast":
		analysis = f"--analysis-mode fast "
	else:
		analysis = f""
   # check final versions for pangolin
	subprocess.check_output(["pangolin", "--all-versions"])
	
	output_dir = Path(f"pangolin_tmp_{time.time()}")
	subprocess.check_output(f"pangolin {analysis}{input_genomes} -t {threads} "
							f"-o {str(output_dir)}".split(),
							stderr=subprocess.DEVNULL)

	output_path = output_dir / "lineage_report.csv"
	if not output_path.exists():
		raise FileNotFoundError(f"{str(output_path)} not created, check "
								"pangolin install")

	pangolin_df = pd.read_csv(str(output_path), sep=',')

	# tidy up the dataframe
	try:
		pangolin_df = pangolin_df.rename(columns={'taxon': 'isolate',
												'lineage': 'pango_lineage',
												'status': 'pangolin_qc',
												'note': 'pangolin_note',
												'conflict': 'pangolin_conflict',
												'ambiguity_score': 'pangolin_ambiguity_score'}, errors='raise')
	except KeyError: # likely associated with Pangolin v4+
		pangolin_df = pangolin_df.rename(columns={'taxon': 'isolate',
												'lineage': 'pango_lineage',
												'qc_status': 'pangolin_qc',
												'qc_notes': 'pangolin_qc_note',
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

	try: # adjust for pangolin v3 and nextclade v1.1.0
		merged_df = merged_df[['isolate', 'pango_lineage',
							'pangolin_conflict', 'pangolin_ambiguity_score',
							'pangolin_note', 'scorpio_call', 'scorpio_support',
							'scorpio_conflict',
							'pangolin_qc', 'nextstrain_clade',
							'nextclade_qc', 'nextclade_errors',
							'totalInsertions', 'totalMissing',
							'totalNonACGTNs','totalPcrPrimerChanges',
							'substitutions', 'deletions', 'insertions',
							'missing', 'nonACGTNs',
							'pcrPrimerChanges', 'aaSubstitutions',
							'totalAminoacidSubstitutions',
							'aaDeletions', 'totalAminoacidDeletions',
							'alignmentStart', 'alignmentEnd', 'alignmentScore',
							'pangolin_version', 'pango_version',
							'pangoLEARN_version', 'pango_version', 'nextclade_version']]
	except KeyError: # adjust for pangolin v4 and nextclade v1.1.0+
		merged_df = merged_df[['isolate', 'pango_lineage',
							'pangolin_conflict', 'pangolin_ambiguity_score',
							'pangolin_note', 'scorpio_call', 'scorpio_support',
							'scorpio_conflict', 'pangolin_qc', 'pangolin_qc_note',
							'nextstrain_clade',
							'nextclade_qc', 'nextclade_errors',
							'totalSubstitutions', 'totalDeletions', 'totalInsertions', 'totalMissing',
							'totalNonACGTNs', 'totalPcrPrimerChanges',
							'substitutions', 'deletions', 'insertions',
							'missing', 'nonACGTNs',
							'pcrPrimerChanges', 'aaSubstitutions',
							'totalAminoacidSubstitutions',
							'aaDeletions', 'totalAminoacidDeletions',
							'alignmentStart', 'alignmentEnd', 'alignmentScore', 
							'version','pangolin_version', 
							'scorpio_version', 'constellation_version', 'nextclade_version']]
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
	parser.add_argument("-n", "--nextclade_ver", type=check_file, required=False, default=None,
						help="Input file containing version information for Nextclade tools")
	parser.add_argument('--mode', default='accurate', required=False, help="Pangolin analysis mode. Either 'accurate' for Usher or 'fast' for pangolearn")
	parser.add_argument("--skip", action="store_true", help="Skip updates to pangolin and nextclade")
	args = parser.parse_args()

	if args.skip is False:
		if args.pangolin_ver is None: 
			update_latest_pangolin()
		else:
			update_pangolin(args.pangolin_ver)
		nextclade_dataset, nextclade_version = update_nextclade_dataset(args.nextclade_ver, False)
	else:
		nextclade_dataset, nextclade_version = update_nextclade_dataset(args.nextclade_ver, True)

	print("\nRunning Pangolin...")
	pangolin = run_pangolin(args.input_genomes, args.threads, args.mode)
	print("Running Nextclade...")
	nextclade = run_nextclade(args.input_genomes, nextclade_dataset, args.threads, nextclade_version)

	print("Collating data...")
	collate_output(nextclade, pangolin, args.output)

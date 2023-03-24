#!/usr/bin/env python3

# v1.6.0+
# signalexe.py assumes Snakefile is in current working directory (i.e., SIGNAL root)

import signal
import argparse
import subprocess, os, sys
import re
from pathlib import Path
import platform

# for compatibility between platforms
if platform.system() != 'Linux':
	signal.SIGHUP = 1

def create_parser():
	allowed = {'install': False, 'all': False, 'postprocess': False, 'ncov_tools': False}

	parser = argparse.ArgumentParser(prog='signalexe.py', description="SARS-CoV-2 Illumina GeNome Assembly Line (SIGNAL) aims to take Illumina short-read sequences and perform consensus assembly + variant calling for ongoing surveillance and research efforts towards the emergent coronavirus: Severe Acute Respiratory Syndrome Coronavirus 2 (SARS-CoV-2).")
	parser.add_argument('all', nargs='*',
						help="Run SIGNAL with all associated assembly rules. Does not include postprocessing '--configfile' or '--directory' required. The latter will automatically generate a configuration file and sample table. If both provided, then '--configfile' will take priority")
	parser.add_argument('postprocess', nargs='*',
						help="Run SIGNAL postprocessing on completed SIGNAL run. '--configfile' is required but will be generated if '--directory' is provided")
	parser.add_argument('ncov_tools', nargs='*',
						help="Generate configuration file and filesystem setup required and then execute ncov-tools quality control assessment. Requires 'ncov-tools' submodule! '--configfile' is required but will be generated if '--directory' is provided")
	parser.add_argument('install', nargs='*',
						help="Install individual rule environments and ensure SIGNAL is functional. The only parameter operable will be '--data'. Will override other operations!")
	parser.add_argument('-c', '--configfile', type=check_file, default=None,
						help="Configuration file (i.e., config.yaml) for SIGNAL analysis")
	parser.add_argument('-d', '--directory', type=check_directory, default=None,
						help="Path to directory containing reads. Will be used to generate sample table and configuration file")
	parser.add_argument('--cores', type=int, default=1, help="Number of cores. Default = 1")
	parser.add_argument('--config-only', action='store_true', help="Generate sample table and configuration file (i.e., config.yaml) and exit. '--directory' required")
	parser.add_argument('--remove-freebayes', action='store_false', help="Configuration file generator parameter. Set flag to DISABLE freebayes variant calling (improves overall speed)")
	parser.add_argument('--add-breseq', action='store_true', help="Configuration file generator parameter. Set flag to ENABLE optional breseq step (will take more time for analysis to complete)")
	parser.add_argument('-neg', '--neg-prefix', default=None, help="Configuration file generator parameter. Comma-separated list of negative control sample name(s) or prefix(es). For example, 'Blank' will cover Blank1, Blank2, etc. Recommended if running ncov-tools. Will be left empty, if not provided")
	parser.add_argument('--dependencies', action='store_true', help="Download data dependencies (under a created 'data' directory) required for SIGNAL analysis and exit. Note: Will override other parameters! (~10 GB storage required)")
	parser.add_argument('--data', default='data', help="SIGNAL install and data dependencies parameter. Set location for data dependancies. If '--dependancies' is run, a folder will be created in the specified directory. If '--config-only' or '--directory' is used, the value will be applied to the configuration file. (Upcoming feature): When used with 'SIGNAL install', any tests run will use the dependencies located at this directory. Default = 'data'")
	#parser.add_argument('--enable-test', action='store_true', help='SIGNAL install parameter. Add SIGNAL testing after environment installation using curated test data')
	parser.add_argument('-ri', '--rerun-incomplete', action='store_true', help="Snakemake parameter. Re-run any incomplete samples from a previously failed run")
	parser.add_argument('-ii', '--ignore-incomplete', action='store_true', help='Snakemake parameter. Do not check for incomplete output files')
	parser.add_argument('--unlock', action='store_true', help="Snakemake parameter. Remove a lock on the working directory after a failed run")
	parser.add_argument('-F', '--forceall', action='store_true', help='Snakemake parameter. Force the re-run of all rules regardless of prior output')
	parser.add_argument('-n', '--dry-run', action='store_true', help='Snakemake parameter. Do not execute anything and only display what would be done')
	parser.add_argument('--quiet', action='store_true', help="Snakemake parameter. Do not output any progress or rule information. If used with '--dry-run`, it will just display a summary of the DAG of jobs")
	parser.add_argument('--verbose', action='store_true', help="Snakemake parameter. Display snakemake debugging output")
	parser.add_argument('-v', '--version', action='store_true', help="Display version number")
	args, unknown = parser.parse_known_args()

	provided = []
	for opt in allowed: # ['install', 'all', 'postprocess', 'ncov_tools']
		if len(getattr(args, opt)) > 0:
			provided = provided + getattr(args, opt)
			getattr(args, opt).clear()
		else:
			continue
	
	for val in provided:
		if val.lower() in allowed:
			allowed[val.lower()] = True
		else:
			print(f"Ignoring unknown command: {val}")

	return args, allowed

def check_frontend():
	try:
		subprocess.check_call(['mamba', 'list'], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
		return 'mamba'
	except subprocess.CalledProcessError:
		return 'conda'

def check_directory(path: str) -> Path:
	"""
	Check an input directory exists and is readable
	"""
	path = Path(path)
	if path.exists() and path.is_dir():
		return path
	elif path is None:
		return None
	else:
		raise argparse.ArgumentTypeError(f"{path} can't be read")

def check_file(path: str) -> Path:
	"""
	Check an input file exists and is readable
	"""
	path = Path(path)
	if path.exists() and path.is_file():
		return path
	elif path is None:
		return None
	else:
		raise argparse.ArgumentTypeError(f"{path} can't be read")

def check_single_replicate_and_resolve_paths(project_directory):
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

def check_submodule(exec_dir):
	print("Checking for ncov-tools!")
	if (not os.path.exists(os.path.join(exec_dir, 'ncov-tools'))) or (not os.path.exists(os.path.join(exec_dir, 'ncov-tools', 'workflow', 'Snakefile'))):
		try:
			print("Updating ncov-tools!")
			subprocess.run(['git', 'submodule', 'update', '--init', '--recursive'])
			return True
		except subprocess.CalledProcessError:
			print("Could not find nor update the required 'ncov-tools' directory! Manually download/update and try again!")
			sys.exit(1)
	else:
		return False
		
def write_sample_table(sample_data, output_table):
	"""
	Generate sample table from parsed folder read data
	"""
	header = ",".join(['sample', 'r1_path', 'r2_path']) + '\n'
	with open(output_table, 'w') as out_fh:
		out_fh.write(header)
		for sample in sample_data:
			out_fh.write(",".join(sample) + '\n')

def download_dependencies(dir_name):
	script = os.path.join(script_path, 'scripts', 'get_data_dependencies.sh')
	subprocess.run(['bash', script, '-d', dir_name, '-a', 'MN908947.3'])

def generate_sample_table(project_directory, project_name):
	"""
	Bash shell script generalizes the search for R1/R2 files, allowing multiple replicates that SIGNAL supports."
	"""
	script=os.path.join(script_path, 'scripts', "generate_sample_table.sh")
	out_table = project_name + "_sample_table.csv"
	subprocess.run(['bash', script, '-d', project_directory, '-n', out_table])

def write_config_file(run_name, config_file, data_directory, opt_tasks):
### opt_tasks = [args.breseq, args.freebayes, [args.neg_prefix]] - latter only applies to SIGNAL v1.5.8 and earlier

	config = f"""# This file contains a high-level summary of pipeline configuration and inputs.
# It is ingested by the Snakefile, and also intended to be human-readable.

# Sample table (can be created using example_sample_table.csv)
samples: "{run_name}_sample_table.csv"

# Folder to place all results (can be relative or absolute path)
result_dir: "{run_name}_results_dir"

# Minimum quality threshold, used in 'trim_galore' and 'ivar trim'
min_qual: 20

# Minimum read length to retain after trimming, used in 'trim_galore' and 'ivar trim'
min_len: 20

# Path from snakemake dir to .bed file defining amplicon primer scheme
scheme_bed: 'resources/primer_schemes/artic_v3/nCoV-2019.bed'

# Path from snakemake dir to bwa indexed human + viral reference genome
composite_reference: "{data_directory}/composite_human_viral_reference.fna"

# Used as bwa reference genome when removing host sequences.
# Also used as 'ivar' reference genome in variant detection + consensus.
# Used as -r,-g arguments to 'quast'
# contig needed for hostremoval filtering script
viral_reference_contig_name: 'MN908947.3'
viral_reference_genome: "{data_directory}/MN908947.3.fasta"
viral_reference_feature_coords: "{data_directory}/MN908947.3.gff3"

# breseq_reference must be defined if run_breseq == True
run_breseq: {opt_tasks[0]}
# Used as --reference argument to 'breseq'
breseq_reference: "{data_directory}/MN908947.3.gbk"

# run freebayes for variant and consensus calling (as well as ivar)
run_freebayes: {opt_tasks[1]}

# Used as --db argument to 'kraken2'
kraken2_db: "{data_directory}/Kraken2/db"

# For Ivar's amplicon filter 
# https://github.com/andersen-lab/ivar/commit/7027563fd75581c78dabc6040ebffdee2b24abe6
# must be set to nothing if you are not wanting to use this setting
# and "-f primer_pairs.tsv" with the correct file path if you do wish to use it
primer_pairs_tsv: 

# Consensus and variant calling ivar/samtools params from https://github.com/connor-lab/ncov2019-artic-nf/blob/master/conf/illumina.config
mpileup_depth: 100000
# ivar/freebayes frequency threshold to build consensus
var_freq_threshold: 0.75
# Minimum coverage depth to call variant
var_min_coverage_depth: 10
# iVar/freebayes frequency threshold to call variant (ivar variants: -t )
var_min_freq_threshold: 0.25
# iVar/freebayes minimum mapQ to call variant (ivar variants: -q)
var_min_variant_quality: 20

# Toggle faster Pangolin analysis at the cost of accuracy (uses Pangolearn instead of Usher)
# Use for significantly larger datasets
pangolin_fast: False

# Versions of software related to lineage calling (use numbers only, i.e., 3.1.1). Dates (YYYY-mm-dd) are accepted for pangolearn. Leave blank for latest version(s).
pangolin: 
constellations:
scorpio:

# Required for Pangolin <v4.0
pangolearn: 
pango-designation:

# Required for Pangolin v4+
pangolin-data:

# Versions for Nextclade (software & datasets)
# Software version (nextclade) should use numbers only (i.e., 1.11.0)
# Be as specific as possible with the desired dataset tag (nextclade-data). Can accept dates (YYYY-mm-dd) alone, but will assume corresponding timestamp (HH:MM:SS) 
# Typical tag format is YYYY-mm-ddTHH:MM:SSZ
# Leave blank for latest versions
# Setting nextclade-include-recomb to False will download the recombinant-sequence free version of the Nextclade database
nextclade:
nextclade-data:
nextclade-include-recomb: True

# ANYTHING BELOW IS ONLY NEEDED IF USING NCOV-TOOLS SUMMARIES
# Path from snakemake dir to .bed file defining the actual amplicon locations not the primers
amplicon_loc_bed: 'resources/primer_schemes/artic_v3/ncov-qc_V3.scheme.bed'

# fasta of sequences to include with pangolin phylogeny
phylo_include_seqs: "{data_directory}/blank.fasta"

# List of negative control sample names or prefixes (i.e., ['Blank'] will cover Blank1, Blank2, etc.)
negative_control_prefix: {opt_tasks[2]}"""

	with open(config_file, 'w') as fh:
		fh.write(config)

def install_signal(frontend, data='data'):
	"""
	Install SIGNAL dependencies per rule and test using a sample dataset, if desired
	"""
	dep_snakefile = os.path.join(script_path, 'resources', 'dependencies')
	assert os.path.exists(dep_snakefile)
	try:
		subprocess.run(f"snakemake -s {dep_snakefile} --conda-frontend {frontend} --cores 1 --use-conda --conda-prefix=$PWD/.snakemake/conda", shell=True, check=True)
	except subprocess.CalledProcessError: 
		print("Installation of environments failed!")
		sys.exit(1)
	
	### TODO: Test SIGNAL with curated data
	if os.path.exists(data):
		pass
	

if __name__ == '__main__':
	# note: add root_dir to determine the root directory of SIGNAL
	script_path = os.path.join(os.path.abspath(sys.argv[0]).rsplit("/",1)[0])
	args, allowed = create_parser()
	version = 'v1.6.0'
	alt_options = []
	
	if args.version:
		print(f"{version}")
		sys.exit(0)
	
	conda_frontend = check_frontend() # 'mamba' or 'conda'
	
	if allowed['install']:
		install_signal(conda_frontend, args.data)
		print("Installation of environments completed successfully!")
		sys.exit(0)
	
	if args.dependencies:
		print("Downloading necessary reference and dependency files!")
		download_dependencies(args.data)
		print("Complete!")
		sys.exit(0)
		
	if args.configfile is None:
		assert args.directory is not None, "Please provide '--directory' to proceed! ('--configfile' if a configuration file already exists!)"
		run_name = args.directory.name
		generate_sample_table(args.directory, run_name)
		config_file = run_name + "_config.yaml"
		if args.neg_prefix is not None:
			neg = [pre.replace(" ","") for pre in args.neg_prefix.split(",")]
		else:
			neg = [args.neg_prefix]
		write_config_file(run_name, config_file, args.data, [args.add_breseq, args.remove_freebayes, neg])
		if args.config_only:
			print("Configuration file and sample table generated!")
			sys.exit(0)
	else:
		config_file = args.configfile
	
	if not any([allowed[x] for x in allowed]):
		print("No task specified! Please provide at least one of 'all', 'postprocess', or 'ncov_tools'! See 'signalexe.py -h' for details!")
		sys.exit(1)
	else:
		if args.verbose: alt_options.append('--verbose')
		if args.quiet: alt_options.append('--quiet')
		if args.unlock: alt_options.append('--unlock')
		if args.forceall: alt_options.append('--forceall')
		if args.dry_run: alt_options.append('--dry-run')
		if args.rerun_incomplete: alt_options.append('--rerun-incomplete')
		if args.ignore_incomplete: alt_options.append('--ignore-incomplete')
		opt = " ".join(alt_options)
		for task in allowed:
			if (allowed[task] is True) and (task != 'install'):
				print(f"Running SIGNAL {task}!")
				try:
					subprocess.run(f"snakemake --configfile {config_file} --cores={args.cores} {task} -kp {opt}", shell=True, check=True)
				except subprocess.CalledProcessError:
					if task == "ncov_tools":
						mod = check_submodule(os.getcwd())
						if mod:
							if opt.split(" ")[-1] == '--rerun-incomplete': # remove redundant flag
								opt = " ".join(opt.split(" ")[:-1])
							try:
								print("Retrying...ncov-tools!")
								subprocess.run(f"snakemake --configfile {config_file} --cores={args.cores} {task} -kp --rerun-incomplete {opt}", shell=True, check=True)
							except subprocess.CalledProcessError:
								print(f"Some jobs failed while running SIGNAL {task}! Check snakemake logs and the ncov-tools directory for additional details!")
								continue
						else:
							print(f"Some jobs failed while running SIGNAL {task}! Check snakemake logs and the ncov-tools directory for additional details!")
							continue
					elif task == 'all':
						print(f"Some jobs failed while running SIGNAL {task}! This does NOT necessarily mean your run was erroneous! Samples that failed assembly can be found in 'failed_samples.log'! If no such file exists or is blank, check your inputs and logs and try again!")
						continue
					elif task == 'postprocess':
						print(f"Some jobs failed while running SIGNAL {task}! Some output files may be missing! Check SIGNAL results and logs and try again!")
						continue
					else:
						print(f"Some jobs failed while running SIGNAL {task}! Check SIGNAL inputs, logs, and results and try again!")
	
	print("SIGNAL run complete! Check corresponding snakemake logs for any details!")
	sys.exit(0)

import os
import shutil
import subprocess
import fileinput
import glob

def set_up():
	print("Writing config for ncov to ncov-tools/config.yaml")

	exec_dir = snakemake.params['exec_dir']
	result_dir = os.path.basename(snakemake.params['result_dir']) # basename of SIGNAL result directory

	# result_root = os.path.abspath(os.path.join(exec_dir, result_dir, "ncov-tools-results"))
	# if os.path.exists(result_root):
		# shutil.rmtree(result_root)
	# ncov_output = ["plots", "lineages", "qc_analysis"]
	# for dir in ncov_output:
		# os.makedirs(os.path.join(result_root, dir))

### Create data directory within ncov-tools
	data_root = os.path.abspath(os.path.join(exec_dir, 'ncov-tools', "%s" %(result_dir)))
	if os.path.exists(data_root):
		shutil.rmtree(data_root)
	os.mkdir(data_root)

	# snakemake_dir = os.path.join(exec_dir, 'ncov-tools', '.snakemake')
	# if os.path.exists(snakemake_dir):
		# shutil.rmtree(snakemake_dir)

### Pull negative samples (based on common identifiers)
	#neg_names = ("Negative", "NEG", "PCR-NEG", "UTM", "BLANK", "Blank", "blank")
	neg_names = tuple(snakemake.params['negative_control_prefix'])
	neg_samples = set()
	with open(snakemake.params['sample_csv_filename']) as fh:
		fh.readline()
		for line in fh:
			id = line.split(",")[0]
			if id.startswith(neg_names):
				neg_samples.add(id)
	neg_list = list(neg_samples)
	print("Negative control samples found include: %s" %(neg_list))


### config.yaml parameters
	config = {'data_root': f"'{data_root}'",
			  'run_name': f"'{result_dir}'", # name ncov-tools output files with name of SIGNAL results directory (default: "default")
			  'amplicon_bed': f"'{snakemake.params['amplicon_bed']}'", #grab from signal snakemake config
			  'reference_genome': f"'{snakemake.params['viral_reference_genome']}'", #grab from signal snakemake config
			  'platform': 'illumina',
			  'primer_bed': f"'{snakemake.params['primer_bed']}'",
			  'bed_type': "unique_amplicons",
			  'offset': 0,
			  'completeness_threshold': 0.9,
			  'bam_pattern': "'{data_root}/{sample}.bam'", # symlink files following this
			  'primer_trimmed_bam_pattern': "'{data_root}/{sample}.mapped.primertrimmed.sorted.bam'",
			  'consensus_pattern': "'{data_root}/{sample}.consensus.fasta'", # symlink files following this
			  'variants_pattern': "'{data_root}/{sample}.variants.tsv'",
			  'assign_lineages': 'true',
			  'tree_include_consensus': f"'{snakemake.params['phylo_include_seqs']}'",
			  'negative_control_samples': f"{neg_list}",
			  'mutation_set': 'spike_mutations',
			  'output_directory': f"{result_dir}_ncovresults"}

	with open(os.path.join(exec_dir, 'ncov-tools', 'config.yaml'), 'w') as fh:
			for key, value in config.items():
					fh.write(f"{key}: {value}\n")

	print("Linking files to ncov")
	for bam in snakemake.input['bams']:
		sample = bam.split('/')[0]
		ln_path = f"{data_root}/{sample}.bam"
		if not os.path.exists(ln_path):
				os.link(bam, ln_path)


	for primer_trimmed_bam in snakemake.input['primertrimmed_bams']:
		sample = primer_trimmed_bam.split('/')[0]
		ln_path = f"{data_root}/{sample}.mapped.primertrimmed.sorted.bam"
		if not os.path.exists(ln_path):
			os.link(primer_trimmed_bam, ln_path)

	for variants in snakemake.input['variants']:
		sample = variants.split('/')[0]
		ln_path = f"{data_root}/{sample}.variants.tsv"
		if not os.path.exists(ln_path):
			os.link(variants, ln_path)

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

	# os.chdir(os.path.join(exec_dir, 'ncov-tools'))
	#return exec_dir, result_root, result_dir
	return exec_dir, result_dir

def run_all():
	os.system(f"snakemake -s workflow/Snakefile --cores {snakemake.threads} all")

def move(cwd, dest, prefix):
	if os.path.exists(os.path.join(cwd, "ncov-tools", "plots")):
		extensions = ["_amplicon_coverage_heatmap.pdf", "_amplicon_covered_fraction.pdf", "_depth_by_position.pdf", "_tree_snps.pdf"]
		for ext in extensions:
			for file in glob.glob(os.path.join(cwd, "ncov-tools", "plots", "default"+ext)):
				try:
					shutil.copy(file, os.path.join(dest, "plots", prefix+ext))
				except IOError:
					print("Missing file %s" %(prefix+ext))
	else:
		print("Missing ncov-tools 'plots' directory")

	if os.path.exists(os.path.join(cwd, "ncov-tools", "lineages")):
		for file in glob.glob(os.path.join(cwd, "ncov-tools", "lineages", "default_lineage_report.csv")):
			try:
				shutil.copy(file, os.path.join(dest, "lineages", prefix+"_lineage_report.csv"))
			except IOError:
				print("Missing file %s_lineage_report.csv" %(prefix))
	else:
		print("Missing ncov-tools 'lineages' directory")

	if os.path.exists(os.path.join(cwd, "ncov-tools", "qc_analysis")):
		extensions = ["_aligned-delim.fasta.log", "_aligned-delim.iqtree.log", "_aligned.fasta", "_aligned.fasta.insertions.csv", \
					"_aligned.fasta.log", "_alleles.tsv", "_tree.nwk", "_tree_raw.nwk", "_consensus.fasta"]
		for ext in extensions:
			for file in glob.glob(os.path.join(cwd, "ncov-tools", "qc_analysis", "default"+ext)):
				try:
					shutil.copy(file, os.path.join(dest, "qc_analysis", prefix+ext))
				except IOError:
					print("Missing file %s" %(prefix+ext))
	else:
		print("Missing ncov-tools 'qc_analysis' directory")

if __name__ == '__main__':
	exec_dir, result_dir = set_up()
	print("Don't forget to update the config.yaml file as needed prior to running ncov-tools.")
	#exec_dir, result_root, result_dir = set_up()
	#run_all()
	#move(exec_dir, result_root, result_dir)

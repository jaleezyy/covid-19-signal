import os
import shutil
import subprocess
import fileinput
import glob

def set_up():
	print("Writing config for ncov to ncov-tools/config.yaml")

	exec_dir = snakemake.params['exec_dir']
	result_dir = snakemake.params['result_dir'] # basename of SIGNAL result directory
	
	result_root = os.path.abspath(os.path.join(exec_dir, result_dir, "ncov-tools-results"))
	if os.path.exists(result_root):
		shutil.rmtree(result_root)
	ncov_output = ["plots", "lineages", "qc_analysis"]
	for dir in ncov_output:
		os.makedirs(os.path.join(result_root, dir))

	data_root = os.path.abspath(os.path.join(exec_dir, 'ncov-tools', 'data'))
	if os.path.exists(data_root):
		shutil.rmtree(data_root)
	os.mkdir(data_root)

	snakemake_dir = os.path.join(exec_dir, 'ncov-tools', '.snakemake')
	if os.path.exists(snakemake_dir):
		shutil.rmtree(snakemake_dir)

	config = {'data_root': f"'{data_root}'",
			  #'run_name': f"'{result_dir}'", # name ncov-tools output files with name of SIGNAL results directory (default: "default")
			  'amplicon_bed': f"'{snakemake.params['amplicon_bed']}'", #grab from signal snakemake config
			  'reference_genome': f"'{snakemake.params['viral_reference_genome']}'", #grab from signal snakemake config
			  'bam_pattern': "'{data_root}/{sample}.bam'", # symlink files following this
			  'consensus_pattern': "'{data_root}/{sample}.consensus.fasta'", # symlink files following this
			  'tree_include_consensus': f"'{snakemake.params['phylo_include_seqs']}'",
			  #'metadata': snakemake.config['samples'], # causes error get example
			  'assign_lineages': 'true'}

	with open(os.path.join(exec_dir, 'ncov-tools', 'config.yaml'), 'w') as fh:
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

	os.chdir(os.path.join(exec_dir, 'ncov-tools'))
	return exec_dir, result_root, result_dir

def run_all_qc_sequencing():
	os.system(f"snakemake -s qc/Snakefile --cores {snakemake.threads} all_qc_sequencing")

def run_all_qc_analysis():
	os.system(f"snakemake -s qc/Snakefile --cores {snakemake.threads} all_qc_analysis")

def run_all_qc_summary():
	os.system(f"snakemake -s qc/Snakefile --cores {snakemake.threads} all_qc_summary")

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
	exec_dir, result_root, result_dir = set_up()
	run_all_qc_sequencing()
	run_all_qc_analysis()
	# seems to include hard-coded issues
	#run_all_qc_summary()
	move(exec_dir, result_root, result_dir)

#!/usr/bin/env python

import os, sys
import shutil
import subprocess
import fileinput
import glob

def link_ivar(root, neg, failed, replace=False):
	print("Linking iVar files to ncov-tools!")

	for variants in snakemake.input['variants']:
		sample = variants.split('/')[0]
		if (sample in failed) and (sample not in neg):
			continue
		ln_path = f"{root}/{sample}.variants.tsv"
		try:
			if (not os.path.exists(ln_path)) or (replace is True):
				# remove equivalent link for Freebayes VCF (namely for replace=True)
				if os.path.exists(f"{root}/{sample}.variants.vcf"):
					os.remove(f"{root}/{sample}.variants.vcf")
				os.link(variants, ln_path)
		except FileExistsError: # equal link ~ error
			if replace:
				os.remove(ln_path)
				os.link(variants, ln_path)
			else:
				pass

	for consensus in snakemake.input['consensus']:
		sample = consensus.split('/')[0]
		if (sample in failed) and (sample not in neg):
			continue
		ln_path = f"{root}/{sample}.consensus.fasta"
		try:
			if (not os.path.exists(ln_path)) or (replace is True):
				os.link(consensus, ln_path)
		except FileExistsError: # more likely given same link path
			if replace:
				os.remove(ln_path)
				os.link(variants, ln_path)
			else:
				pass
		for line in fileinput.input(ln_path, inplace=True):
			if line.startswith(">"):
				new_header = str(">"+sample)
				new_line = line.replace(line, new_header)
				print(new_line, end='\n')
			else:
				print(line, end='\n')

# take sample name from iVar results, redirect to where corresponding FreeBayes should be
# if FreeBayes file cannot be found, break from loop, replace all with iVar
def link_freebayes(root, neg, failed):
	print("Linking FreeBayes files to ncov-tools!")

	for variants in snakemake.input['variants']:
		sample = variants.split('/')[0]
		if (sample in failed) and (sample not in neg):
			continue
		expected_path = os.path.join(sample, 'freebayes', sample+'.variants.norm.vcf')
		if not os.path.exists(expected_path):
			print("Missing FreeBayes variant file! Switching to iVar input!")
			link_ivar(root, neg, failed, replace=True)
			vcf_method = False
			break
		else:
			vcf_method = True
			ln_path = f"{root}/{sample}.variants.vcf"
			# redundant check however, it may play a roll in re-runs
			if os.path.exists(f"{root}/{sample}.variants.tsv"):
				os.remove(f"{root}/{sample}.variants.tsv")
			if not os.path.exists(ln_path):
				os.link(expected_path, ln_path)

	if vcf_method:
		for consensus in snakemake.input['consensus']:
			sample = consensus.split('/')[0]
			if (sample in failed) and (sample not in neg):
				continue
			expected_path = os.path.join(sample, 'freebayes', sample+'.consensus.fasta')
			if not os.path.exists(expected_path):
				print("Missing FreeBayes variant file! Switching to iVar input!")
				link_ivar(root, neg, failed, replace=True)
				break
			else:
				ln_path = f"{root}/{sample}.consensus.fasta"
				if not os.path.exists(ln_path):
					os.link(expected_path, ln_path)
				for line in fileinput.input(ln_path, inplace=True):
					if line.startswith(">"):
						new_header = str(">"+sample)
						new_line = line.replace(line, new_header)
						print(new_line, end='\n')
					else:
						print(line, end='\n')
					
	return vcf_method

def set_up():
	print("Writing config.yaml for ncov-tools to ncov-tools/config.yaml")
	#vcf_variant = True # set default

	exec_dir = snakemake.params['exec_dir']
	result_dir = os.path.basename(snakemake.params['result_dir']) # basename of SIGNAL result directory

	### Pull pangolin version number (i.e., version '3' or '4')
	pangolin = str(snakemake.params['pangolin']).split(".")[0].lower().strip("v")
	try:
		assert pangolin == "3" or pangolin == "4" # directly supported versions
	except AssertionError:
		# import urllib.request as web
		# commit_url = web.urlopen(f"https://github.com/cov-lineages/pangolin/releases/latest").geturl()
		# pangolin = commit_url.split("/")[-1].split(".")[0].lower().strip("v") 
		# latest version (should ensure temporary compatibility)
		installed_versions = subprocess.run(["pangolin", "--all-versions"],
								check=True,
								stdout=subprocess.PIPE)
		installed_versions = installed_versions.stdout.decode('utf-8')
		installed_ver_dict = {}
		for dep_ver in map(str.strip, installed_versions.split('\n')):
		# skip empty line at end
			if len(dep_ver) == 0:
				continue
			try:
				dependency, version = dep_ver.split(': ')
			except ValueError:
				continue
			if dependency == 'pangolin':
				pangolin = str(version).split(".",1)[0].strip('v')
				break

### Create data directory within ncov-tools
	data_root = os.path.abspath(os.path.join(exec_dir, 'ncov-tools', "%s" %(result_dir)))
	if os.path.exists(data_root):
		shutil.rmtree(data_root)
	os.mkdir(data_root)

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

### Pull failed samples (SIGNAL log file: failed_samples.log)
	if os.path.exists(snakemake.params['failed']):
		with open(snakemake.params['failed']) as fail:
			failed_list = [i.strip() for i in fail.readlines()[1:]]
	else:
		failed_list = []
	print("Failed samples found include: %s" %(failed_list))

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
			  'output_directory': f"{result_dir}_ncovresults",
			  'pangolin_version': f"'{pangolin}'",
			  'pango_analysis_mode': f"'{snakemake.params['mode']}'"
			  }

	print("Linking alignment BAMs to ncov-tools!")
	for bam in snakemake.input['bams']:
		sample = bam.split('/')[0]
		# if sample failed and not a negative, skip linking
		if (sample in failed_list) and (sample not in neg_list):
			continue
		ln_path = f"{data_root}/{sample}.bam"
		if not os.path.exists(ln_path):
			os.link(bam, ln_path)

	for primer_trimmed_bam in snakemake.input['primertrimmed_bams']:
		sample = primer_trimmed_bam.split('/')[0]
		if (sample in failed_list) and (sample not in neg_list):
			continue
		ln_path = f"{data_root}/{sample}.mapped.primertrimmed.sorted.bam"
		if not os.path.exists(ln_path):
			os.link(primer_trimmed_bam, ln_path)
			
	if snakemake.params['freebayes_run']:
		vcf_variant = link_freebayes(data_root, neg_list, failed_list)
		if vcf_variant:
			config['variants_pattern'] = "'{data_root}/{sample}.variants.vcf'"
		else:
			pass # keep as TSV
	else:
		link_ivar(data_root, neg_list, failed_list, replace=False)

	with open(os.path.join(exec_dir, 'ncov-tools', 'config.yaml'), 'w') as fh:
		for key, value in config.items():
			fh.write(f"{key}: {value}\n")
			
	return exec_dir, result_dir, data_root

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
	exec_dir, result_dir, data_root = set_up()
	run_script = os.path.join(exec_dir, 'scripts', 'run_ncov_tools.sh')
	#print("Don't forget to update the config.yaml file as needed prior to running ncov-tools.")
	print("Running ncov-tools using %s cores!" %(snakemake.threads))

	subprocess.run([run_script, '-c', str(snakemake.threads), '-s', str(result_dir)])
	
	# clean up
	shutil.rmtree(data_root)
	
	#run_all()
	#move(exec_dir, result_root, result_dir)

# To run the pipeline, do:
#
#    snakemake -kp --cores=NCORES --use-conda --conda-prefix=$HOME/.snakemake all
#    snakemake -p --cores=1 postprocess
#
# Note that the pipeline postprocessing ('snakemake postprocess') is separated from
# the rest of the pipeline ('snakemake all').  This is because in a multi-sample run,
# it's likely that at least one pipeline stage will fail.  The postprocessing script
# should handle failed pipeline stages gracefully, by substituting placeholder values
# when expected pipeline output files are absent.  However, this confuses snakemake's
# dependency tracking, so there seems to be no good alternative to separating piepline
# processing and postprocessing into 'all' and 'postprocess' targets.
#
# Related: because pipeline stages can fail, we recommend running 'snakemake all'
# with the -k flag ("Go on with independent jobs if a job fails").


####################################################################################################

from snakemake.utils import validate
import pandas as pd
import os, sys

# The config file contains a high-level summary of pipeline configuration and inputs.
# It is ingested by the Snakefile, and also intended to be human-readable.
# For an example config file, see pipeline/example_config.yaml in the covid-19-sequencing repo.

# read and validate config.yaml
if '--configfile' in sys.argv:
    config_filename = os.path.abspath(sys.argv[sys.argv.index('--configfile')+1])
    # arguments don't line up
    if not os.path.exists(config_filename):
        print("Invalid filepath for configfile. Looking for default config.yaml")
        configfile: "config.yaml"
        config_filename = os.path.join(os.getcwd(), "config.yaml")
    else:
        configfile: config_filename
else:
    configfile: "config.yaml"
    config_filename = os.path.join(os.getcwd(), "config.yaml")

validate(config, 'resources/config.schema.yaml')

# read and validate sample table specified in config.schema.yaml
samples = pd.read_table(config['samples'], sep=',')
validate(samples, 'resources/fasta.sample.schema.yaml')

# set output directory
exec_dir = os.getcwd()
workdir: os.path.abspath(config['result_dir'])

# get sample names 
sample_names = sorted(samples['sample'].drop_duplicates().values)

# get lineage calling versions
versions = {'pangolin': config['pangolin'],
            'pangolearn': config['pangolearn'],
            'constellations': config['constellations'],
            'scorpio': config['scorpio'],
            'pango-designation': config['pango-designation'],
            'pangolin-data': config['pangolin-data'],
            'nextclade': config['nextclade'],
            'nextclade-data': config['nextclade-data'],
            'nextclade-recomb': config['nextclade-include-recomb']
            }

# Determine Pangolin analysis mode
if config['pangolin_fast']:
    pango_speed = 'fast'
else:
    pango_speed = 'accurate'

def get_input_fastq_files(sample_name, r):
    sample_fastqs = samples[samples['sample'] == sample_name]
    if r == '1':
        relpath = sample_fastqs['r1_path'].values[0]
    elif r == '2':
        relpath = sample_fastqs['r2_path'].values[0]

    return os.path.abspath(os.path.join(exec_dir, relpath))

def get_pooled_fastq_files(sample_name, r):
    sample_fastqs = samples[samples['sample'] == sample_name]
    if r == '1':
        relpath = sample_fastqs['r1_path'].values
    elif r == '2':
        relpath = sample_fastqs['r2_path'].values

    return [ os.path.abspath(os.path.join(exec_dir, r)) for r in relpath ]

def get_fasta_files(sample_name):
    sample_fasta = samples[samples['sample'] == sample_name]
    relpath = sample_fasta['fasta_path'].values[0]
    

    return os.path.abspath(os.path.join(exec_dir, relpath))

######################################   High-level targets   ######################################

rule consensus:
    input: expand('{sn}/input_data/{sn}.consensus.fa', sn=sample_names)

rule quast:
    input: 
        expand('{sn}/quast/{sn}_quast_report.html', sn=sample_names),
        expand('{sn}/quast/{sn}_quast_report.tsv', sn=sample_names)

rule lineages:
    input: 
        'input_pangolin_versions.txt',
        'input_nextclade_versions.txt',
        'lineage_assignments.tsv'

rule config_sample_log:
    input: 
        config_filename,
        config['samples']

rule all:
    input:
        rules.consensus.input,
        rules.config_sample_log.input,
        rules.lineages.input,
        rules.quast.input

rule postprocess:
    conda: 
        'conda_envs/postprocessing.yaml'
    params:
        #sample_csv_filename = os.path.join(exec_dir, config['samples']),
        input_lineage = 'lineage_assignments.tsv',
        postprocess_script_path = os.path.join(exec_dir, 'scripts', 'signal_partial_postprocess.py'),
        stats = expand('{sn}/input_data/{sn}_stats.txt', sn=sample_names),
        quast = expand('{sn}/quast/{sn}_quast_report.tsv', sn=sample_names)
    input:
        #quast = expand('{sn}/quast/{sn}_quast_report.tsv', sn=sample_names)
        #stats = expand('{sn}', sn=sample_names)
        #stats = expand('{sn}/input_data/{sn}_stats.txt', sn=sample_names)
        
    shell:
        '{params.postprocess_script_path} -i {params.input_lineage} -q {params.stats} {params.quast}'


rule ncov_tools:
    # can't use the one in the ncov-tool dir as it has to include snakemake
    conda:
        'ncov-tools/workflow/envs/environment.yml'
    params:
        exec_dir = exec_dir,
        sample_csv_filename = os.path.join(exec_dir, config['samples']),
        result_dir = os.path.basename(config['result_dir']),
        amplicon_bed = os.path.join(exec_dir, config['amplicon_loc_bed']),
        primer_bed = os.path.join(exec_dir, config['scheme_bed']),
        viral_reference_genome = os.path.join(exec_dir, config['viral_reference_genome']),
        phylo_include_seqs = os.path.join(exec_dir, config['phylo_include_seqs']),
        negative_control_prefix = config['negative_control_prefix']
    input:
        consensus = expand('{sn}/core/{sn}.consensus.fa', sn=sample_names),
        primertrimmed_bams = expand("{sn}/core/{sn}_viral_reference.mapping.primertrimmed.sorted.bam", sn=sample_names),
        bams = expand("{sn}/core/{sn}_viral_reference.mapping.bam", sn=sample_names),
        variants = expand("{sn}/core/{sn}_ivar_variants.tsv", sn=sample_names)
    script: "scripts/ncov-tools.py"
        
        
################################# Copy config and sample table to output folder ##################

rule copy_config_sample_log:
    output: 
        config = os.path.basename(config_filename),
        sample_table=config["samples"]
    input:
        origin_config = os.path.join(exec_dir, os.path.relpath(config_filename, exec_dir)),
        origin_sample_table = os.path.join(exec_dir, config['samples'])
    shell:
        """
        cp {input.origin_config} {output.config}
        cp {input.origin_sample_table} {output.sample_table}
        """

#################################   Based on scripts/assemble.sh   #################################

rule link_raw_data:
    priority: 4
    output:
        '{sn}/input_data/{sn}.consensus.fa'
    input:
        lambda wildcards: get_fasta_files(wildcards.sn)
    params:
        sample = '{sn}',
        outdir = '{sn}/input_data'
    shell:
        "sed 's/>consensus sequence_BreakPoint_1/>{params.sample}/g' {input} > {output} && "
        "cp $(dirname {input})/{params.sample}_*StatInfo.txt {params.outdir}/{params.sample}_stats.txt 2>/dev/null || echo 'WARNING: No stats file found!'"
        #TODO: correct format, dummy approach for initial testing
        #'ln -s {input} {output}'


##################################  Based on scripts/quast.sh   ####################################

rule run_lineage_assignment:
    threads: 4
    conda: 'conda_envs/assign_lineages.yaml'
    output:
        pango_ver_out = 'input_pangolin_versions.txt',
        nextclade_ver_out = 'input_nextclade_versions.txt',
        lin_out = 'lineage_assignments.tsv'
    input:
        expand('{sn}/input_data/{sn}.consensus.fa', sn=sample_names)
    params:
        pangolin_ver = versions['pangolin'],
        pangolearn_ver = versions['pangolearn'],
        constellations_ver = versions['constellations'],
        scorpio_ver = versions['scorpio'],
        designation_ver = versions['pango-designation'],
        data_ver = versions['pangolin-data'],
        nextclade_ver = versions['nextclade'],
        nextclade_data = versions['nextclade-data'],
        nextclade_recomb = versions['nextclade-recomb'],
        analysis_mode = pango_speed,
        assignment_script_path = os.path.join(exec_dir, 'scripts', 'assign_lineages.py')
    shell:
        "echo -e 'pangolin: {params.pangolin_ver}\nconstellations: {params.constellations_ver}\nscorpio: {params.scorpio_ver}\npangolearn: {params.pangolearn_ver}\npango-designation: {params.designation_ver}\npangolin-data: {params.data_ver}' > {output.pango_ver_out} && "
        "echo -e 'nextclade: {params.nextclade_ver}\nnextclade-dataset: {params.nextclade_data}\nnextclade-include-recomb: {params.nextclade_recomb}' > {output.nextclade_ver_out} && "
        'cat {input} > all_genomes.fa && '
        '{params.assignment_script_path} -i all_genomes.fa -t {threads} -o {output.lin_out} -p {output.pango_ver_out} -n {output.nextclade_ver_out} --mode {params.analysis_mode}'

rule run_quast:
    threads: 1
    conda: 'conda_envs/assembly_qc.yaml'
    output:
         '{sn}/quast/{sn}_quast_report.html',
         '{sn}/quast/{sn}_quast_report.tsv'
    input:
         '{sn}/input_data/{sn}.consensus.fa'
    log:
         '{sn}/quast/{sn}_quast.log'
    benchmark:
        "{sn}/benchmarks/{sn}_run_quast.benchmark.tsv"
    params:
         outdir = '{sn}/quast',
         genome = os.path.join(exec_dir, config['viral_reference_genome']),
         fcoords = os.path.join(exec_dir, config['viral_reference_feature_coords']),
         sample_name = '{sn}_quast_report',
         unlabelled_reports = '{sn}/quast/report.*'
    shell:
         'quast {input} -r {params.genome} -g {params.fcoords} --output-dir {params.outdir} --threads {threads} >{log} && '
         'for f in {params.unlabelled_reports}; do mv $f ${{f/report/{params.sample_name}}}; done'


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
validate(samples, 'resources/sample.schema.yaml')

# manual assignment of breseq reference
try:
    if os.path.exists(config['breseq_reference']):
        breseq_ref = config['breseq_reference']
    else:
        breseq_ref = ""
except TypeError:
    breseq_ref = ""

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

# determine raw FASTQ handling
# if duplicate sample names in table, run legacy concat_and_sort
if samples['sample'].duplicated().any():
    print("Duplicate sample names in sample table. Assuming multi-lane samples exist")
    ruleorder: concat_and_sort > link_raw_data
else:
    ruleorder: link_raw_data > concat_and_sort

# Determine Pangolin analysis mode
if config['pangolin_fast']:
    pango_speed = 'fast'
else:
    pango_speed = 'accurate'

######################################   High-level targets   ######################################
rule raw_read_data_symlinks:
    input: expand('{sn}/raw_fastq/{sn}_R{r}.fastq.gz', sn=sample_names, r=[1,2])

rule remove_adapters:
    input: expand('{sn}/adapter_trimmed/{sn}_R{r}_val_{r}.fq.gz', sn=sample_names, r=[1,2]),
           expand('{sn}/adapter_trimmed/{sn}_R{r}_val_{r}_posttrim_filter.fq.gz', sn=sample_names, r=[1,2]),

rule host_removed_raw_reads:
    input: expand('{sn}/host_removal/{sn}_R{r}.fastq.gz', sn=sample_names, r=[1,2]),

rule fastqc:
    input: expand('{sn}/raw_fastq/{sn}_R{r}_fastqc.html', sn=sample_names, r=[1,2]),
           expand('{sn}/adapter_trimmed/{sn}_R{r}_val_{r}_fastqc.html', sn=sample_names, r=[1,2]),
           expand('{sn}/mapped_clean_reads/{sn}_R{r}_fastqc.html', sn=sample_names, r=[1,2])

rule clean_reads:
    input:
       expand("{sn}/core/{sn}_viral_reference.mapping.primertrimmed.bam", sn=sample_names),
       expand('{sn}/mapped_clean_reads/{sn}_R{r}.fastq.gz', sn=sample_names, r=[1,2])

rule consensus:
    input: expand('{sn}/core/{sn}.consensus.fa', sn=sample_names)

rule ivar_variants:
    input: expand('{sn}/core/{sn}_ivar_variants.tsv', sn=sample_names)

rule breseq:
    input: expand('{sn}/breseq/output/index.html', sn=sample_names)

rule freebayes:
    input: 
        expand('{sn}/freebayes/{sn}.consensus.fasta', sn=sample_names),
        expand('{sn}/freebayes/{sn}.variants.norm.vcf', sn=sample_names),
        'freebayes_lineage_assignments.tsv',
        expand('{sn}/freebayes/quast/{sn}_quast_report.html', sn=sample_names),
        expand('{sn}/freebayes/{sn}_consensus_compare.vcf', sn=sample_names)

    
rule coverage:
    input: expand('{sn}/coverage/{sn}_depth.txt', sn=sample_names)

rule coverage_plot:
    input: expand('{sn}/coverage/{sn}_coverage_plot.png', sn=sample_names)

rule kraken2:
    input: expand('{sn}/kraken2/{sn}_kraken2.out', sn=sample_names)

rule quast:
    input: expand('{sn}/quast/{sn}_quast_report.html', sn=sample_names)

rule lineages:
    input:
        'input_pangolin_versions.txt',
        'input_nextclade_versions.txt',
        'lineage_assignments.tsv'

rule config_sample_log:
    input: 
        config_filename,
        config['samples']

# to handle different options in variant calling
if config['run_breseq'] and config['run_freebayes']:
    if breseq_ref == "":
        print("Invalid BreSeq reference (paramter: breseq_reference) in config file. Please double check and restart")
        exit(1)
    rule variant_calling:
        input:
            rules.breseq.input,
            rules.ivar_variants.input,
            rules.consensus.input,
            rules.freebayes.input
elif config['run_breseq'] and not config['run_freebayes']:
    if breseq_ref == "":
        print("Invalid BreSeq reference (paramter: breseq_reference) in config file. Please double check and restart")
        exit(1)
    rule variant_calling:
        input:
            rules.breseq.input,
            rules.ivar_variants.input,
            rules.consensus.input
elif not config['run_breseq'] and config['run_freebayes']:
    rule variant_calling:
        input:
            rules.freebayes.input,
            rules.ivar_variants.input,
            rules.consensus.input
else:
    rule variant_calling:
        input:
            rules.ivar_variants.input,
            rules.consensus.input

rule all:
    input:
        rules.raw_read_data_symlinks.input,
        rules.host_removed_raw_reads.input,
        rules.remove_adapters.input,
        rules.fastqc.input,
        rules.clean_reads.input,
        rules.coverage.input,
        rules.coverage_plot.input,
        rules.kraken2.input,
        rules.quast.input,
        rules.config_sample_log.input,
        rules.variant_calling.input,
        rules.lineages.input

rule postprocess:
    conda: 
        'conda_envs/postprocessing.yaml'
    params:
        sample_csv_filename = os.path.join(exec_dir, config['samples']),
        postprocess_script_path = os.path.join(exec_dir, 'scripts', 'signal_postprocess.py')
    shell:
        '{params.postprocess_script_path} {params.sample_csv_filename}'


rule ncov_tools:
    conda:
        'ncov-tools/workflow/envs/environment.yml'
    threads: workflow.cores
    params:
        exec_dir = exec_dir,
        sample_csv_filename = os.path.join(exec_dir, config['samples']),
        result_dir = os.path.basename(config['result_dir']),
        amplicon_bed = os.path.join(exec_dir, config['amplicon_loc_bed']),
        primer_bed = os.path.join(exec_dir, config['scheme_bed']),
        viral_reference_genome = os.path.join(exec_dir, config['viral_reference_genome']),
        phylo_include_seqs = os.path.join(exec_dir, config['phylo_include_seqs']),
        negative_control_prefix = config['negative_control_prefix'],
        freebayes_run = config['run_freebayes'],
        pangolin = versions['pangolin'],
        mode = pango_speed
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
        '{sn}/raw_fastq/{sn}_R{r}.fastq.gz'
    input:
        lambda wildcards: get_input_fastq_files(wildcards.sn, wildcards.r)
    shell:
        'ln -s {input} {output}'

rule concat_and_sort:
    priority: 4
    output:
        '{sn}/raw_fastq/{sn}_R{r}.fastq.gz'
    input:
        lambda wildcards: get_pooled_fastq_files(wildcards.sn, wildcards.r)
    benchmark:
        "{sn}/benchmarks/{sn}_concat_and_sort_R{r}.benchmark.tsv"
    shell:
        'if [ $(echo {input} | wc -w) -gt 1 ]; then zcat -f {input} | paste - - - - | sort -k1,1 -t " " | tr "\\t" "\\n" | gzip > {output}; else ln -s {input} {output}; fi'

rule run_raw_fastqc:
    conda: 
        'conda_envs/trim_qc.yaml'
    output:
        r1_fastqc = '{sn}/raw_fastq/{sn}_R1_fastqc.html',
        r2_fastqc = '{sn}/raw_fastq/{sn}_R2_fastqc.html'
    input:
        r1 = '{sn}/raw_fastq/{sn}_R1.fastq.gz',
        r2 = '{sn}/raw_fastq/{sn}_R2.fastq.gz'
    benchmark:
        '{sn}/benchmarks/{sn}_raw_fastqc.benchmark.tsv'
    params:
        output_prefix = '{sn}/raw_fastq'
    log:
        '{sn}/raw_fastq/{sn}_fastqc.log'
    shell:
        """
        fastqc -o {params.output_prefix} {input} 2> {log}
        """

########################## Human Host Removal ################################

rule raw_reads_composite_reference_bwa_map:
    threads: 2
    conda: 
        'conda_envs/snp_mapping.yaml'
    output:
        '{sn}/host_removal/{sn}_viral_and_nonmapping_reads.bam',
    input:
        raw_r1 = '{sn}/raw_fastq/{sn}_R1.fastq.gz',
        raw_r2 = '{sn}/raw_fastq/{sn}_R2.fastq.gz'
    benchmark:
        "{sn}/benchmarks/{sn}_composite_reference_bwa_map.benchmark.tsv"
    log:
        '{sn}/host_removal/{sn}_human_read_mapping.log'
    params:
       composite_index = os.path.join(exec_dir, config['composite_reference']),
       script_path = os.path.join(exec_dir, "scripts", "filter_non_human_reads.py"),
       viral_contig_name = config['viral_reference_contig_name']
    shell:
        '(bwa mem -t {threads} {params.composite_index} '
        '{input.raw_r1} {input.raw_r2} | '
        '{params.script_path} -c {params.viral_contig_name} > {output}) 2> {log}'

rule get_host_removed_reads:
    threads: 2
    conda: 'conda_envs/snp_mapping.yaml'
    output:
        r1 = '{sn}/host_removal/{sn}_R1.fastq.gz',
        r2 = '{sn}/host_removal/{sn}_R2.fastq.gz',
        s = '{sn}/host_removal/{sn}_singletons.fastq.gz',
        bam = '{sn}/host_removal/{sn}_viral_and_nonmapping_reads_filtered_sorted.bam'
    input:
        '{sn}/host_removal/{sn}_viral_and_nonmapping_reads.bam',
    benchmark:
        "{sn}/benchmarks/{sn}_get_host_removed_reads.benchmark.tsv"
    log:
        '{sn}/host_removal/{sn}_samtools_fastq.log'
    shell:
        """
        samtools view -b {input} | samtools sort -n -@{threads} > {output.bam} 2> {log}
        samtools fastq -1 {output.r1} -2 {output.r2} -s {output.s} {output.bam} 2>> {log} 
        """

###### Based on github.com/connor-lab/ncov2019-artic-nf/blob/master/modules/illumina.nf#L124 ######

rule run_trimgalore:
    threads: 2
    priority: 2
    conda: 
        'conda_envs/trim_qc.yaml'
    output:
        '{sn}/adapter_trimmed/{sn}_R1_val_1.fq.gz',
        '{sn}/adapter_trimmed/{sn}_R2_val_2.fq.gz',
        '{sn}/adapter_trimmed/{sn}_R1_val_1_fastqc.html',
        '{sn}/adapter_trimmed/{sn}_R2_val_2_fastqc.html'
    input:
        raw_r1 = '{sn}/host_removal/{sn}_R1.fastq.gz',
        raw_r2 = '{sn}/host_removal/{sn}_R2.fastq.gz'
    log:
        '{sn}/adapter_trimmed/{sn}_trim_galore.log'
    benchmark:
        "{sn}/benchmarks/{sn}_trimgalore.benchmark.tsv"
    params:
        min_len = config['min_len'],
        min_qual = config['min_qual'],
        output_prefix = '{sn}/adapter_trimmed'
    shell:
        'trim_galore --quality {params.min_qual} --length {params.min_len} '
        ' -o {params.output_prefix} --cores {threads} --fastqc '
        '--paired {input.raw_r1} {input.raw_r2} 2> {log} || touch {output}'

rule run_filtering_of_residual_adapters:
    threads: 2
    priority: 2
    conda: 
        'conda_envs/snp_mapping.yaml'
    input:
        r1 = '{sn}/adapter_trimmed/{sn}_R1_val_1.fq.gz',
        r2 = '{sn}/adapter_trimmed/{sn}_R2_val_2.fq.gz'
    output:
        '{sn}/adapter_trimmed/{sn}_R1_val_1_posttrim_filter.fq.gz',
        '{sn}/adapter_trimmed/{sn}_R2_val_2_posttrim_filter.fq.gz'
    params:
        script_path = os.path.join(exec_dir, "scripts", "filter_residual_adapters.py")
    shell:
        """
        python {params.script_path} --input_R1 {input.r1} --input_R2 {input.r2}
        """
       
rule viral_reference_bwa_build:
    conda: 
        'conda_envs/snp_mapping.yaml'
    output:
        '{sn}/core/viral_reference.bwt'
    input:
        reference = os.path.join(exec_dir, config['viral_reference_genome']),
    log:
        '{sn}/core/{sn}_viral_reference_bwa-build.log'
    benchmark:
        "{sn}/benchmarks/{sn}_reference_bwa_build.benchmark.tsv"
    params:
        output_prefix = "{sn}/core/viral_reference"
    shell:
        'bwa index -p {params.output_prefix} {input} >{log} 2>&1'


rule viral_reference_bwa_map:
    threads: 2
    conda: 
        'conda_envs/snp_mapping.yaml'
    output:
        '{sn}/core/{sn}_viral_reference.bam'
    input:
        r1 = '{sn}/adapter_trimmed/{sn}_R1_val_1_posttrim_filter.fq.gz',
        r2 = '{sn}/adapter_trimmed/{sn}_R2_val_2_posttrim_filter.fq.gz',
        ref = '{sn}/core/viral_reference.bwt'
    benchmark:
        "{sn}/benchmarks/{sn}_viral_reference_bwa_map.benchmark.tsv"
    log:
        '{sn}/core/{sn}_viral_reference_bwa.log'
    params:
       ref_prefix = '{sn}/core/viral_reference'
    shell:
        '(bwa mem -t {threads} {params.ref_prefix} '
        '{input.r1} {input.r2} | '
        'samtools view -bS | samtools sort -@{threads} -o {output}) 2> {log}'


rule run_bed_primer_trim:
    conda: 
        'conda_envs/ivar.yaml'
    input:
        "{sn}/core/{sn}_viral_reference.bam"
    output:
        sorted_trimmed_mapped_bam = "{sn}/core/{sn}_viral_reference.mapping.primertrimmed.sorted.bam",
        trimmed_mapped_bam = "{sn}/core/{sn}_viral_reference.mapping.primertrimmed.bam",
        mapped_bam = "{sn}/core/{sn}_viral_reference.mapping.bam"
    benchmark:
        "{sn}/benchmarks/{sn}_bed_primer_trim.benchmark.tsv"
    log:
        "{sn}/core/{sn}_ivar_trim.log"
    params:
        scheme_bed = os.path.join(exec_dir, config['scheme_bed']),
        ivar_output_prefix = "{sn}/core/{sn}_viral_reference.mapping.primertrimmed",
        min_len = config['min_len'],
        min_qual = config['min_qual'],
        primer_pairs = config['primer_pairs_tsv']
    shell:
        'samtools view -F4 -o {output.mapped_bam} {input}; '
        'samtools index {output.mapped_bam}; '
        'ivar trim -e -i {output.mapped_bam} -b {params.scheme_bed} '
        '-m {params.min_len} -q {params.min_qual} '
        '{params.primer_pairs} '
        '-p {params.ivar_output_prefix} 2> {log}; '
        'samtools sort -o {output.sorted_trimmed_mapped_bam} '
        '{output.trimmed_mapped_bam}'


rule run_fastqc_on_mapped_reads:
    conda: 'conda_envs/trim_qc.yaml'
    output:
        r1_fastqc = '{sn}/mapped_clean_reads/{sn}_R1_fastqc.html',
        r2_fastqc = '{sn}/mapped_clean_reads/{sn}_R2_fastqc.html'
    input:
        r1 = '{sn}/mapped_clean_reads/{sn}_R1.fastq.gz',
        r2 = '{sn}/mapped_clean_reads/{sn}_R2.fastq.gz'
    benchmark:
        '{sn}/benchmarks/{sn}_clean_fastqc.benchmark.tsv'
    params:
        output_prefix = '{sn}/mapped_clean_reads'
    log:
        '{sn}/mapped_clean_reads/{sn}_fastqc.log'
    shell:
        """
        fastqc -o {params.output_prefix} {input} 2> {log}
        """

rule get_mapping_reads:
    priority: 2
    conda: 'conda_envs/snp_mapping.yaml'
    output:
        r1 = '{sn}/mapped_clean_reads/{sn}_R1.fastq.gz',
        r2 = '{sn}/mapped_clean_reads/{sn}_R2.fastq.gz',
        s = '{sn}/mapped_clean_reads/{sn}_singletons.fastq.gz',
        bam = '{sn}/mapped_clean_reads/{sn}_sorted_clean.bam'
    input:
        "{sn}/core/{sn}_viral_reference.mapping.primertrimmed.bam",
    benchmark:
        "{sn}/benchmarks/{sn}_get_mapping_reads.benchmark.tsv"
    log:
        '{sn}/mapped_clean_reads/{sn}_samtools_fastq.log'
    shell:
        """
        samtools sort -n {input} -o {output.bam} 2> {log}
        samtools fastq -1 {output.r1} -2 {output.r2} -s {output.s} {output.bam} 2>> {log} 
        """

rule run_ivar_consensus:
    conda: 
        'conda_envs/ivar.yaml'
    output:
        '{sn}/core/{sn}.consensus.fa'
    input:
        "{sn}/core/{sn}_viral_reference.mapping.primertrimmed.sorted.bam"
    log:
        '{sn}/core/{sn}_ivar_consensus.log'
    benchmark:
        "{sn}/benchmarks/{sn}_ivar_consensus.benchmark.tsv"
    params:
        mpileup_depth = config['mpileup_depth'],
        ivar_min_coverage_depth = config['var_min_coverage_depth'],
        ivar_freq_threshold = config['var_freq_threshold'],
        output_prefix = '{sn}/core/{sn}.consensus',
    shell:
        '(samtools mpileup -aa -A -d {params.mpileup_depth} -Q0 {input} | '
        'ivar consensus -t {params.ivar_freq_threshold} '
        '-m {params.ivar_min_coverage_depth} -n N -p {params.output_prefix}) '
        '2>{log}'

rule index_viral_reference:
    # from @jts both mpileup and ivar need a reference .fai file and will create 
    # it when it doesn't exist. 
    # When they're run in a pipe like mpileup | ivar there's a race condition 
    # that causes the error
    conda: 
        'conda_envs/ivar.yaml'
    output:
        os.path.join(exec_dir, config['viral_reference_genome']) + ".fai"
    input:
        os.path.join(exec_dir, config['viral_reference_genome']),
    shell:
        'samtools faidx {input}'


rule run_ivar_variants:
    conda: 
        'conda_envs/ivar.yaml'
    output:
        '{sn}/core/{sn}_ivar_variants.tsv'
    input:
        reference = os.path.join(exec_dir, config['viral_reference_genome']),
        indexed_reference = os.path.join(exec_dir, config['viral_reference_genome']) + ".fai",
        read_bam = "{sn}/core/{sn}_viral_reference.mapping.primertrimmed.sorted.bam",
        viral_reference_gff = os.path.join(exec_dir, config['viral_reference_feature_coords'])
    log:
        '{sn}/core/{sn}_ivar_variants.log'
    benchmark:
        "{sn}/benchmarks/{sn}_ivar_variants.benchmark.tsv"
    params:
        output_prefix = '{sn}/core/{sn}_ivar_variants',
        ivar_min_coverage_depth = config['var_min_coverage_depth'],
        ivar_min_freq_threshold = config['var_min_freq_threshold'],
        ivar_min_variant_quality = config['var_min_variant_quality'],
    shell:
        '(samtools mpileup -aa -A -d 0 --reference {input.reference} -B '
            '-Q 0 {input.read_bam} | '
        'ivar variants -r {input.reference} -m {params.ivar_min_coverage_depth} '
        '-p {params.output_prefix} -q {params.ivar_min_variant_quality} '
        '-t {params.ivar_min_freq_threshold} -g {input.viral_reference_gff}) 2> {log}'


################################   Based on scripts/breseq.sh   ####################################

rule run_breseq:
    threads: 4
    priority: 1
    conda: 'conda_envs/snp_mapping.yaml'
    output:
        '{sn}/breseq/output/index.html'
    input:
        expand('{{sn}}/mapped_clean_reads/{{sn}}_R{r}.fastq.gz', r=[1,2])
    log:
        '{sn}/breseq/{sn}_breseq.log',
    benchmark:
        "{sn}/benchmarks/{sn}_run_breseq.benchmark.tsv"
    params:
        ref = os.path.join(exec_dir, breseq_ref),
        outdir = '{sn}/breseq'
    shell:
        """
        breseq --reference {params.ref} --num-processors {threads} --polymorphism-prediction --brief-html-output --output {params.outdir} {input} > {log} 2>&1 || touch {output}
        """

################## Based on https://github.com/jts/ncov2019-artic-nf/blob/be26baedcc6876a798a599071bb25e0973261861/modules/illumina.nf ##################

rule run_freebayes:
    threads: 1
    priority: 1
    conda: 'conda_envs/freebayes.yaml'
    output:
        consensus = '{sn}/freebayes/{sn}.consensus.fasta',
        variants = '{sn}/freebayes/{sn}.variants.norm.vcf'
    input:
        reference = os.path.join(exec_dir, config['viral_reference_genome']),
        read_bam = "{sn}/core/{sn}_viral_reference.mapping.primertrimmed.sorted.bam"
    params:
        out = '{sn}/freebayes/work/{sn}',
        freebayes_min_coverage_depth = config['var_min_coverage_depth'],
        freebayes_min_freq_threshold = config['var_min_freq_threshold'],
        freebayes_min_variant_quality = config['var_min_variant_quality'],
        freebayes_freq_threshold = config['var_freq_threshold'],
        script_path = os.path.join(exec_dir, "scripts", "process_gvcf.py")
    shell:
        """
        mkdir -p $(dirname {params.out})
        # the sed is to fix the header until a release is made with this fix
        # https://github.com/freebayes/freebayes/pull/549
        freebayes -p 1 \
                  -f {input.reference} \
                  -F 0.2 \
                  -C 1 \
                  --pooled-continuous \
                  --min-coverage {params.freebayes_min_coverage_depth} \
                  --gvcf --gvcf-dont-use-chunk true {input.read_bam} | sed s/QR,Number=1,Type=Integer/QR,Number=1,Type=Float/ > {params.out}.gvcf

        # make depth mask, split variants into ambiguous/consensus
        # NB: this has to happen before bcftools norm or else the depth mask misses any bases exposed during normalization
        python {params.script_path} -d {params.freebayes_min_coverage_depth} \
                        -l {params.freebayes_min_freq_threshold} \
                        -u {params.freebayes_freq_threshold} \
                        -m {params.out}.mask.txt \
                        -v {params.out}.variants.vcf \
                        -c {params.out}.consensus.vcf {params.out}.gvcf 

        # normalize variant records into canonical VCF representation
        bcftools norm -f {input.reference} {params.out}.variants.vcf > {output.variants} 
        bcftools norm -f {input.reference} {params.out}.consensus.vcf > {params.out}.consensus.norm.vcf

        # split the consensus sites file into a set that should be IUPAC codes and all other bases, using the ConsensusTag in the VCF
        for vt in "ambiguous" "fixed"; do
            cat {params.out}.consensus.norm.vcf | awk -v vartag=ConsensusTag=$vt '$0 ~ /^#/ || $0 ~ vartag' > {params.out}.$vt.norm.vcf
            bgzip -f {params.out}.$vt.norm.vcf  
            tabix -f -p vcf {params.out}.$vt.norm.vcf.gz 
        done
        
        # apply ambiguous variants first using IUPAC codes. this variant set cannot contain indels or the subsequent step will break
        bcftools consensus -f {input.reference} -I {params.out}.ambiguous.norm.vcf.gz > {params.out}.ambiguous.fasta 
        # apply remaninng variants, including indels
        bcftools consensus -f {params.out}.ambiguous.fasta -m {params.out}.mask.txt {params.out}.fixed.norm.vcf.gz | sed s/MN908947\.3.*/{wildcards.sn}/ > {output.consensus}
        """

rule consensus_compare:
    threads: 1
    priority: 1
    conda: 'conda_envs/freebayes.yaml'
    output:
        '{sn}/freebayes/{sn}_consensus_compare.vcf'
    input:
        ivar = '{sn}/core/{sn}.consensus.fa',
        freebayes = '{sn}/freebayes/{sn}.consensus.fasta'
    params:
        script_path = os.path.join(exec_dir, "scripts", "quick_align.py")
    shell:
        """
        python {params.script_path} -g {input.ivar} -r {input.freebayes} -o vcf > {output}
        """

##################  Based on scripts/hisat2.sh and scripts/coverage_stats_avg.sh  ##################

rule coverage_depth:
    conda: 'conda_envs/snp_mapping.yaml'
    output:
        '{sn}/coverage/{sn}_depth.txt'
    input:
        "{sn}/core/{sn}_viral_reference.mapping.primertrimmed.sorted.bam"
    benchmark:
        "{sn}/benchmarks/{sn}_coverage_depth.benchmark.tsv"
    shell:
        'bedtools genomecov -d -ibam {input} > {output}'

rule generate_coverage_plot:
    conda: 'conda_envs/postprocessing.yaml'
    output: 
        '{sn}/coverage/{sn}_coverage_plot.png' 
    input:
        '{sn}/coverage/{sn}_depth.txt'
    params:
        script_path = os.path.join(exec_dir, "scripts", "generate_coverage_plot.py")
    shell:
        "python {params.script_path} {input} {output}"

################################   Based on scripts/kraken2.sh   ###################################

rule run_kraken2:
    threads: 1
    conda: 'conda_envs/trim_qc.yaml'
    output:
        '{sn}/kraken2/{sn}_kraken2.out'
    input:
        r1 = '{sn}/adapter_trimmed/{sn}_R1_val_1_posttrim_filter.fq.gz',
        r2 = '{sn}/adapter_trimmed/{sn}_R2_val_2_posttrim_filter.fq.gz'
    log:
        '{sn}/kraken2/{sn}_kraken2.log'
    benchmark:
        "{sn}/benchmarks/{sn}_run_kraken2.benchmark.tsv"
    params:
        outdir = '{sn}/kraken2',
        db = os.path.join(exec_dir, config['kraken2_db']),
        labelled_output = '{sn}_kraken2.out',
        labelled_report = '{sn}_kraken2.report',
        labelled_unclassified_reads = '{sn}_kraken2_unclassified_reads#',
        labelled_classified_reads = '{sn}_kraken2_classified_reads#'
    shell:
        'cd {params.outdir} '
        '&& kraken2'
            ' --db {params.db}'
            ' --threads {threads}'
            ' --quick --unclassified-out "{params.labelled_unclassified_reads}"'
            ' --classified-out "{params.labelled_classified_reads}"'
            ' --output {params.labelled_output}'
            ' --paired --gzip-compressed'
            ' ../../{input.r1} ../../{input.r2}'
            ' --report {params.labelled_report}'
            ' 2>../../{log} && (cd ../.. && touch {output})'
            # kraken2 also fails if empty input is provided which will happen
            # if there are no valid reads e.g., very clean negative control


##################################  Based on scripts/quast.sh   ####################################


rule run_quast:
    threads: 1
    conda: 'conda_envs/assembly_qc.yaml'
    output:
         '{sn}/quast/{sn}_quast_report.html'
    input:
         '{sn}/core/{sn}.consensus.fa'
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

rule run_quast_freebayes:
    threads: 1
    conda: 'conda_envs/assembly_qc.yaml'
    output:
         '{sn}/freebayes/quast/{sn}_quast_report.html'
    input:
         '{sn}/freebayes/{sn}.consensus.fasta'
    log:
         '{sn}/freebayes/quast/{sn}_quast.log'
    benchmark:
        "{sn}/benchmarks/{sn}_run_quast.benchmark.tsv"
    params:
         outdir = '{sn}/freebayes/quast',
         genome = os.path.join(exec_dir, config['viral_reference_genome']),
         fcoords = os.path.join(exec_dir, config['viral_reference_feature_coords']),
         sample_name = '{sn}_quast_report',
         unlabelled_reports = '{sn}/freebayes/quast/report.*'
    shell:
         'quast {input} -r {params.genome} -g {params.fcoords} --output-dir {params.outdir} --threads {threads} >{log} && '
         'for f in {params.unlabelled_reports}; do mv $f ${{f/report/{params.sample_name}}}; done'

rule run_lineage_assignment:
    threads: 4
    conda: 'conda_envs/assign_lineages.yaml'
    output:
        pango_ver_out = 'input_pangolin_versions.txt',
        nextclade_ver_out = 'input_nextclade_versions.txt',
        lin_out = 'lineage_assignments.tsv'
    input:
        expand('{sn}/core/{sn}.consensus.fa', sn=sample_names)
    params:
        pangolin_ver = versions['pangolin'],
        pangolearn_ver = versions['pangolearn'],
        constellations_ver = versions['constellations'],
        scorpio_ver = versions['scorpio'],
        designation_ver = versions['pango-designation'],
        data_ver = versions['pangolin-data'],
        #accession = config['viral_reference_contig_name'],
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

rule run_lineage_assignment_freebayes:
    threads: 4
    conda: 'conda_envs/assign_lineages.yaml'
    output:
        'freebayes_lineage_assignments.tsv'
    input:
        p_vers = 'input_pangolin_versions.txt',
        n_vers = 'input_nextclade_versions.txt',
        consensus = expand('{sn}/freebayes/{sn}.consensus.fasta', sn=sample_names)
    params:
        analysis_mode = pango_speed,
        assignment_script_path = os.path.join(exec_dir, 'scripts', 'assign_lineages.py')
    shell:
        'cat {input.consensus} > all_freebayes_genomes.fa && '
        '{params.assignment_script_path} -i all_freebayes_genomes.fa -t {threads} -o {output} -p {input.p_vers} -n {input.n_vers} --mode {params.analysis_mode} --skip'

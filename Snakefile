from snakemake.utils import validate
import pandas as pd

# Note: this workflow uses specialized conda envs, run with 'snakemake --use-conda'.

# TODO: needs final debug iteration, after the dust settles

# The config file contains a high-level summary of pipeline configuration and inputs.
# It is ingested by the Snakefile, and also intended to be human-readable.
# For an example config file, see pipeline/example_config.yaml in the covid-19-sequencing repo.

# read and validate config.yaml
configfile: "config.yaml"
validate(config, 'resources/config.schema.yaml')

# read and validate sample table specified in config.schema.yaml
samples = pd.read_table(config['samples'], sep=',')
validate(samples, 'resources/sample.schema.yaml')

# get sample names 
sample_names = sorted(samples['sample'].drop_duplicates().values)

def get_input_fastq_files(sample_name, r):
    sample_fastqs = samples[samples['sample'] == sample_name]
    if r == '1':
        relpaths = sample_fastqs['r1_path'].values
    elif r == '2':
        relpaths = sample_fastqs['r2_path'].values

    return [ os.path.join(sample_name, r) for r in relpaths ]


######################################   High-level targets   ######################################


rule all:
    input: 'summary.html'

rule sort:
    input: expand('{sn}/fastq_sorted/R{r}_fastqc.html', sn=sample_names, r=[1,2])
    
rule remove_primers:
    input: expand('{sn}/fastq_primers_removed/R{r}.fastq.gz', sn=sample_names, r=[1,2])

rule trim:
    input: expand('{sn}/fastq_trimmed/R{r}_{s}_fastqc.html', sn=sample_names, r=[1,2], s=['paired','unpaired'])

rule hostremove:
    input:
       expand('{sn}/host_removed/both_ends_mapped_lsorted.bam', sn=sample_names),
       expand('{sn}/host_removed/R{r}.fastq.gz', sn=sample_names, r=[1,2])

rule consensus:
    input: expand('{sn}/consensus/virus.consensus.fa', sn=sample_names)

rule ivar_variants:
    input: expand('{sn}/ivar_variants/ivar_variants.tsv', sn=sample_names)

rule breseq:
    input: expand('{sn}/breseq/output/index.html', sn=sample_names)
    
rule coverage:
    input: expand('{sn}/coverage/depth.txt', sn=sample_names)

rule kraken2:
    input: expand('{sn}/kraken2/kraken2.out', sn=sample_names)

rule lmat:
    input: expand('{sn}/lmat/parseLMAT_output.txt', sn=sample_names)

rule quast:
    input: expand('{sn}/quast/report.html', sn=sample_names)


#################################   Based on scripts/assemble.sh   #################################


rule concat_and_sort:
    priority: 4
    output:
        '{sn}/fastq_sorted/R{r}.fastq.gz'
    input:
        lambda wildcards: get_input_fastq_files(wildcards.sn, wildcards.r)
    shell:
        'zcat {input} | paste - - - - | sort -k1,1 -t " " | tr "\\t" "\\n" | gzip > {output}'


rule run_fastqc:
    conda: 'conda_envs/trim_qc.yaml'
    output: expand('{{s}}_fastqc.{ext}', ext=['html','zip'])
    input: '{s}.fastq.gz'
    log: '{s}_fastqc.log',
    shell: 'fastqc {input} 2>{log}'


rule run_cutadapt:
    threads: 16
    priority: 3
    conda: 'conda_envs/trim_qc.yaml'
    output:
        expand('{{sn}}/fastq_primers_removed/R{r}.fastq.gz', r=[1,2])
    input:
        expand('{{sn}}/fastq_sorted/R{r}.fastq.gz', r=[1,2])
    log:
        '{sn}/fastq_primers_removed/cutadapt.log'
    params:
        primer_fw = config['primer_fw'],
        primer_rc = config['primer_rc']
    shell:
        'cutadapt -j {threads}'
	' -a file:{params.primer_rc} -A file:{params.primer_fw}'   # primers
	' -o {output[0]} -p {output[1]}'     # output files
	' {input}'                           # input files
        ' >{log}'                            # log file

# Note: expand()-statements in 'output:' and 'input:' have been written so that the ordering
# of their outputs is consistent with the ordering of trimmomatic's command-line arguments.

rule run_trimmomatic:
    threads: 2
    priority: 2
    conda: 'conda_envs/trim_qc.yaml'
    output:
        expand('{{sn}}/fastq_trimmed/R{r}_{s}.fastq.gz', r=[1,2], s=['paired','unpaired'])
    input:
        expand('{{sn}}/fastq_primers_removed/R{r}.fastq.gz', r=[1,2])
    log:
        '{sn}/fastq_trimmed/trim.log'
    params:
        targs = config['trimmomatic_args']
    shell:
        'trimmomatic PE -threads {threads} {input} {output} {params.targs} 2>{log}'


############################  Based on scripts/remove_host_sequences.sh  ###########################


rule hostremove_hisat2_build:
    priority: 2
    conda: 'conda_envs/snp_mapping.yaml'
    output:
        '{sn}/host_removed/sars-cov-2.1.ht2'
    log:
        '{sn}/host_removed/hisat2-build.log'
    params:
        reference = config['viral_reference_genome'],
	genome = '{sn}/host_removed/sars-cov-2'
    shell:
        'hisat2-build {params.reference} {params.genome} >{log} 2>&1'


rule hostremove_hisat2:
    threads: 2
    priority: 2
    conda: 'conda_envs/snp_mapping.yaml'
    output:
        '{sn}/host_removed/mapped_and_unmapped.sam'
    input:
        '{sn}/fastq_trimmed/R1_paired.fastq.gz',
        '{sn}/fastq_trimmed/R2_paired.fastq.gz',
        '{sn}/host_removed/sars-cov-2.1.ht2'
    log:
        '{sn}/host_removed/hisat2.log'
    params:
        genome = '{sn}/host_removed/sars-cov-2',
	summary_file = '{sn}/fastq_hist_removed/hisat2_summary.txt'
    shell:
        'hisat2 --threads {threads}'
        ' -x {params.genome}'
	' -1 {input[0]} -2 {input[1]}'
	' --summary-file {params.summary_file}'
	' -S {output}'
	' 2>{log}'


rule hostremove_sam_to_bam:
    priority: 2
    conda: 'conda_envs/snp_mapping.yaml'
    output:
        '{sn}/host_removed/mapped_and_unmapped.bam'
    input:
        '{sn}/host_removed/mapped_and_unmapped.sam'
    shell:
        'samtools view -bS {input} > {output}'


rule hostremove_map_pairs:
    priority: 2
    conda: 'conda_envs/snp_mapping.yaml'
    output:
        '{sn}/host_removed/both_ends_mapped.bam'
    input:
        '{sn}/host_removed/mapped_and_unmapped.bam'
    shell:
        'samtools view -b -f 3 -F 4 {input} > {output}'


rule hostremove_lsort:
    conda: 'conda_envs/snp_mapping.yaml'
    output:
        '{sn}/host_removed/both_ends_mapped_lsorted.bam'
    input:
        '{sn}/host_removed/both_ends_mapped.bam'
    shell:
        'samtools sort {input} -o {output}'


rule hostremove_nsort:
    priority: 2
    conda: 'conda_envs/snp_mapping.yaml'
    output:
        '{sn}/host_removed/both_ends_mapped_nsorted.bam'
    input:
        '{sn}/host_removed/both_ends_mapped.bam'
    shell:
        'samtools sort -n {input} -o {output}'


rule hostremove_fastq:
    priority: 2
    conda: 'conda_envs/snp_mapping.yaml'
    output:
        '{sn}/host_removed/R1.fastq',
        '{sn}/host_removed/R2.fastq'
    input:
        '{sn}/host_removed/both_ends_mapped_nsorted.bam'
    shell:
        'bedtools bamtofastq -i {input} -fq {output[0]} -fq2 {output[1]}'


rule hostremove_fastq_gzip:
    priority: 2
    output:
        '{sn}/host_removed/R{r}.fastq.gz'
    input:
        '{sn}/host_removed/R{r}.fastq'
    shell:
        'gzip {input}'


###### Based on github.com/connor-lab/ncov2019-artic-nf/blob/master/modules/illumina.nf#L124 ######


rule run_consensus:
    conda: 'conda_envs/ivar.yaml'
    output:
        '{sn}/consensus/virus.consensus.fa'
    input:
        '{sn}/host_removed/both_ends_mapped_lsorted.bam'
    log:
        '{sn}/consensus/ivar.log'
    params:
        mpileup_depth = config['mpileup_depth'],
        ivar_min_coverage_depth = config['ivar_min_coverage_depth'],
        ivar_freq_threshold = config['ivar_freq_threshold'],
        output_prefix = '{sn}/consensus/virus.consensus'
    shell:
        """
        samtools mpileup -A -d {params.mpileup_depth} -Q0 {input} | \
        ivar consensus -t {params.ivar_freq_threshold} -m {params.ivar_min_coverage_depth} \
                                    -n N -p {params.output_prefix} 2>{log}
        """


rule run_ivar_variants:
    conda: 'conda_envs/ivar.yaml'
    output:
        '{sn}/ivar_variants/ivar_variants.tsv'
    input:
        reference = config['viral_reference_genome'],
        read_bam = '{sn}/host_removed/both_ends_mapped_lsorted.bam'
    params:
        output_prefix = '{sn}/ivar_variants/ivar_variants',
        ivar_min_coverage_depth = config['ivar_min_coverage_depth'],
        ivar_min_freq_threshold = config['ivar_min_freq_threshold'],
        ivar_min_variant_quality = config['ivar_min_variant_quality']
    shell:
        """
        samtools mpileup -A -d 0 --reference {input.reference} -B -Q 0 {input.read_bam} |\
        ivar variants -r {input.reference} -m {params.ivar_min_coverage_depth} -p {params.output_prefix} -q {params.ivar_min_variant_quality} -t {params.ivar_min_freq_threshold}
        """


################################   Based on scripts/breseq.sh   ####################################


rule run_breseq:
    threads: 1
    priority: 1
    conda: 'conda_envs/snp_mapping.yaml'
    output:
        '{sn}/breseq/output/index.html'
    input:
        expand('{{sn}}/host_removed/R{r}.fastq.gz', r=[1,2])
    log:
        '{sn}/breseq/breseq.log',
    params:
        ref = config['breseq_reference'],
	outdir = '{sn}/breseq'
    shell:
        'breseq --reference {params.ref} --num-processors {threads} --polymorphism-prediction --brief-html-output --output {params.outdir} {input} >{log} 2>&1'


##################  Based on scripts/hisat2.sh and scripts/coverage_stats_avg.sh  ##################


rule coverage_hisat2_build:
    conda: 'conda_envs/snp_mapping.yaml'
    output:
        '{sn}/coverage/genome.1.ht2'
    input:
        '{sn}/consensus/virus.consensus.fa'
    log:
        '{sn}/coverage/hisat2-build.log'
    params:
        genome = '{sn}/coverage/genome'
    shell:
        'hisat2-build {input} {params.genome} >{log} 2>&1'


rule coverage_hisat2:
    threads: 2
    conda: 'conda_envs/snp_mapping.yaml'
    output:
        '{sn}/coverage/output.sam'
    input:
        '{sn}/host_removed/R1.fastq.gz',
        '{sn}/host_removed/R2.fastq.gz',
        '{sn}/coverage/genome.1.ht2'
    log:
        '{sn}/coverage/hisat2.log'
    params:
        genome = '{sn}/coverage/genome',
        summary_file = '{sn}/coverage/hisat2_summary.txt'
    shell:
        'hisat2 --threads {threads}'
        ' -x {params.genome}'
	' -1 {input[0]} -2 {input[1]}'
	' --summary-file {params.summary_file}'
	' -S {output}'
	' 2>{log}'


rule coverage_sam_to_bam:
    conda: 'conda_envs/snp_mapping.yaml'
    output:
        '{sn}/coverage/output.bam'
    input:
        '{sn}/coverage/output.sam'
    shell:
        'samtools view -b {input} | samtools sort > {output}'


rule coverage_depth:
    conda: 'conda_envs/snp_mapping.yaml'
    output:
        '{sn}/coverage/depth.txt'
    input:
        '{sn}/coverage/output.bam'
    shell:
        'bedtools genomecov -d -ibam {input} >{output}'


################################   Based on scripts/kraken2.sh   ###################################


rule run_kraken2:
    threads: 1
    conda: 'conda_envs/trim_qc.yaml'
    output:
        '{sn}/kraken2/kraken2.out'
    input:
        expand('{{sn}}/fastq_trimmed/R{r}_paired.fastq.gz', r=[1,2])
    log:
        '{sn}/kraken2/kraken2.log'
    params:
        outdir = '{sn}/kraken2',
	db = config['kraken2_db']
    shell:
        'cd {params.outdir} '
	'&& kraken2'
	' --db {params.db}'
	' --threads {threads}'
	' --quick --unclassified-out unclassified-sequences# --classified-out classified-sequences#'
	' --output kraken2.out'
	' --paired --gzip-compressed'
	' ../../{input[0]} ../../{input[1]}'
	' --report report'
        ' 2>../../{log}'


##################################   Based on scripts/lmat.sh   ####################################


rule lmat_pretile:
    output:
        '{sn}/lmat/assembly.tiled.fasta'
    input:
        '{sn}/consensus/virus.consensus.fa'
    params:
        fsize = config['lmat_fragment_size']
    shell:
        'perl scripts/fatile {input} {params.fsize} > {output}'


rule run_lmat:
    threads: 1
    conda: 'conda_envs/lmat.yaml'
    output:
        '{sn}/lmat/assembly.tiled.fasta.{db}.lo.rl_output0.out'
    input:
        '{sn}/lmat/assembly.tiled.fasta'
    params:
        outdir = '{sn}/lmat',
        lmat_basedir = config['lmat_basedir'],
        lmat_db = '{db}'
    shell:
        'LMAT_DIR={params.lmat_basedir}/runtime_inputs '
	'bash $(which run_rl.sh)'
	' --db_file={params.lmat_basedir}/data/{params.lmat_db}'
	' --query_file={input}'
	' --odir={params.outdir}'
	' --overwrite --verbose'
	' --threads={threads}'


rule lmat_postprocess:
    output:
        '{sn}/lmat/parseLMAT_output.txt'
    input:
        expand('{{sn}}/lmat/assembly.tiled.fasta.{db}.lo.rl_output0.out', db=[config['lmat_db']])
    params:
        outdir = '{sn}/lmat'
    shell:
        'cd {params.outdir} && perl ../../scripts/parseLMAT > parseLMAT_output.txt'


##################################  Based on scripts/quast.sh   ####################################


rule run_quast:
    threads: 1
    conda: 'conda_envs/assembly_qc.yaml'
    output:
         '{sn}/quast/report.html'
    input:
         '{sn}/consensus/virus.consensus.fa'
    log:
         '{sn}/quast/quast.log'
    params:
         outdir = '{sn}/quast',
         genome = config['viral_reference_genome'],
         fcoords = config['viral_reference_feature_coords']
    shell:
         'quast {input} -r {params.genome} -g {params.fcoords} --output-dir {params.outdir} --threads {threads} >{log}'


########################################   Postprocessing   ########################################


rule postprocess:
    output: 'summary.html'
    input:
        rules.sort.input,
        rules.remove_primers.input,
        rules.trim.input,
        rules.hostremove.input,
        rules.consensus.input,
        rules.ivar_variants.input,
        rules.breseq.input,
        rules.coverage.input,
        rules.kraken2.input,
        rules.lmat.input,
        rules.quast.input
    params:
        sample_csv_filename = config['samples']
    shell:
        'scripts/c19_postprocess.py {params.sample_csv_filename}'

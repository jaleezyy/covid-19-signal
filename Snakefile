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
import os

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


rule sort:
    input: expand('{sn}/fastq_sorted/R{r}_fastqc.html', sn=sample_names, r=[1,2])
    
rule remove_primers:
    input: expand('{sn}/fastq_sequencing_adapter_trimming/R{r}_paired.fastq.gz', sn=sample_names, r=[1,2])

rule trim:
    input: expand('{sn}/fastq_trimmed/R{r}_paired_fastqc.html', sn=sample_names, r=[1,2])

rule hostremove:
    input:
       expand('{sn}/host_removed/both_ends_mapped_lsorted.bam', sn=sample_names),
       expand('{sn}/host_removed/R{r}.fastq.gz', sn=sample_names, r=[1,2])

rule consensus:
    input: expand('{sn}/core/virus.consensus.fa', sn=sample_names)

rule ivar_variants:
    input: expand('{sn}/core/ivar_variants.tsv', sn=sample_names)

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


rule all:
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


rule postprocess:
    params:
        sample_csv_filename = config['samples']
    shell:
        'scripts/c19_postprocess.py {params.sample_csv_filename}'


#################################   Based on scripts/assemble.sh   #################################


rule concat_and_sort:
    priority: 4
    output:
        '{sn}/fastq_sorted/R{r}.fastq.gz'
    input:
        lambda wildcards: get_input_fastq_files(wildcards.sn, wildcards.r)
    benchmark:
        "{sn}/benchmarks/concat_and_sort_R{r}.benchmark.tsv"
    shell:
        'zcat {input} | paste - - - - | sort -k1,1 -t " " | tr "\\t" "\\n" | gzip > {output}'


rule run_fastqc:
    conda: 'conda_envs/trim_qc.yaml'
    output: expand('{{s}}_fastqc.{ext}', ext=['html','zip'])
    input: '{s}.fastq.gz'
    benchmark:
        "{s}/benchmarks/fastqc.benchmark.tsv"
    log: '{s}_fastqc.log',
    shell: 'fastqc {input} 2>{log} > /dev/null'

# Note: expand()-statements in 'output:' and 'input:' have been written so that the ordering
# of their outputs is consistent with the ordering of trimmomatic's command-line arguments.
rule run_trimmomatic:
    threads: 2
    priority: 2
    conda: 'conda_envs/trim_qc.yaml'
    output:
        expand('{{sn}}/fastq_sequencing_adapter_trimming/R{r}_{s}.fastq.gz', r=[1,2], s=['paired','unpaired'])
    input:
        expand('{{sn}}/fastq_sorted/R{r}.fastq.gz', r=[1,2])
    log:
        '{sn}/fastq_sequencing_adapter_trimming/trim.log'
    benchmark:
        "{sn}/benchmarks/trimmomatic.benchmark.tsv"
    params:
        targs = config['trimmomatic_args']
    shell:
        'trimmomatic PE -threads {threads} {input} {output} {params.targs} 2>{log}'


rule run_cutadapt:
    threads: 16
    priority: 3
    conda: 'conda_envs/trim_qc.yaml'
    output:
        expand('{{sn}}/fastq_trimmed/R{r}_paired.fastq.gz', r=[1,2])
    input:
        expand('{{sn}}/fastq_sequencing_adapter_trimming/R{r}_paired.fastq.gz', r=[1,2])
    log:
        '{sn}/fastq_trimmed/cutadapt.log'
    benchmark:
        "{sn}/benchmarks/cutadapt.benchmark.tsv"
    params:
        primer_fw = config['primer_fw'],
        primer_rc = config['primer_rc']
    shell:
        'cutadapt -j {threads}'
	' -a file:{params.primer_rc} -A file:{params.primer_fw}'   # primers
	' -o {output[0]} -p {output[1]}'     # output files
	' {input}'                           # input files
        ' >{log}'                            # log file


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
    benchmark:
        "{sn}/benchmarks/hostremove_hisat2_build.benchmark.tsv"
    shell:
        'hisat2-build {params.reference} {params.genome} >{log} 2>&1'


rule hostremove_hisat2:
    threads: 2
    priority: 2
    conda: 'conda_envs/snp_mapping.yaml'
    output:
        '{sn}/host_removed/both_ends_mapped.bam'
    input:
        '{sn}/fastq_trimmed/R1_paired.fastq.gz',
        '{sn}/fastq_trimmed/R2_paired.fastq.gz',
        '{sn}/host_removed/sars-cov-2.1.ht2'
    log:
        '{sn}/host_removed/hisat2.log'
    benchmark:
        "{sn}/benchmarks/hostremove_hisat2.benchmark.tsv"
    params:
        genome = '{sn}/host_removed/sars-cov-2',
	summary_file = '{sn}/fastq_hist_removed/hisat2_summary.txt'
    shell:
        'hisat2 --threads {threads}'
        ' -x {params.genome}'
	    ' -1 {input[0]} -2 {input[1]}'
	    ' --summary-file {params.summary_file} 2>{log} | '
        'samtools view -bS | samtools view -b -f 3 -F 4 > {output}'


rule hostremove_lsort:
    conda: 'conda_envs/snp_mapping.yaml'
    output:
        '{sn}/host_removed/both_ends_mapped_lsorted.bam'
    input:
        '{sn}/host_removed/both_ends_mapped.bam'
    benchmark:
        "{sn}/benchmarks/hostremove_lsort.benchmark.tsv"
    shell:
        'samtools sort {input} -o {output}'


rule hostremove_nsort:
    priority: 2
    conda: 'conda_envs/snp_mapping.yaml'
    output:
        '{sn}/host_removed/both_ends_mapped_nsorted.bam'
    input:
        '{sn}/host_removed/both_ends_mapped.bam'
    benchmark:
        "{sn}/benchmarks/hostremove_nsort.benchmark.tsv"
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
    benchmark:
        "{sn}/benchmarks/hostremove_fastq.benchmark.tsv"
    log:
        '{sn}/host_removed/bamtofastq.log'
    shell:
        'bedtools bamtofastq -i {input} -fq {output[0]} -fq2 {output[1]} 2> {log}'


rule hostremove_fastq_gzip:
    priority: 2
    output:
        '{sn}/host_removed/R{r}.fastq.gz'
    input:
        '{sn}/host_removed/R{r}.fastq'
    benchmark:
        "{sn}/benchmarks/hostremove_fastq_gzip_R{r}.benchmark.tsv"
    shell:
        'gzip {input}'


###### Based on github.com/connor-lab/ncov2019-artic-nf/blob/master/modules/illumina.nf#L124 ######

rule reference_bwa_build:
    conda: 
        'conda_envs/snp_mapping.yaml'
    output:
        '{sn}/core/reference.bwt'
    input:
        reference = config['viral_reference_genome'],
    log:
        '{sn}/core/bwa-build.log'
    benchmark:
        "{sn}/benchmarks/reference_bwa_build.benchmark.tsv"
    params:
        output_prefix = "{sn}/core/reference"
    shell:
        'bwa index -p {params.output_prefix} {input} >{log} 2>&1'

rule reference_bwa_map:
    threads: 2
    conda: 
        'conda_envs/snp_mapping.yaml'
    output:
        '{sn}/core/reference.bam'
    input:
        r1 = '{sn}/host_removed/R1.fastq.gz',
        r2 = '{sn}/host_removed/R2.fastq.gz',
        ref = '{sn}/core/reference.bwt'
    benchmark:
        "{sn}/benchmarks/reference_bwa_map.benchmark.tsv"
    log:
        '{sn}/core/bwa.log'
    params:
       ref_prefix = '{sn}/core/reference'
    shell:
        '(bwa mem -t {threads} {params.ref_prefix} '
        '{input.r1} {input.r2} | '
        'samtools view -bS | samtools sort -@{threads} -o {output}) 2> {log}'

rule run_bed_primer_trim:
    conda: 
        'conda_envs/ivar.yaml'
    input:
        "{sn}/core/reference.bam"
    output:
        sorted_trimmed_mapped_bam = "{sn}/core/reference.mapped.primertrimmed.sorted.bam",
        trimmed_mapped_bam = "{sn}/core/reference.mapped.primertrimmed.bam",
        mapped_bam = "{sn}/core/reference.mapped.bam"
    benchmark:
        "{sn}/benchmarks/bed_primer_trim.benchmark.tsv"
    log:
        "{sn}/core/ivar_trim.log"
    params:
        scheme_bed = config['scheme_bed'],
        ivar_output_prefix = "{sn}/core/reference.mapped.primertrimmed",
        min_len = config['min_len'],
        min_qual = config['min_qual'],
    shell:
        """
        samtools view -F4 -o {output.mapped_bam} {input}
        samtools index {output.mapped_bam}
        ivar trim -e -i {output.mapped_bam} -b {params.scheme_bed} \
            -m {params.min_len} -q {params.min_qual} \
            -p {params.ivar_output_prefix} 2> {log}
        samtools sort -o {output.sorted_trimmed_mapped_bam} \
            {output.trimmed_mapped_bam} 
        """

rule run_ivar_consensus:
    conda: 
        'conda_envs/ivar.yaml'
    output:
        '{sn}/core/virus.consensus.fa'
    input:
        "{sn}/core/reference.mapped.primertrimmed.sorted.bam"
    log:
        '{sn}/core/ivar_consensus.log'
    benchmark:
        "{sn}/benchmarks/ivar_consensus.benchmark.tsv"
    params:
        mpileup_depth = config['mpileup_depth'],
        ivar_min_coverage_depth = config['ivar_min_coverage_depth'],
        ivar_freq_threshold = config['ivar_freq_threshold'],
        output_prefix = '{sn}/core/virus.consensus'
    shell:
        '(samtools mpileup -A -d {params.mpileup_depth} -Q0 {input} | '
        'ivar consensus -t {params.ivar_freq_threshold} '
        '-m {params.ivar_min_coverage_depth} -n N -p {params.output_prefix}) '
        '2>{log}'

rule run_ivar_variants:
    conda: 
        'conda_envs/ivar.yaml'
    output:
        '{sn}/core/ivar_variants.tsv'
    input:
        reference = config['viral_reference_genome'],
        read_bam = "{sn}/core/reference.mapped.primertrimmed.sorted.bam"
    log:
        '{sn}/core/ivar_variants.log'
    benchmark:
        "{sn}/benchmarks/ivar_variants.benchmark.tsv"
    params:
        output_prefix = '{sn}/core/ivar_variants',
        ivar_min_coverage_depth = config['ivar_min_coverage_depth'],
        ivar_min_freq_threshold = config['ivar_min_freq_threshold'],
        ivar_min_variant_quality = config['ivar_min_variant_quality']
    shell:
        '(samtools mpileup -A -d 0 --reference {input.reference} -B '
            '-Q 0 {input.read_bam} | '
        'ivar variants -r {input.reference} -m {params.ivar_min_coverage_depth} '
        '-p {params.output_prefix} -q {params.ivar_min_variant_quality} '
        '-t {params.ivar_min_freq_threshold}) 2> {log}'


################################   Based on scripts/breseq.sh   ####################################


rule run_breseq:
    threads: 4
    priority: 1
    conda: 'conda_envs/snp_mapping.yaml'
    output:
        '{sn}/breseq/output/index.html'
    input:
        expand('{{sn}}/host_removed/R{r}.fastq.gz', r=[1,2])
    log:
        '{sn}/breseq/breseq.log',
    benchmark:
        "{sn}/benchmarks/run_breseq.benchmark.tsv"
    params:
        ref = config['breseq_reference'],
	outdir = '{sn}/breseq'
    shell:
        'breseq --reference {params.ref} --num-processors {threads} --polymorphism-prediction --brief-html-output --output {params.outdir} {input} >{log} 2>&1'


##################  Based on scripts/hisat2.sh and scripts/coverage_stats_avg.sh  ##################


rule coverage_bwa_build:
    conda: 'conda_envs/snp_mapping.yaml'
    output:
        '{sn}/coverage/genome.bwt'
    input:
        '{sn}/core/virus.consensus.fa'
    log:
        '{sn}/coverage/bwa-build.log'
    benchmark:
        "{sn}/benchmarks/coverage_bwa_build.benchmark.tsv"
    params:
        genome = '{sn}/coverage/genome'
    shell:
        'bwa index -p {params.genome} {input} >{log} 2>&1'


rule coverage_bwa:
    threads: 2
    conda: 'conda_envs/snp_mapping.yaml'
    output:
        '{sn}/coverage/output.bam'
    input:
        r1 = '{sn}/host_removed/R1.fastq.gz',
        r2 = '{sn}/host_removed/R2.fastq.gz',
        ref = '{sn}/coverage/genome.bwt'
    benchmark:
        "{sn}/benchmarks/coverage_bwa.benchmark.tsv"
    log:
        '{sn}/coverage/bwa.log'
    params:
        genome = '{sn}/coverage/genome'
    shell:
        'bwa mem -t {threads} {params.genome} '
        '{input.r1} {input.r2} 2> {log} | '
        'samtools view -bS | samtools sort -@{threads} -o {output}'


rule coverage_depth:
    conda: 'conda_envs/snp_mapping.yaml'
    output:
        '{sn}/coverage/depth.txt'
    input:
        '{sn}/coverage/output.bam'
    benchmark:
        "{sn}/benchmarks/coverage_depth.benchmark.tsv"
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
    benchmark:
        "{sn}/benchmarks/run_kraken2.benchmark.tsv"
    params:
        outdir = '{sn}/kraken2',
	    db = os.path.abspath(config['kraken2_db'])
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
        '{sn}/core/virus.consensus.fa'
    benchmark:
        "{sn}/benchmarks/lmat_pretile.benchmark.tsv"
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
    benchmark:
        "{sn}/benchmarks/run_lmat_{db}.benchmark.tsv"
    log:
        "{sn}/lmat/lmat_{db}.log"
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
	    ' --threads={threads} 2> {log}'


rule lmat_postprocess:
    output:
        '{sn}/lmat/parseLMAT_output.txt'
    input:
        expand('{{sn}}/lmat/assembly.tiled.fasta.{db}.lo.rl_output0.out', db=[config['lmat_db']])
    benchmark:
        "{sn}/benchmarks/lmat_postprocess.benchmark.tsv"
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
         '{sn}/core/virus.consensus.fa'
    log:
         '{sn}/quast/quast.log'
    benchmark:
        "{sn}/benchmarks/run_quast.benchmark.tsv"
    params:
         outdir = '{sn}/quast',
         genome = config['viral_reference_genome'],
         fcoords = config['viral_reference_feature_coords']
    shell:
         'quast {input} -r {params.genome} -g {params.fcoords} --output-dir {params.outdir} --threads {threads} >{log}'

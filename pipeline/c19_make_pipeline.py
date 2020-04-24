#!/usr/bin/env python3

"""
c19_make_pipeline.py: creates a snakemake-based pipeline directory, given a set of input .fastq.gz files.
The input .fastq.gz filenames must contain substrings such as '_R1_' or '_R2_' indicating the read direction.

Usage:

    # First create the directory $HOME/data, which should contain a copy of (or be a symlink to)
    # galaxylab:/home/kmsmith/data (warning: 61 GB!)

    # Create pipeline
    ./c19_make_pipeline.py -o Iran1 $HOME/data/MT-swab-Iran-Liverpool*.fastq.gz

    # Run pipeline (cacheing conda envs in $HOME/.snakemake)
    cd Iran1/   # directory created by 'c19_make_pipeline.py'
    snakemake -p --cores=16 --use-conda --conda-prefix=$HOME/.snakemake all

TODO: this script only creates single-sample pipelines, but the Snakemake workflow supports multiple samples.
"""

import os
import re
import sys
import shutil
import argparse
from collections import OrderedDict


class Pipeline:
    """
    The default constructor initializes a Pipeline from command-line arguments in sys.argv.
    The following members are initialized:

       self.outdir                    pipeline directory to be created
       self.original_fastq_files_R1   list of "out-of-tree" filenames specified on command line
       self.original_fastq_files_R2   list of "out-of-tree" filenames specified on command line
       self.input_fastq_files_R1      list of "in-tree" filenames, relative to sample subdir (not toplevel dir)
       self.input_fastq_files_R2      list of "in-tree" filenames, relative to sample subdir (not toplevel dir)
       self.copies                    for debugging: number of redundant copies of sample to be analyzed in parallel
    """

    # TODO datadir currently hardcoded as '$HOME/data', should introduce command-line flag to override default
    datadir = os.path.join(os.environ['HOME'], 'data')
    
    # TODO lots of hardcoded parameters here, is this what we want?
    
    # Used as -a,-A arguments to 'cutadapt'
    primer_rc = os.path.join(datadir, 'wuhan_primers_28.01.20_trim_RC.fa')
    primer_fw = os.path.join(datadir, 'wuhan_primers_28.01.20_trim_FW.fa')

    # Last arguments on 'trimmomatic' command line (after input, output files)
    trimmomatic_args = 'ILLUMINACLIP:/home/kmsmith/data/NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:4:20'

    # Used as hisat2 reference genome when removing host sequences.
    # Also used as 'ivar' reference genome in variant detection + consensus.
    hostremove_reference = os.path.join(datadir, 'MN908947_3.fasta')

    # Used as --reference argument to 'breseq'
    breseq_reference = os.path.join(datadir, 'MN908947_primer_annotated_prot_clinical.gb')

    # Used as --db argument to 'kraken2'
    kraken2_db = os.path.join(datadir, 'Kraken2/db')

    # Consensus and variant calling ivar/samtools params from https://github.com/connor-lab/ncov2019-artic-nf/blob/master/conf/illumina.config
    mpileup_depth = 100000
    # ivar frequency threshold to build consensus
    ivar_freq_threshold = 0.75
    # Minimum coverage depth to call variant 
    ivar_min_coverage_depth = 10
    # iVar frequency threshold to call variant (ivar variants: -t )
    ivar_min_freq_threshold = 0.25
    # iVar minimum mapQ to call variant (ivar variants: -q)
    ivar_min_variant_quality = 20

    # lmat_fragment_size: size of fragments (in bp) analyzed by 'lmat'
    # Absolute pathname of the LMAT DB is {lmat_basedir}/data/{lmat_db}.
    # LMAT's expected "runtime inputs" (e.g. 'ncbi_taxid_to_rank.txt') should be in {lmat_basedir}/runtime_inputs.
    lmat_fragment_size = 250
    lmat_basedir = os.path.join(datadir, 'LMAT-1.2.6')
    lmat_db = 'kML+Human.v4-14.20.g10.db'

    # Used as -r,-g arguments to 'quast'
    quast_reference_genome = os.path.join(datadir, 'MN908947_3.fasta')
    quast_feature_coords = os.path.join(datadir, 'MN908947_3.gff3')
    
    
    def __init__(self):
        parser = argparse.ArgumentParser()
        parser.add_argument('-o', '--outdir', required=True, help='pipeline directory (will be created, must not already exist)')
        parser.add_argument('--copies', type=int, default=1, help='for debugging: number of redundant copies of sample to be analyzed in parallel')
        parser.add_argument('input_fastq_files', nargs='+', help='list of .fastq.gz input files')

        args = parser.parse_args()
        assert args.copies >= 1

        # Currently 'prefix' is deduced from the output directory name.  (This behavior could
        # be overridden with optional command-line argument, if this seems like a useful feature.)
        
        self.outdir = args.outdir
        self.ncopies = args.copies
        
        # We fail with an error if the pipeline output directory already exists.
        # This behavior is heavy-handed, but ensures that a pipeline user can never
        # inadvertantly overwrite or modify previous runs.
        # TODO: add command-line flag to override this behavior?

        if os.path.exists(self.outdir):
            self._die(f"output directory {self.outdir} already exists; must be deleted or renamed before making a new pipeline.")

        # The rest of this method mostly consists of sanity-checking code, to verify that
        # the input filenames "pair" nicely into (R1,R2) pairs.  If this sanity check fails,
        # then we currently fail with an error -- is this the right thing to do?
        #
        # An example of a well-formed (R1,R2) filename pair:
        #   MT-swab-Iran-Liverpool-pool1_S3_L001_R1_001.fastq.gz
        #   MT-swab-Iran-Liverpool-pool1_S3_L001_R2_001.fastq.gz
            
        # Hash (prefix, suffix) -> (R1_filename, R2_filename)
        pair_analysis = OrderedDict()
        
        for input_file in args.input_fastq_files:
            # Currently we require all input files to be .fastq.gz files.
            # TODO: allow input files to be either .fastq or .fastq.gz instead?
            if not input_file.endswith('.fastq.gz'):
                self._die(f"expected input filename {input_file} to end in .fastq.gz")
                
            if not os.path.exists(input_file):
                self._die(f"input file {input_file} not found")

            # Regex pattern for identifying a substring such as '_R1_' or '_R2_' indicating read direction.
            # We allow a little flexibility here (case-insensitive, underscore can be replaced by hyphen or period).
            # This pattern is just a guess and we may want to rethink it later!
            pattern = r'([-_.][rR][12][-_.])'
            
            b = os.path.basename(input_file)
            
            m = re.search(pattern, b)
            if m is None:
                self._die(f"expected input filename {input_file} to contain substring such as '_R1_' or '_R2_' indicating read direction")
                
            assert m.group(0)[2] in ['1','2']
            r = int(m.group(0)[2])    # read direction (either 1 or 2)
            x = b[:m.start()]         # prefix (part of basename preceding regex pattern match)
            y = b[m.end():]           # suffix (part of basename following regex pattern match)
                  
            if re.search(pattern, y) is not None:
                self._die(f"input filename {input_file} contains multiple substrings such as '_R1_' or '_R2_' indicating read direction")
                  
            if (x,y) not in pair_analysis:
                pair_analysis[x,y] = [None, None]
            if pair_analysis[x,y][r-1] is not None:
                self._die(f"confused by similar filenames {pair_analysis[x,y][r-1]} and {input_file}")
                  
            pair_analysis[x,y][r-1] = input_file

        # Reorganize the input filenames into (R1,R2) pairs.
        
        self.original_fastq_files_R1 = [ ]
        self.original_fastq_files_R2 = [ ]
        pair_analysis_ok = True
        
        for ((x,y),(f1,f2)) in pair_analysis.items():
            self.original_fastq_files_R1.append(f1)
            self.original_fastq_files_R2.append(f2)

            if (f1 is None) or (f2 is None):
                pair_analysis_ok = False

        # If the input filenames don't "pair" nicely into (R1,R2) pairs,
        # we currently fail with an error and let the user sort it out.
        
        if not pair_analysis_ok:
            print("Fatal: couldn't pair input filenames into (R1,R2) pairs", file=sys.stderr)
            print("Filename pair analysis follows", file=sys.stderr)

            for (f1,f2) in self.original_fastq_file_pairs:
                for f in (f1,f2):
                    t = f if (f is not None) else "[** missing **]"
                    print(f"  {t}", file=sys.stderr, end='')
                    
            print(file=sys.stderr)
            sys.exit(1)

        # Assign in-tree filenames to the input files.  (Each input file is represented by
        # an "original" filename which was specified on the command line, and an "in-tree"
        # or "input" copy which will be used in the actual pipeline.  The original files
        # are copied to their in-tree counterparts in self.copy_input_fastq_files().)
            
        self.input_fastq_files_R1 = [ ]
        self.input_fastq_files_R2 = [ ]
        
        for f in self.original_fastq_files_R1:
            self.input_fastq_files_R1.append(os.path.join('fastq_inputs', os.path.basename(f)))

        for f in self.original_fastq_files_R2:
            self.input_fastq_files_R2.append(os.path.join('fastq_inputs', os.path.basename(f)))

    
    def write(self):
        self.write_config_yaml()
        self.copy_workflow_files()
        self.copy_input_fastq_files()

    
    def write_config_yaml(self):
        """Writes {pipeline_output_dir}/config.yaml."""
        
        filename = os.path.join(self.outdir, 'config.yaml')
        self._mkdir_for_file(filename)
        
        print(f"Writing {filename}")

        with open(filename,'w') as f:
            print(f"# Autogenerated by c19_make_pipeline.py -- do not edit!", file=f)
            print(file=f)
            
            print(f"# This file contains a high-level summary of pipeline configuration and inputs.", file=f)
            print(f"# It is ingested by the Snakefile, and also intended to be human-readable.", file=f)
            print(file=f)
            
            print("# Used as -a,-A arguments to 'cutadapt'", file=f)
            print(f"primer_rc: {repr(self.primer_rc)}", file=f)
            print(f"primer_fw: {repr(self.primer_fw)}", file=f)
            print(file=f)
            
            print(f"# Last arguments on 'trimmomatic' command line (after input, output files)", file=f)
            print(f"trimmomatic_args: {repr(self.trimmomatic_args)}", file=f)
            print(file=f)

            print("# Used as hisat2 reference genome when removing host sequences.", file =f)
            print("# Also used as 'ivar' reference genome in variant detection + consensus.", file =f)
            print(f"hostremove_reference: {repr(self.hostremove_reference)}", file=f)
            print(file=f)
            
            print(f"# Used as --reference argument to 'breseq'", file=f)
            print(f"breseq_reference: {repr(self.breseq_reference)}", file=f)
            print(file=f)
            
            print(f"# Used as --db argument to 'kraken2'", file=f)
            print(f"kraken2_db: {repr(self.kraken2_db)}", file=f)
            print(file=f)

            print(f"# Consensus and variant calling ivar/samtools params from https://github.com/connor-lab/ncov2019-artic-nf/blob/master/conf/illumina.config", file=f)
            print(f"mpileup_depth: {repr(self.mpileup_depth)}", file=f)
            print(f"# ivar frequency threshold to build consensus", file=f)
            print(f"ivar_freq_threshold: {repr(self.ivar_freq_threshold)}", file=f)
            print(f"# Minimum coverage depth to call variant", file=f)
            print(f"ivar_min_coverage_depth: {repr(self.ivar_min_coverage_depth)}", file=f)
            print(f"# iVar frequency threshold to call variant (ivar variants: -t )", file=f)
            print(f"ivar_min_freq_threshold: {repr(self.ivar_min_freq_threshold)}", file=f)
            print(f"# iVar minimum mapQ to call variant (ivar variants: -q)", file=f)
            print(f"ivar_min_variant_quality: {repr(self.ivar_min_variant_quality)}", file=f)
            print(file=f)
            
            print(f"# lmat_fragment_size: size of fragments (in bp) analyzed by 'lmat'", file=f)
            print(f"# Absolute pathname of the LMAT DB is {{lmat_basedir}}/data/{{lmat_db}}.", file=f)
            print(f"# LMAT's expected \"runtime inputs\" (e.g. 'ncbi_taxid_to_rank.txt') should be in {{lmat_basedir}}/runtime_inputs.", file=f)
            print(f"lmat_fragment_size: {repr(self.lmat_fragment_size)}", file=f)
            print(f"lmat_basedir: {repr(self.lmat_basedir)}", file=f)
            print(f"lmat_db: {repr(self.lmat_db)}", file=f)
            print(file=f)

            print(f"# Used as -r,-g arguments to 'quast'", file=f)
            print(f"quast_reference_genome: {repr(self.quast_reference_genome)}", file=f)
            print(f"quast_feature_coords: {repr(self.quast_feature_coords)}", file=f)
            print(file=f)

            print(f"# List of samples to be analyzed follows.", file=f)
            print(f"# Dictionary is keyed by sample name (=directory name).", file=f)
            print(f"# Input filenames are relative to sample directory (or absolute).", file=f)

            if self.ncopies > 1:
                print(f"# In this artificial example, we analyze the same sample {self.ncopies} times in parallel.", file=f)
            
            print(file=f)            
            print(f"samples:", file=f)

            for s in range(1, self.ncopies+1):
                print(f"  'sample{s}':", file=f)
                
                for r in [1, 2]:
                    print(f"    input_fastq_files_R{r}:", file=f)
                    
                    for filename in getattr(self, f"input_fastq_files_R{r}"):
                        print(f"      - {filename}", file=f)

            
    def copy_workflow_files(self):
        """Writes {pipeline_output_dir}/Snakefile, and {pipeline_output_dir}/conda_envs/*.yaml."""

        # List of (src_relpath, dst_relpath) pairs
        # Source filenames are relative to the directory containing the c19_make_pipeline.py script
        # Destination filenames are relative to the toplevel pipeline directory.
        todo = [ ('Snakefile.master', 'Snakefile'),
                 ('../scripts/fatile', 'fatile'),
                 ('../scripts/parseLMAT', 'parseLMAT') ]

        for conda_envname in [ 'trim_qc', 'assembly', 'assembly_qc', 'snp_mapping', 'ivar', 'lmat' ]:
            filename = f'conda_envs/{conda_envname}.yaml'
            todo.append((filename, filename))
            
        for (src_relpath, dst_relpath) in todo:
            # TODO this 'src_filename' is OK if we're running 'c19_make_pipeline.py' out of the git repository,
            # but will fail if the script is installed anywhere.
            src_filename = os.path.join(os.path.dirname(__file__), src_relpath)
            dst_filename = os.path.join(self.outdir, dst_relpath)
            
            self._copyfile(src_filename, dst_filename)

    
    def copy_input_fastq_files(self):
        """Copies the input .fastq.gz files into the pipeline directory."""

        # The "original" filenames are source files specified on the command line,
        # and the "input" filenames are copies in the pipeline directory.  We copy
        # each original file to the corresponding input.
        
        src_list = self.original_fastq_files_R1 + self.original_fastq_files_R2
        dst_list = self.input_fastq_files_R1 + self.input_fastq_files_R2

        for (src, dst) in zip(src_list, dst_list):
            for i in range(1, self.ncopies+1):
                dst2 = os.path.join(self.outdir, f"sample{i}", dst)
                self._copyfile(src, dst2)

            
    def _die(self, msg):
        """Helper method which prints an error message and exits."""
        
        print(f"Fatal: {msg}", file=sys.stderr)
        sys.exit(1)

    
    def _mkdir_for_file(self, filename):
        dirname = os.path.dirname(filename)
        
        if not os.path.isdir(dirname):
            print(f"Creating directory {dirname}")
            os.makedirs(dirname)

            
    def _copyfile(self, src_filename, dst_filename):
        self._mkdir_for_file(dst_filename)
        
        print(f"Copying {src_filename} -> {dst_filename}")
        shutil.copyfile(src_filename, dst_filename)
        
        

if __name__ == '__main__':
    p = Pipeline()
    p.write()

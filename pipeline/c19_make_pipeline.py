#!/usr/bin/env python3

"""
c19_make_pipeline.py: creates a snakemake-based pipeline directory, given a set of input .fastq.gz files.
The input .fastq.gz filenames must contain substrings such as '_R1_' or '_R2_' indicating the read direction.

Currently the pipeline is a placeholder: just concatenates, sorts, and runs fastqc (for each read direction).
This placeholder design is intended to be a prototype for comments on the overall software design.

Usage (on galaxylab):

    # Create pipeline
    c19_make_pipeline.py -o Iran1 MT-swab-Iran-Liverpool*.fastq.gz

    # Run pipeline
    cd Iran1/   # directory created by 'c19_make_pipeline.py'
    snakemake --cores=16 all

After the pipeline finishes, you'll see the following files:

   Iran1/Snakefile
   Iran1/fastq_inputs/*.fastq.gz                     original fastq files (copied for completeness)
   Iran1/fastq_sorted/Iran1_R[12].fastq.gz           fastq file for each read direction, after concatenate/sort
   Iran1/fastq_sorted/Iran1_R[12].fastqc.[html,zip]  fastqc output files for each read direction

TODO: add a file 'metadata.json' which contains:
  - docker image ID of the pipeline container (when the pipeline is dockerized!)
  - md5 hash of each input .fastq.gz file
"""

import os
import re
import sys
import shutil
import argparse
from collections import OrderedDict


class Pipeline:
    def __init__(self):
        self.parse_args()
        self.sanity_check_args()
        self.assign_filenames()

        
    def write(self):
        self.create_directory_layout()
        self.write_snakefile()
        self.copy_input_fastq_files()

        
    def parse_args(self):
        """
        Parses command line and initializes the following members.

        self.outdir              pipeline directory to be created
        self.input_fastq_files   list of .fastqc.gz filenames
        self.prefix              prefix prepended to many output filenames (e.g. fastq_sorted/{prefix}_R1.fastq.gz)
        """
        
        parser = argparse.ArgumentParser(description='This text will appear at the top of parser.print_help()',
                                         epilog='This text will appear at the bottom of parser.print_help()')

        parser.add_argument('-o', '--outdir', required=True, help='pipeline directory (will be created, must not already exist)')
        parser.add_argument('input_fastq_files', nargs='+', help='list of .fastq.gz input files')

        args = parser.parse_args()

        self.outdir = args.outdir
        self.input_fastq_files = args.input_fastq_files

        # Currently the output filename prefix is deduced from the output directory.
        # (Could be overridden with optional command-line argument, if this seems like a useful feature.)
        self.prefix = os.path.basename(args.outdir)

        
    def sanity_check_args(self):
        """
        Does some initial sanity checks and initializes the following member.

        self.input_fastq_file_pairs: list of .fastqc.gz filenames, reorganized into list of "paired" (R1,R2) filenames.
        """

        # We fail with an error if the pipeline output directory already exists.
        # This behavior is heavy-handed, but ensures that a pipeline user can never inadvertantly overwrite or modify previous runs.
        # TODO: add command-line flag to override this behavior? I found it inconvenient while debugging the c19_make_pipeline script!
        
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
        
        for input_file in self.input_fastq_files:
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
            i = 0 if (m.group(0)[2] == '1') else 1    # either 0 or 1 depending on read direction
            x = b[:m.start()]                         # prefix (part of basename preceding regex pattern match)
            y = b[m.end():]                           # suffix (part of basename following regex pattern match)
                  
            if re.search(pattern, y) is not None:
                self._die(f"input filename {input_file} contains multiple substrings such as '_R1_' or '_R2_' indicating read direction")
                  
            if (x,y) not in pair_analysis:
                pair_analysis[x,y] = [None, None]
            if pair_analysis[x,y][i] is not None:
                self._die(f"confused by similar filenames {pair_analysis[x,y][i]} and {input_file}")
                  
            pair_analysis[x,y][i] = input_file

        pair_analysis_ok = True
        self.input_fastq_file_pairs = [ ]

        # This loop reorganizes the input filenames into (R1,R2) pairs.
        for ((x,y),(f1,f2)) in pair_analysis.items():
            self.input_fastq_file_pairs.append((f1,f2))
            if (f1 is None) or (f2 is None):
                pair_analysis_ok = False

        # If the input filenames don't "pair" nicely into (R1,R2) pairs,
        # we currently fail with an error and let the user sort it out.
        
        if not pair_analysis_ok:
            print("Fatal: couldn't pair input filenames into (R1,R2) pairs", file=sys.stderr)
            print("Filename pair analysis follows", file=sys.stderr)

            for (f1,f2) in self.input_fastq_file_pairs:
                for f in (f1,f2):
                    t = f if (f is not None) else "[** missing **]"
                    print(f"  {t}", file=sys.stderr, end='')
                    
            print(file=sys.stderr)
            sys.exit(1)

            
    def assign_filenames(self):
        """
        This method (called by the constructor) generates filenames which will be used
        throughout the pipeline, but does not actually create the files.

        In this placeholder pipeline, the following members are initialized:

           self.unsorted_fastq_files    # list of (R1,R2) pairs
           self.sorted_fastq_files      # single (R1,R2) pair
        """

        self.unsorted_fastq_files = []

        for (src1,src2) in self.input_fastq_file_pairs:
            dst1 = os.path.join('fastq_inputs', os.path.basename(src1))
            dst2 = os.path.join('fastq_inputs', os.path.basename(src2))
            self.unsorted_fastq_files.append((dst1,dst2))
        
        self.sorted_fastq_files = [
            os.path.join('fastq_sorted', f"{self.prefix}_R1.fastq.gz"),
            os.path.join('fastq_sorted', f"{self.prefix}_R2.fastq.gz")
        ]


    def create_directory_layout(self):
        """Creates pipeline output directory and a few subdirectories."""
        
        self._mkdir(self.outdir)

        for subdir in ['fastq_inputs', 'fastq_sorted']:
            self._mkdir(os.path.join(self.outdir, subdir))

            
    def write_snakefile(self):
        """Writes {pipeline_output_dir}/Snakefile."""
        
        filename = os.path.join(self.outdir, "Snakefile")
        print(f"Writing {filename}")

        with open(filename, 'w') as f:
            self._write_snakefile_rule_toplevel(f)
            
            for ix in [0,1]:
                self._write_snakefile_rule_concat_and_sort(f, ix, fastqc=True)

    
    def copy_input_fastq_files(self):
        """Copies the input .fastq.gz files into the pipeline directory."""
        
        for ((src1,src2),(dst1,dst2)) in zip(self.input_fastq_file_pairs, self.unsorted_fastq_files):
            self._copyfile(src1, dst1)
            self._copyfile(src2, dst2)

            
    def _die(self, msg):
        """Helper method which prints an error message and exits."""
        
        print(f"Fatal: {msg}", file=sys.stderr)
        sys.exit(1)

    
    def _mkdir(self, dirname):
        """Like os.makedirs(), but noisier."""
        
        print(f"Creating directory {dirname}")
        os.makedirs(dirname)


    def _copyfile(self, src, dst):
        """Like shutil.copyfile(), but noisier, and assumes 'dst' is relative to the pipeline output directory."""
        
        dst = os.path.join(self.outdir, dst)
        print(f"Copying {src} -> {dst}")
        shutil.copyfile(src, dst)


    def _fastqc_outfile(self, fastq_filename):
        """Returns FASTQC html output file, given a .fastq.gz filename."""

        assert fastq_filename.endswith('.fastq.gz')
        return fastq_filename[:-9] + '_fastqc.html'


    def _write_snakefile_rule_toplevel(self, f):
        """
        Writes the Snakefile rule for 'snakemake all'.
        The 'f' argument is the Snakefile file object.
        """

        # In this placeholder pipeline, the only targets are the FASTQC output files from the sorted R1/R2 fastq data.
        targets = [ self._fastqc_outfile(f) for f in self.sorted_fastq_files ]
        
        self._write_snakefile_rule(f, 'all', targets, None, None)
        
        
    def _write_snakefile_rule_concat_and_sort(self, f, ix, fastqc=True):
        """
        Writes the Snakefile rule which concatenates and sorts the input fastq's corresponding
        to a single read direction.

        The 'f' argument is the Snakefile file object.
        The 'ix' argument is the read direction (0=R1, 1=R2).
        """

        name = f"concat_and_sort_R{ix+1}"
        input_files = [ t[ix] for t in self.unsorted_fastq_files ]
        output_file = self.sorted_fastq_files[ix]
        shell_command = f'zcat {" ".join(input_files)} | paste - - - - | sort -k1,1 -t " " | tr "\\t" "\\n" | gzip > {output_file}'
        self._write_snakefile_rule(f, name, input_files, output_file, shell_command)

        if fastqc:
            name = f"fastqc_sorted_R{ix+1}"
            self._write_snakefile_rule_fastqc(f, name, output_file)
        

    def _write_snakefile_rule_fastqc(self, f, name, fastq_filename):
        """
        Writes the Snakefile rule for a single FASTQC run.
        The 'f' argument is the Snakefile file object.
        """
        
        shell_command = f'fastqc {fastq_filename}'
        outfile = self._fastqc_outfile(fastq_filename)
        self._write_snakefile_rule(f, name, [fastq_filename], outfile, shell_command)

    
    def _write_snakefile_rule(self, f, name, input_files, output_file, shell_command):
        """
        Helper method which writes a single Snakefile rule.
        The 'f' argument is the Snakefile file object.        
        """
        
        assert len(input_files) > 0
        
        print(f"{name}:", file=f)
        print(f"    input:", file=f)

        for (i,x) in enumerate(input_files):
            y = ',' if (i < len(input_files)-1) else ''
            print(f"        {repr(x)}{y}", file=f)

        if output_file is not None:
            print(f"    output:", file=f)
            print(f"        {repr(output_file)}", file=f)

        if shell_command is not None:
            print(f"    shell:", file=f)
            print(f"        {repr(shell_command)}", file=f)
        
        print(file=f)



if __name__ == '__main__':
    p = Pipeline()
    p.write()

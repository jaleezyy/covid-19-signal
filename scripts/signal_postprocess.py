#!/usr/bin/env python3

import os
import re
import sys
import glob
import json
import zipfile
import contextlib
import html.parser

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

long_git_id = '$Id: b158164f87c79271ddc9d1083e64e4be1fc26d8e $'

assert long_git_id.startswith('$Id: ')
#short_git_id = long_git_id[5:12]
short_git_id = "v1.5.6"

# Suppresses matplotlib warning (https://github.com/jaleezyy/covid-19-signal/issues/59)
# Creates a small memory leak, but it's nontrivial to fix, and won't be a practical concern!
plt.rcParams.update({'figure.max_open_warning': 0})
plt.style.use('seaborn-whitegrid')


########################    Helper functions/classes for text file parsing   #######################


def file_is_missing(filename, allow_missing=True):
    """
    This helper function is called in several places where we want
    to detect a missing file, and either print a warning or throw
    an exception, depending on whether 'allow_missing' is True.
    """

    if os.path.exists(filename):
        return False
    if not allow_missing:
        raise RuntimeError(f"File {filename} does not exist")
    print(f"Warning: file {filename} does not exist")
    return True


def read_file(filename, allow_missing=True, zname=None):
    """
    Generator which yields all lines in specified file.

    If allow_missing=True, then missing file is treated as empty file
    (but a warning is printed).

    If 'zname' is specified, then 'filename' should be a .zip file,
    and 'zname' should be the name of the constituent file to be read.
    """

    if file_is_missing(filename, allow_missing):
        pass
    elif zname is None:
        with open(filename) as f:
            for line in f:
                yield line
    else:
        with zipfile.ZipFile(filename) as z:
            with z.open(zname) as f:
                for line in f:
                    yield line.decode('ascii')

class QUASTParser(html.parser.HTMLParser):

    def __init__(self):
        html.parser.HTMLParser.__init__(self)
        self.extracting = False
        self.quast_data = ""


    def handle_starttag(self, tag, attr):
        if tag == 'div' and attr == [('id', 'total-report-json')]:
            self.extracting = True
        else:
            return

    def handle_data(self, data):
        if self.extracting == True:
            self.quast_data += data

    def handle_endtag(self, tag):
        if tag == 'div' and self.extracting == True:
            self.extracting = False

    def convert_data_to_json(self):
        # convert to json and then simplify as each report has only one sample
        data = json.loads(self.quast_data.strip())

        data = data['report']

        # simplify report by creating metric: metric_value pairs
        simplified_report = {}
        for report_group in data:
            if report_group[1] != []:
                for metric in report_group[1]:
                    metric_name = metric['metricName'].strip()
                    metric_value = metric['values'][0]

                    if metric_name in simplified_report:
                        print(f"{metric_name} collision in report")
                    else:
                        simplified_report[metric_name] = metric_value
        return simplified_report


class TextFileParser:
    """
    This helper class is used in several places to parse text files for
    multiple 'fields', each field specified by a regular expression.
    """

    def __init__(self):
        self._field_names = [ ]
        self._field_details = [ ]  # List of 4-tuples (regexp_pattern, regexp_group, dtype, required)

    def add_field(self, field_name, regexp_pattern, regexp_group=1, dtype=str, required=True, reducer=None):
        assert field_name not in self._field_names
        self._field_names.append(field_name)
        self._field_details.append((regexp_pattern, regexp_group, dtype, required, reducer))

    def parse_file(self, filename, allow_missing=True, zname=None):
        """
        Parses the specified file and returns a dict (field_name) -> (parsed_value).

        If allow_missing=True, then missing file is treated as empty file
        (but a warning is printed).

        If 'zname' is specified, then 'filename' should be a .zip file,
        and 'zname' should be the name of the constituent file to be read.

        If a regexp fails to match, then either an exception is thrown, or the
        corresponding parsed_value is set to None, depending on whether 'required'
        was True when add_field() was called.
        """

        if file_is_missing(filename, allow_missing):
            return { name: None for name in self._field_names }

        ret = { name: [] for name in self._field_names }

        for line in read_file(filename, allow_missing, zname):
            for (name, (regexp_pattern, regexp_group, dtype, _, _)) in zip(self._field_names, self._field_details):
                m = re.match(regexp_pattern, line)
                if m is not None:
                    val = dtype(m.group(regexp_group))
                    ret[name].append(val)

        for (name, (_,_,_,required,reducer)) in zip(self._field_names, self._field_details):
            if required and len(ret[name]) == 0:
                raise RuntimeError(f"{filename}: failed to parse field '{name}'")
            if reducer is not None:
                ret[name] = reducer(ret[name])
            elif len(ret[name]) > 1:
                raise RuntimeError(f"{filename}: field '{name}' parsed more than once")
            elif len(ret[name]) == 1:
                ret[name] = ret[name][0]
            else:
                ret[name] = None

        return ret


def comma_separated_int(s):
    """
    Used as the 'dtype' argument to TextFileParser.add_field(), to parse integer fields
    which appear in text files with comma separators.
    """

    s = s.replace(',','')
    return int(s)


########################    Helper functions/classes for HTML file parsing   #######################


class SimpleHTMLTableParser(html.parser.HTMLParser):
    """Helper class for parse_html_tables()."""

    def __init__(self):
        html.parser.HTMLParser.__init__(self)

        # List of list of list of strings
        # self.tables[table_index][row_index][col_index] -> string
        self.tables = [ ]

        # State machine
        self.in_table = False
        self.in_tr = False
        self.in_td = False

    def handle_starttag(self, tag, attrs):
        if tag == 'table':
            assert not self.in_table
            self.in_table = True
            self.tables.append([])
        elif tag == 'tr':
            assert self.in_table
            assert not self.in_tr
            self.in_tr = True
            self.tables[-1].append([])
        elif tag in ['td','th']:
            assert self.in_tr
            assert not self.in_td
            self.tables[-1][-1].append('')
            self.in_td = True

    def handle_endtag(self, tag):
        if tag == 'table':
            assert self.in_table
            assert not self.in_tr
            self.in_table = False
        elif tag == 'tr':
            assert self.in_tr
            assert not self.in_td
            self.in_tr = False
        elif tag in ['td','th']:
            assert self.in_td
            self.in_td = False

    def handle_data(self, data):
        if self.in_td:
            s = self.tables[-1][-1][-1]
            sep = ' ' if (len(s) > 0) else ''
            self.tables[-1][-1][-1] = f"{s}{sep}{data}"

def parse_html_tables(html_filename):
    """
    Reads all tables from an HTML file and returns a 3-d array
    indexed by (table-index, row-index, col-index).

    Limited in scope! Assumes non-nested tables, and perfectly-formed
    HTML, e.g. each <tr> and <td> tag must be matched by a </tr> or </td>.
    If anything goes wrong, an exception will be thrown!
    """

    with open(html_filename) as f:
        p = SimpleHTMLTableParser()
        p.feed(f.read())
        return p.tables


def show_html_tables(html_tables):
    """Pretty-prints the return value of parse_html_tables(). Intended for debugging."""

    for (it,t) in enumerate(html_tables):
        print(f"Table {it}")
        for (ir,r) in enumerate(t):
            print(f"  Row {ir}")
            for (ic,c) in enumerate(r):
                print(f"    Col {ic}: {c}")


#########################    Helper functions for handling missing data    #########################


def binop(x, y, op):
    """
    Returns op(x,y), where either x or y can be None, to indicate 'missing data'.
    The op() argument is a binary callable, e.g. op=min or op=(lambda x,y:x+y).
    """

    if (x is None) and (y is None):
        return None

    x = x if (x is not None) else 0
    y = y if (y is not None) else 0
    return op(x,y)


def xround(x, ndigits):
    """Rounds x to the specified number of digits, where x can be None to indicate 'missing data'."""

    return round(x,ndigits) if (x is not None) else None


########################    Parsing functions for pipeline output files   ##########################


def parse_trim_galore_log(filename, allow_missing=True):
    """Returns dict (field_name) -> (parsed_value), see code for list of field_names."""

    t = TextFileParser()
    t.add_field('read_pairs_processed', r'Total reads processed:\s+([0-9,]+)', dtype=comma_separated_int, reducer=min)
    t.add_field('read_pairs_written', r'Reads written \(passing filters\):\s+([0-9,]+)\s+', dtype=comma_separated_int, reducer=min)
    t.add_field('base_pairs_processed', r'Total basepairs processed:\s+([0-9,]+)\s+', dtype=comma_separated_int, reducer=sum)
    t.add_field('base_pairs_written', r'Total written \(filtered\):\s+([0-9,]+)\s+', dtype=comma_separated_int, reducer=sum)

    return t.parse_file(filename, allow_missing)


def parse_fastqc_output(zip_filename, allow_missing=True):
    """Returns dict (field_name) -> (parsed_value), see code for list of field_names."""

    assert zip_filename.endswith('_fastqc.zip')
    zname_data = f"{os.path.basename(zip_filename[:-4])}/fastqc_data.txt"
    zname_summ = f"{os.path.basename(zip_filename[:-4])}/summary.txt"

    t = TextFileParser()
    t.add_field('total_sequences', r'Total Sequences\s+(\d+)', dtype=int)
    t.add_field('flagged_sequences', r'Sequences flagged as poor quality\s+(\d+)', dtype=int)

    ret = t.parse_file(zip_filename, allow_missing, zname_data)
    ret['summary'] = { }   # dict (text -> flavor) pairs, where flavor is in ['PASS','WARN','FAIL']

    for line in read_file(zip_filename, allow_missing, zname_summ):
        (flavor, text, _) = line.split('\t')
        assert flavor in ['PASS','WARN','FAIL']
        assert text not in ret['summary']
        ret['summary'][text] = flavor

    return ret


def parse_fastqc_pair(zip_filename1, zip_filename2, allow_missing=True):
    """Returns dict (field_name) -> (parsed_value), see code for list of field_names."""

    fastqc_r1 = parse_fastqc_output(zip_filename1)
    fastqc_r2 = parse_fastqc_output(zip_filename2)

    seq_tot = binop(fastqc_r1['total_sequences'], fastqc_r2['total_sequences'], lambda x,y:x+y)
    flagged_tot = binop(fastqc_r1['flagged_sequences'], fastqc_r2['flagged_sequences'], lambda x,y:x+y)
    read_pairs = binop(fastqc_r1['total_sequences'], fastqc_r2['total_sequences'], min)

    summary = { }
    for text in list(fastqc_r1['summary'].keys()) + list(fastqc_r2['summary'].keys()):
        flavor1 = fastqc_r1['summary'].get(text, 'PASS')
        flavor2 = fastqc_r2['summary'].get(text, 'PASS')

        summary[text] = 'PASS'
        if 'WARN' in [flavor1,flavor2]:
            summary[text] = 'WARN'
        if 'FAIL' in [flavor1,flavor2]:
            summary[text] = 'FAIL'

    if fastqc_r1['total_sequences'] != fastqc_r2['total_sequences']:
        summary['R1/R2 read count mismatch'] = 'FAIL'

    if (flagged_tot is not None) and (flagged_tot > 0):
        summary[f'{flagged_tot} sequences flagged as poor quality'] = 'WARN'

    return { 'total_sequences': seq_tot,
             'flagged_sequences': flagged_tot,
             'read_pairs': read_pairs,
             'summary': summary }


def parse_kraken2_report(report_filename, allow_missing=True):
    """Returns dict (field_name) -> (parsed_value), see code for list of field_names."""

    t = TextFileParser()
    t.add_field('sars_cov2_percentage', r'\s*([\d\.]*)\s+.*Severe acute respiratory syndrome coronavirus 2', dtype=float)
    try:
        t_dict = t.parse_file(report_filename, allow_missing)
    except RuntimeError:
        t_dict = {'sars_cov2_percentage': 0}
    return t_dict
    #return t.parse_file(report_filename, allow_missing)


def parse_hostremove_hisat2_log(log_filename, allow_missing=True):
    """
    Returns dict (field_name) -> (parsed_value), see code for list of field_names.
    No longer used, but kept around in case it's useful to resurrect it some day.
    """

    t = TextFileParser()
    t.add_field('alignment_rate', r'([\d\.]*)%\s+overall alignment rate', dtype=float)
    return t.parse_file(log_filename, allow_missing)


def parse_quast_report(report_filename, allow_missing=True):
    """Returns dict (field_name) -> (parsed_value), see code for list of field_names."""

    # unfortunately only the quast report html contains all the fields
    # need for summaries, fortunately, it is all encoded in easily extractable
    # json
    q = QUASTParser()
    try:
        with open(report_filename) as fh:
            report = fh.read()
        q.feed(report)
        quast_report = q.convert_data_to_json()
    except FileNotFoundError:
        print("Warning: file %s does not exist" %(report_filename))
        report = None
        quast_report ={'Total length (>= 0 bp)': 0, "# N's per 100 kbp": 0}


    ret = {}
    ret['genome_length'] = float(quast_report['Total length (>= 0 bp)'])
    ret['Ns_per_100_kbp'] = float(quast_report["# N's per 100 kbp"])

    # if genome fails to align to reference these all fail thus the try/except
    try:
        ret['genomic_features'] = str(quast_report['# genomic features'])
        ret['mismatches'] = float(quast_report['# mismatches'])
        ret['mismatches_per_100_kbp'] = float(quast_report['# mismatches per 100 kbp'])
        ret['indels'] = float(quast_report['# indels'])
        ret['indels_per_100_kbp'] = float(quast_report['# indels per 100 kbp'])
        ret['genome_fraction'] = float(quast_report['Genome fraction (%)'])
    except KeyError:
        ret['genomic_features'] = "Failure to align to reference"
        ret['mismatches'] = "Failure to align to reference"
        ret['mismatches_per_100_kbp'] = "Failure to align to reference"
        ret['indels'] = "Failure to align to reference"
        ret['indels_per_100_kbp'] = "Failure to align to reference"
        ret['genome_fraction'] = 0.0

    gfrac = ret['genome_fraction']
    ret['qc_gfrac'] = "PASS" if ((gfrac is not None) and (gfrac >= 90)) else "FAIL"

    # to add a failure mode if the reference fails to align
    indels = ret['indels']
    if indels == 0:
        ret['qc_indel'] = "PASS"
    elif type(indels) == float:
        if indels > 0:
            ret['qc_indel'] = "WARN"
    else:
        ret['qc_indel'] = "FAIL"


    return ret

def parse_consensus_assembly(fasta_filename, allow_missing=True):
    """Returns dict (field_name) -> (parsed_value), see code for list of field_names."""

    if file_is_missing(fasta_filename, allow_missing):
        return { 'N5prime': None, 'N3prime': None }

    lines = open(fasta_filename).readlines()
    if len(lines) != 2:
        assert "".join(lines).count(">") == 1
        lines = [lines[0], "".join(i.strip("\n") for i in lines[1:])]

    line = lines[1].rstrip()
    prime5 = prime3 = 0
    n = len(line)

    # Count leading N's
    for i in range(n):
        if line[i] != 'N':
            break
        prime5 = i+1

    # Count trailing N's
    for i in range(n):
        if line[n-1-i] != 'N':
            break
        prime3 = i+1

    return { 'N5prime': prime5, 'N3prime': prime3 }


def parse_coverage(depth_filename, allow_missing=True):
    """Returns dict (field_name) -> (parsed_value), see code for list of field_names."""

    delims = [ 0, 10, 100, 1000, 2000, 10000]
    nbins = len(delims)+1

    bin_labels = ['0'] + [f"{delims[i-1]+1}x-{delims[i]}x" for i in range(1,nbins-1)] + [f"> {delims[-1]}x"]
    bin_labels = [ f"Fraction with {l} coverage" for l in bin_labels ]

    ret = {
        'bin_labels': bin_labels,
        'bin_fractions': [ None for b in range(nbins) ],
        'mean_coverage': None,
        'qc_meancov': 'FAIL',
        'qc_cov100': 'FAIL',
        'qc_cov1000': 'FAIL',
        'cov100': 0
    }

    if file_is_missing(depth_filename, allow_missing):
        return ret

    coverage = []
    for line in open(depth_filename):
        t = line.split('\t')
        assert len(t) == 3
        coverage.append(int(float(t[2].strip("\n"))))

    coverage = np.array(coverage)
    bin_assignments = np.searchsorted(np.array(delims), coverage, side='left')
    bin_fractions = np.bincount(bin_assignments, minlength=nbins) / float(len(coverage))
    assert bin_fractions.shape == (nbins,)


    ret['cov100'] = np.mean(coverage >= 100)
    ret['bin_fractions'] = [ xround(f,3) for f in bin_fractions ]
    ret['mean_coverage'] = xround(np.mean(coverage), 1)
    ret['qc_meancov'] = "PASS" if (np.mean(coverage) >= 2000) else "FAIL"
    ret['qc_cov100'] = "PASS" if (np.mean(coverage >= 100) >= 0.9) else "FAIL"
    ret['qc_cov1000'] = "PASS" if (np.mean(coverage >= 1000) >= 0.9) else "WARN"

    return ret


def parse_lmat_output(lmat_dirname, allow_missing=True):
    """
    Returns dict (field_name) -> (parsed_value), see code for list of field_names.

    No longer used (LMAT isn't part of pipeline any more), but kept around in case
    it's useful to resurrect it some day.
    """

    # Represent each taxon by a 4-tuple (nreads, score, rank, name)
    taxa = [ ]
    nreads_tot = 0

    for filename in glob.glob(f"{lmat_dirname}/*.fastsummary"):
        for line in open(filename):
            line = line.rstrip('\r\n')

            t = line.split('\t')
            assert len(t) == 4

            i = t[3].rfind(',')
            assert i >= 0

            (score, nreads, ncbi_id) = (float(t[0]), int(t[1]), int(t[2]))
            (rank, name) = (t[3][:i], t[3][(i+1):])
            assert nreads > 0

            taxon = (nreads, score, rank, name)
            taxa.append(taxon)
            nreads_tot += nreads

    if (not allow_missing) and (nreads_tot == 0):
        raise RuntimeError(f"couldn't find fastsummary files in lmat dir '{lmat_dirname}'")

    top_taxa = [ ]
    top_taxa_ann = [ ]
    nreads_cumul = 0

    for (nreads, score, rank, name) in reversed(sorted(taxa)):
        # Roundoff-tolerant version of (nreads_cumul >= 0.9 * nreads_tot)
        if 10*nreads_cumul >= 9*nreads_tot:
            break

        percentage = 100.*nreads/nreads_tot

        top_taxa.append(name)
        top_taxa_ann.append(f"{name} ({rank}, {percentage:.1f}%)")
        nreads_cumul += nreads

    # 'top_taxa_ann' = "top taxa with annotations"
    return { 'top_taxa': top_taxa, 'top_taxa_ann': top_taxa_ann }


def parse_ivar_variants(tsv_filename, allow_missing=True):
    """Returns dict (field_name) -> (parsed_value), see code for list of field_names."""

    if file_is_missing(tsv_filename, allow_missing):
        return { 'variants': [] }

    variants = []

    # Skip first line
    for line in open(tsv_filename).readlines()[1:]:
        t = line.split('\t')
        assert len(t) == 19

        if t[3] != '':
            variants.append(f"{t[2]}{t[1]}{t[3]}")

    return { 'variants': variants }

def parse_freebayes_variants(vcf_filename, allow_missing=True):
    """Returns dict (field_name) -> (parsed_value), see code for list of field_names."""

    if file_is_missing(vcf_filename, allow_missing):
        return { 'variants': [], 'run': False }

    variants = []

    # Only interpret lines that DO NOT start with "#"
    for line in open(vcf_filename):
        if not line.startswith("#"):
            t = line.split('\t')
            assert len(t) == 10

            if t[4] != '':
                variants.append(f"{t[3]}{t[1]}{t[4]}")

    return { 'variants': variants, 'run': True }

def parse_consensus_compare(vcf_filename, allow_missing=True):
    """Returns dict (field_name) -> (parsed_value), see code for list of field_names."""

    if file_is_missing(vcf_filename, allow_missing):
        return { 'positions': [], 'run': False }

    positions = []

    # Only interpret lines that DO NOT start with "#"
    for line in open(vcf_filename):
        if not line.startswith("#"):
            t = line.split('\t')
            assert len(t) == 7

            if t[4] != '':
                positions.append(f"{t[3]}{t[1]}{t[4]}")

    return { 'positions': positions, 'run': True }


def parse_lineage(tsv_filename, sample_names, allow_missing=True):
    """Returns dict (field_name) -> (parsed_value), see code for list of field_names."""

    samples = {}

    if file_is_missing(tsv_filename, allow_missing):
        for name in sample_names:
            samples[name] = { 'lineage' : None,
                              'clade': None,
                              'pangolin_ver': None,
                              'pangodata_ver': None,
                              'nextclade_ver': None }
        return { 'samples': samples }

    lineages = pd.read_table(tsv_filename, sep='\t')
    try:
        df = lineages[['isolate',
                        'pango_lineage',
                        'nextstrain_clade',
                        'pangolin_version',
                        'pangoLEARN_version',
                        'nextclade_version'
                        ]]
    except KeyError:
        df = lineages[['isolate',
                        'pango_lineage',
                        'nextstrain_clade',
                        'pangolin_version',
                        'version',
                        'nextclade_version'
                        ]]

    # Pull each row, identify sid 
    for row in df.itertuples():
        if row.isolate.startswith("Consensus"):
            sid = re.findall("_(.*?)\.", row.isolate)[0]
        else:
            sid = str(row.isolate)

        assert sid in sample_names

    # Pull Pangolin lineage
        lineage = str(row.pango_lineage)
        clade = str(row.nextstrain_clade)
        pangolin = str(row.pangolin_version)
        try:
            pangodata = str(row.pangoLEARN_version)
        except AttributeError:
            pangodata = str(row.version)
        nextclade = str(row.nextclade_version)
        samples[sid] = { 'lineage' : lineage,
                         'clade': clade,
                         'pangolin_ver': pangolin,
                         'pangodata_ver': pangodata,
                         'nextclade_ver': nextclade }

    assert len(samples) == len(sample_names)
    return { 'samples': samples }

def parse_breseq_output(html_filename, allow_missing=True):
    """Returns dict (field_name) -> (parsed_value), see code for list of field_names."""

    if file_is_missing(html_filename, allow_missing):
        return { 'variants': [], 'qc_varfreq': 'MISSING', 'qc_orf_frameshift': 'MISSING', 'run': False}

    tables = parse_html_tables(html_filename)

    assert len(tables) >= 2
    assert tables[1][0] in [ ['Predicted mutation'], ['Predicted mutations'] ]
    assert tables[1][1] == [ 'evidence', 'position', 'mutation', 'freq', 'annotation', 'gene', 'description']

    for t in tables[2:]:
        assert t[0] in [ ['Unassigned missing coverage evidence'], ['Unassigned new junction evidence'] ]

    variants = [ ]
    qc_varfreq = 'PASS'
    qc_orf_frameshift = 'PASS'

    for row in tables[1][2:]:
        assert len(row) == 7
        (evi, pos, mut, freq, ann, gene, desc) = row

        assert freq.endswith('%')
        freq = freq[:-1]

        # TODO: improve breseq html parsing
        # Currently, the "annotation" is unreadable.
        # The "description" is sometimes readable and sometimes not (e.g. it can contain embedded javascript!)

        # Ad hoc improvement of html parsing for 'gene', may need revisiting
        gene = gene.replace('\xa0','')       # remove cosmetic html '&nbsp;'
        gene = gene.replace('\u2011', '-')   # replace unicode underscore by vanilla underscore

        ann = ann .replace('\xa0','')       # remove cosmetic html '&nbsp;'
        ann = ann.replace('\u2011', '-')   # replace unicode underscore by vanilla underscore
        ann = ann.replace(' ', '').replace('(', ' (')

        if ann.count('(') > 1:
            ann = ann.replace(')', ') & ', 1)
            gene = gene.replace('→', '→ &', 1)

        variant = f"{evi};\t{pos};\t{mut};\t({freq}%);\t'{ann}';\t'{gene}'"
        variants.append(variant)

        if float(freq[:-1]) < 90:
            qc_varfreq = 'WARN'

        # check for deletions i.e. delta and insertions i.e. +TTT in coding regions
        if "coding" in ann:
            # +TTTT
            letter_ins_check = re.search(r"\+([A-Za-z]+)", mut)
            if letter_ins_check:
                # check if insertion is %3 in size i.e. not a frameshift and just
                # a warning
                if len(letter_ins_check.group(0)) % 3 == 0:
                    qc_orf_frameshift = 'WARN'
                else:
                    qc_orf_frameshift = 'FAIL'

            # +136 bp
            size_ins_check = re.search(r"\+([0-9]+) bp", mut)
            if size_ins_check:
                # check if insertion is %3 in size i.e. whole codon
                if int(size_ins_check.group(1)) % 3 == 0:
                    qc_orf_frameshift = 'WARN'
                else:
                    qc_orf_frameshift = 'FAIL'

            # Δ22 bp
            del_check = re.search(r"Δ([0-9]+)", mut)
            if del_check:
                if int(del_check.group(1)) % 3 == 0:
                    qc_orf_frameshift = 'WARN'
                else:
                    qc_orf_frameshift = 'FAIL'

    return { 'variants': variants, 'qc_varfreq': qc_varfreq,
            'qc_orf_frameshift': qc_orf_frameshift, 'run': True}


########  Base classes for writing summary files, see WriterBase docstring for explanation  ########


class WriterBase:
    """
    The postprocessing script writes a bunch of output files with similar contents.

    It's convenient to represent this by one "writer" class per output file, with the
    following class hierarchy:

        WriterBase
           SampleTextWriter        writes {sample_name}_sample.txt (single-sample)
           HTMLWriterBase
              SampleHTMLWriter     writes {sample_name}_sample.html (single-sample)
              SummaryHTMLWriter    writes summary.html (multi-sample)

    TODO: add single-sample PDF, multi-sample CSV.

    For each file, the write_sample() argument processes one sample. It will be called
    once for "single-sample" outputs, and multiple times for "multi-sample" outputs.

    Note that write_sample() is implemented in the base class, but is factored into
    multiple methods, in case a subclass wants to override part of the logic (e.g.
    SampleHTMLWriter overrides write_breseq()).
    """


    def __init__(self, filename, unabridged):
        self.filename = filename
        self.unabridged = unabridged
        self.pipeline_name = f'SARS-CoV-2 Illumina GeNome Assembly Line (SIGNAL), version {short_git_id}'
        self.pipeline_url = 'https://github.com/jaleezyy/covid-19-signal'
        self.pipeline_note = "Note: Asterisks (*) indicates a discrepancy between iVar (default) and FreeBayes (if run)"

        self.f = open(filename, 'w')


    def start_sample(self, s):
        """
        The 's' argument should be an instance of type Sample (defined later in this file).

        Each call to start_sample() is followed by:
           - one or more calls to write_lines()
           - one or more key/value blocks, delimited by start_kv_pairs..end_kv_pairs,
             and containing one or more calls to write_kv_pair()
           - one call to end_sample()
        """
        raise RuntimeError('To be overridden by subclass')


    def start_kv_pairs(self, title, link_filenames=[]):
        raise RuntimeError('To be overridden by subclass')


    def write_kv_pair(self, key, val, indent=0, qc=False):
        """
        The 'key' argument is a string, and often contains newlines.  The newlines
        are intended to be helpful for explicit formatting in narrow columns, but
        some subclasses ignore them via "val = val.replace('\n','').

        The 'val' argument can be None, to indicate missing data.

        The 'qc' boolean flag can be used to enable special formatting for QC flags
        (e.g. color-coding in HTML tables).
        """
        raise RuntimeError('To be overridden by subclass')


    def end_kv_pairs(self):
        raise RuntimeError('To be overridden by subclass')


    def write_lines(self, title, lines, coalesce=False):
        """
        The 'lines' argument is a list of strings.
        The 'coalesce' argument tells the subclass that coalescing multiple lines is okay.
        """
        raise RuntimeError('To be overridden by subclass')


    def end_sample(self, s):
        raise RuntimeError('To be overridden by subclass')


    def close(self):
        if self.f is not None:
            self.f.close()
            self.f = None
            print(f"Wrote {self.filename}")


    def write_data_volume_summary(self, s):
        self.start_kv_pairs("Data Volume", link_filenames=[f"adapter_trimmed/{s.name}_trim_galore.log"])
        self.write_kv_pair("Raw\nData\n(read\npairs)", s.trim_galore['read_pairs_processed'], indent=1)

        if self.unabridged:
            self.write_kv_pair("Raw Data (base pairs)", s.trim_galore['base_pairs_processed'], indent=1)
            self.write_kv_pair("Post Primer Removal (read pairs)", s.trim_galore['read_pairs_written'], indent=1)
            self.write_kv_pair("Post Primer Removal (base pairs)", s.trim_galore['base_pairs_written'], indent=1)

        self.write_kv_pair("Post\nTrim\n(read\npairs)", s.post_trim_qc['read_pairs'], indent=1)
        # self.write_kv_pair("Post\nhuman\npurge\n(%)", s.hostremove['alignment_rate'], indent=1)
        self.end_kv_pairs()


    def write_qc_flags(self, s):
        self.start_kv_pairs("Quality Control Flags")

        key = "Genome Fraction greater than 90%" if self.unabridged else "Genome\nfraction\n>90%"
        self.write_kv_pair(key, s.quast['qc_gfrac'], indent=1, qc=True)

        key = "No indels detected (maximum length 85bp)" if self.unabridged else "No\nindels"
        self.write_kv_pair(key, s.quast['qc_indel'], indent=1, qc=True)

        key = "Depth of coverage >= 2000x" if self.unabridged else "Depth\n>2000"
        self.write_kv_pair(key, s.coverage['qc_meancov'], indent=1, qc=True)

        key = "All variants with at least 90% frequency among reads" if self.unabridged else "Variants\n>90%"
        if 'MISSING' not in s.breseq['qc_varfreq']:
            self.write_kv_pair(key, s.breseq['qc_varfreq'], indent=1, qc=True)

        key = "Frameshifts in SARS-CoV-2 open reading frames" if self.unabridged else "ORF\nFrameshifts"
        if 'MISSING' not in s.breseq['qc_orf_frameshift']:
            self.write_kv_pair(key, s.breseq['qc_orf_frameshift'], indent=1, qc=True)

        key = "Reads per base sequence quality" if self.unabridged else "Fastqc\nquality"
        val = s.post_trim_qc['summary'].get('Per base sequence quality', 'FAIL')
        self.write_kv_pair(key, val, indent=1, qc=True)

        key = "Sequencing adapter removed" if self.unabridged else "Fastqc\nadapter"
        val = s.post_trim_qc['summary'].get('Adapter Content', 'FAIL')
        self.write_kv_pair(key, val, indent=1, qc=True)

        key = "At least 90% of positions have coverage >= 100x" if self.unabridged else "90%\ncov\n>100"
        self.write_kv_pair(key, s.coverage['qc_cov100'], indent=1, qc=True)

        key = "At least 90% of positions have coverage >= 1000x" if self.unabridged else "90%\ncov\n>1000"
        self.write_kv_pair(key, s.coverage['qc_cov1000'], indent=1, qc=True)

        self.end_kv_pairs()


    def write_fastqc_summary(self, s):
        if not self.unabridged:
            return

        self.start_kv_pairs("FASTQC Flags", link_filenames=[f"adapter_trimmed/{s.name}_R{r}_val_{r}_fastqc.html" for r in [1,2]])

        for flavor in [ 'FAIL', 'WARN' ]:
            for (msg,f) in s.post_trim_qc['summary'].items():
                if msg in [ 'Sequence Duplication Levels', 'Overrepresented sequences' ]:
                    # From Natalie Knox (https://github.com/jaleezyy/covid-19-signal/issues/54):
                    # Fastqc HTML report should not include "sequence duplication" or "overrepresented sequence"
                    # flags as these will generally fail with an amplicon-based protocol.
                    continue
                if f == flavor:
                    self.write_kv_pair(msg, f, indent=1, qc=True)

        self.end_kv_pairs()


    def write_kraken2(self, s):
        self.start_kv_pairs("Kraken2", link_filenames=[f"kraken2/{s.name}_kraken2.report"])
        self.write_kv_pair("Reads\nSARS-CoV-2\n(%)", s.kraken2['sars_cov2_percentage'], indent=1)
        self.end_kv_pairs()


    def write_quast(self, s):
        self.start_kv_pairs("QUAST", link_filenames=[f"quast/{s.name}_quast_report.html"])
        self.write_kv_pair("Genome\nLength\n(bp)", s.quast['genome_length'], indent=1)
        self.write_kv_pair("Genome\nFraction\n(%)", s.quast['genome_fraction'], indent=1)
        self.write_kv_pair("N's per\n100 kbp", s.quast['Ns_per_100_kbp'], indent=1)

        if self.unabridged:
            self.write_kv_pair("Genomic Features", s.quast['genomic_features'], indent=1)
            self.write_kv_pair("Mismatches", s.quast['mismatches'], indent=1)
            self.write_kv_pair("Mismatches per 100 kbp", s.quast['mismatches_per_100_kbp'], indent=1)
            self.write_kv_pair("Indels", s.quast['indels'], indent=1)
            self.write_kv_pair("Indels per 100 kbp", s.quast['indels_per_100_kbp'], indent=1)

        self.write_kv_pair("Average\nDepth of\nCoverage", s.coverage['mean_coverage'], indent=1)

        if self.unabridged:
            for (l,f) in zip(s.coverage['bin_labels'], s.coverage['bin_fractions']):
                self.write_kv_pair(l, f, indent=2)
            self.write_kv_pair("5' Ns", s.consensus['N5prime'], indent=1)
            self.write_kv_pair("3' Ns", s.consensus['N3prime'], indent=1)

        self.end_kv_pairs()


    def write_ivar(self, s):
        title = "Variants in Consensus Genome (iVar)" if self.unabridged else "Variants (iVar)"
        self.write_lines(title, s.ivar['variants'], coalesce=True)

    def write_breseq(self, s):
        if s.breseq['run'] == True:
            title = "Variants in Read Alignment (BreSeq)" if self.unabridged else "Variants (BreSeq)"
            self.write_lines(title, s.breseq['variants'])
        else:
            return None

    def write_freebayes(self, s):
        if s.freebayes['run'] == True:
            title = "Unique Variants in Consensus Genome (FreeBayes)" if self.unabridged else "Unique Variants (FreeBayes)"
            self.write_lines(title, s.freebayes['variants'], coalesce=True)
        else:
            return None

    def write_lineage(self, s):
        title = "Pangolin Lineage Assignment" if self.unabridged else "Lineage (Pangolin)"
        self.write_lines(title, [s.lineage['lineage']], coalesce=True)

    def write_clade(self, s):
        title = "Nextclade Clade Assignment" if self.unabridged else "Clade (Nextstrain)"
        self.write_lines(title, [s.lineage['clade']], coalesce=True)


    def write_compare(self, s):
        if s.compare['run'] == True:
            title = "Nucleotide Differences in Consensus Genomes (FreeBayes as reference)" if self.unabridged else "Consensus Nucleotide Differences (FreeBayes as Reference)"
            self.write_lines(title, s.compare['positions'], coalesce=True)
        else:
            return None

    def write_sample(self, s):
        self.start_sample(s)
        self.write_lineage(s)
        self.write_clade(s)
        self.write_data_volume_summary(s)
        self.write_qc_flags(s)
        self.write_fastqc_summary(s)
        self.write_kraken2(s)
        self.write_quast(s)
        self.write_ivar(s)
        self.write_freebayes(s)
        self.write_compare(s)
        self.write_breseq(s)
        self.end_sample(s)


    @staticmethod
    def coalesce_lines(lines, n):
        """Helper method which might be useful in subclass implementation of write_lines()."""

        ret = [ ]

        for line in lines:
            if (len(ret) > 0) and (len(ret[-1])+len(line) < n):
                ret[-1] = f"{ret[-1]} {line}"
            else:
                ret.append(line)

        return ret


class HTMLWriterBase(WriterBase):
    def __init__(self, filename, unabridged):
        WriterBase.__init__(self, filename, unabridged)

        print("<!DOCTYPE html>", file=self.f)
        print("<html>", file=self.f)
        print("<head></head>", file=self.f)
        print("<body>", file=self.f)

        print(f'<h3>{self.pipeline_name}&nbsp;&nbsp;(<a href="{self.pipeline_url}">{self.pipeline_url}</a>)</h3>', file=self.f)
        print(f'<p>{self.pipeline_note}</p>', file=self.f)

    def close(self):
        if self.f is not None:
            print("</body>", file=self.f)
            print("</html>", file=self.f)

        WriterBase.close(self)


####################################################################################################


class SampleTextWriter(WriterBase):
    """Writes single-sample output file {sample_dir}/sample.txt"""

    def __init__(self, filename):
        WriterBase.__init__(self, filename, unabridged=True)

        print(self.pipeline_name, file=self.f)
        print(self.pipeline_url, file=self.f)
        print("", file=self.f)
        print(self.pipeline_note, file=self.f)
        print("", file=self.f)

    def start_sample(self, s):
        print(f"Sample: {s.name}", file=self.f)
        print("", file=self.f)

    def start_kv_pairs(self, title, link_filenames=[]):
        print(title, file=self.f)

    def write_kv_pair(self, key, val=None, indent=0, qc=False):
        s = ' ' * (4*indent)
        k = key.replace('\n',' ')
        v = val if (val is not None) else ''
        t = f'{v}  {k}' if qc else f'{k}: {v}'
        print(f"{s}{t}", file=self.f)

    def write_lines(self, title, lines, coalesce=False):
        if coalesce:
            lines = self.coalesce_lines(lines, 80)

        self.start_kv_pairs(title)
        for line in lines:
            print(f"    {line}", file=self.f)
        self.end_kv_pairs()

    def end_kv_pairs(self):
        print(file=self.f)

    def end_sample(self, s):
        pass


class SampleHTMLWriter(HTMLWriterBase):
    """Writes single-sample output file {sample_dir}/sample.html"""

    def __init__(self, filename):
        HTMLWriterBase.__init__(self, filename, unabridged=True)

    def start_sample(self, s):
        self.sample_name = s.name
        print(f'<h3>Sample: {s.name}</h3>', file=self.f)

        # Start outer table for left/right alignment of summary stats, coverage plot
        print('<p><table><tr>', file=self.f)

        # Start inner table for summary stats
        print('<td style="vertical-align: top"><table>', file=self.f)

    def start_kv_pairs(self, title, link_filenames=[]):
        t = [ ]

        for link_relpath in link_filenames:
            link_abspath = os.path.join(self.sample_name, link_relpath)
            if os.path.exists(link_abspath):
                t.append(f'<a href="{link_relpath}">{os.path.basename(link_relpath)}</a>')
            else:
                t.append(f'{link_relpath} not found')

        if len(t) > 0:
            t = ' | '.join(t)
            title = f'{title} [ {t} ]'

        print(f'<tr><th colspan="2" style="text-align: left">{title}</th>', file=self.f)

    def write_kv_pair(self, key, val=None, indent=0, qc=False):
        k = key.replace('\n',' ')
        v = val if (val is not None) else ''
        td1 = f'<td style="padding-left: {20*indent}px">{k}</td>'
        td2 = f'<td>{v}</td>'
        print(f"<tr>{td1}{td2}</tr>", file=self.f)

    def write_lines(self, title, lines, coalesce=False):
        if coalesce:
            lines = self.coalesce_lines(lines, 80)

        self.start_kv_pairs(title)
        for line in lines:
            print(f'<tr><td colspan="2" style="padding-left: 20px">{line}</td></tr>', file=self.f)
        self.end_kv_pairs()

    def end_kv_pairs(self):
        pass

    def write_breseq(self, s):
        # Override write_breseq(), in favor of including breseq/index.html as an iframe (see below)
        pass

    def end_sample(self, s):
        # End inner table for summary stats
        print('</table></td>', file=self.f)

        # Coverage plot
        print(f'<td style="vertical-align: top"><img src="coverage/{s.name}_coverage_plot.png"></td>', file=self.f)

        # Breseq iframe
        if os.path.exists("breseq/{s.name}_output/index.html"):
            print(f'<iframe src="breseq/{s.name}_output/index.html" width="100%" height="800px" style="border: 0px"></iframe>', file=self.f)


class SummaryHTMLWriter(HTMLWriterBase):
    """Writes multi-sample output file {sample_dir}/summary.html"""

    def __init__(self, filename, maxlines=8):
        HTMLWriterBase.__init__(self, filename, unabridged=False)

        self.first_row = ''
        self.second_row = ''
        self.current_row = ''
        self.num_rows_written = 0
        self.maxlines = maxlines

        self.current_group_text = None
        self.current_group_colspan = 0

        self.css_lborder = 'border-left: 1px solid black;'
        self.css_bborder = 'border-bottom: 1px solid black;'

        # Single-row table containing two summary plots
        print('<p><table><tr>', file=self.f)
        print('<td><img src="summary_ncov2_in_reads_v_genome_fraction.png"></td>', file=self.f)
        print('<td><img src="summary_average_depth_v_genome_fraction.png"></td>', file=self.f)
        print('</tr>', file=self.f)
        print('<td><img src="summary_highly_covered_v_genome_fraction.png"></td>', file=self.f)
        print('</tr></table>', file=self.f)

        # Start long table containing statistics
        print('<p><table style="border-collapse: collapse; border: 1px solid black;">', file=self.f)


    def css_color(self, val=None, qc=False):
        qc_false = ('#dddddd', '#ffffff')
        qc_true = { 'PASS': ('#7ca37c', '#8fbc8f'),
                    'WARN': ('#ddba00', '#ffd700'),
                    'FAIL': ('#dd2a2a', '#ff3030'),
                    'MISSING': ('#808080', '#ffffff'),
                    'PASS*': ('#7ca37c', '#8fbc8f'),
                    'WARN*': ('#ddba00', '#ffd700'),
                    'FAIL*': ('#dd2a2a', '#ff3030') }

        odd = (self.num_rows_written % 2)
        color = qc_true[val][odd] if qc else qc_false[odd]
        return f"background-color: {color};"


    def start_sample(self, s, link_filenames=[]):
        assert self.current_group_text is None

        link_text = s.name

        if os.path.exists(s.name):
            url = f"{os.path.basename(s.name)}/{os.path.basename(s.name)}_sample.html"
            link_text = f'<a href="{url}">{link_text}</a>'

        self.first_row += f'<th>Sample</th>\n'
        self.second_row += f'<th style="{self.css_bborder}"></th>\n'
        self.current_row += f'<td style="{self.css_color()}">&nbsp;<br>{link_text}<br>&nbsp;</td>'


    def start_kv_pairs(self, title, link_filenames=[]):
        assert self.current_group_text is None
        self.current_group_text = title
        self.current_group_colspan = 0


    def write_kv_pair(self, key, val=None, indent=0, qc=False):
        assert self.current_group_text is not None

        css_c = self.css_color(val, qc)
        css_lb = self.css_lborder if (self.current_group_colspan == 0) else ''

        k = key.replace('\n','<br>')
        v = val if (val is not None) else ''
        s = '&nbsp;' if (self.current_group_colspan > 0) else ''

        self.second_row += f'<td style="{self.css_bborder} {css_lb}">{k}</td>\n'
        self.current_row += f'<td style="{css_c} {css_lb}">{s}{v}</td>\n'
        self.current_group_colspan += 1


    def write_lines(self, title, lines, coalesce=False):
        if len(lines) > self.maxlines:
            n = len(lines)
            m = self.maxlines
            last = f'&nbsp;&nbsp;&nbsp; (+ {n-m+1} more)'
            lines = lines[:(m-1)] + [last]

        val = '\n<p style="margin-bottom:0px; margin-top:8px">'.join(lines)
        val = f'<p style="margin-bottom:0px; margin-top:0px"> {val}'

        self.start_kv_pairs(title)
        self.write_kv_pair('', val)
        self.end_kv_pairs()


    def end_kv_pairs(self):
        assert self.current_group_text is not None

        if self.current_group_colspan > 0:
            self.first_row += f'<th colspan="{self.current_group_colspan}" style="{self.css_lborder}">{self.current_group_text}</th>\n'

        self.current_group_text = None
        self.current_group_colspan = 0


    def end_sample(self, s):
        assert self.current_group_text is None

        if self.num_rows_written == 0:
            print(f'<tr>{self.first_row}</tr>', file=self.f)
            print(f'<tr>{self.second_row}</tr>', file=self.f)

        print(f'<tr>{self.current_row}</tr>', file=self.f)
        self.first_row = ''
        self.second_row = ''
        self.current_row = ''
        self.num_rows_written += 1


    def close(self):
        print('</table>', file=self.f)
        HTMLWriterBase.close(self)


####################################################################################################


class Archive:
    """
    Wrapper class around zipfile.ZipFile, used in Pipeline.write_archive().
    Must be run from toplevel pipeline directory.
    """

    def __init__(self, filename, debug=False):
        self.zipfile = zipfile.ZipFile(filename, 'w')
        self.filename = filename
        self.contents = set()
        self.debug = debug


    def add_file(self, filename):
        assert not os.path.isabs(filename)

        if not os.path.exists(filename):
            raise RuntimeError(f'File {filename} not found')

        if filename in self.contents:
            return

        if self.debug:
            print(f'{self.filename}: adding {filename}')

        self.zipfile.write(filename)
        self.contents.add(filename)


    def add_glob(self, pattern):
        for filename in glob.glob(pattern):
            self.add_file(filename)


    def add_dir(self, topdir, allow_missing=True):
        assert not os.path.isabs(topdir)

        if allow_missing and not os.path.exists(topdir):
            return

        for dirname, subdirs, filenames in os.walk(topdir):
            for f in filenames:
                self.add_file(os.path.join(dirname,f))


    def close(self):
        self.zipfile.close()
        print(f'Wrote {self.filename}')


#############################   High-level classes (Pipeline, Sample)   ############################


class Sample:
    """Must be constructed from toplevel pipeline directory."""

    def __init__(self, name, ivarlin, fblin):
        self.name = name

        self.trim_galore = parse_trim_galore_log(f"{name}/adapter_trimmed/{name}_trim_galore.log")
        self.post_trim_qc = parse_fastqc_pair(f"{name}/adapter_trimmed/{name}_R1_val_1_fastqc.zip", f"{name}/adapter_trimmed/{name}_R2_val_2_fastqc.zip")
        self.kraken2 = parse_kraken2_report(f"{name}/kraken2/{name}_kraken2.report")
        self.quast = parse_quast_report(f"{name}/quast/{name}_quast_report.html")
        self.quast_freebayes = parse_quast_report(f"{name}/freebayes/quast/{name}_quast_report.html")
        self.consensus = parse_consensus_assembly(f"{name}/core/{name}.consensus.fa")
        self.consensus_freebayes = parse_consensus_assembly(f"{name}/freebayes/{name}.consensus.fasta")
        self.coverage = parse_coverage(f"{name}/coverage/{name}_depth.txt")
        self.ivar = parse_ivar_variants(f"{name}/core/{name}_ivar_variants.tsv")
        self.freebayes = parse_freebayes_variants(f"{name}/freebayes/{name}.variants.norm.vcf")
        self.compare = parse_consensus_compare(f"{name}/freebayes/{name}_consensus_compare.vcf")
        self.breseq = parse_breseq_output(f"{name}/breseq/output/index.html")


        if ivarlin['lineage'] != fblin['lineage'] and fblin['lineage'] is not None:
            assert ivarlin['pangolin_ver'] == fblin['pangolin_ver']
            assert ivarlin['pangodata_ver'] == fblin['pangodata_ver']
            if ivarlin['clade'] == fblin['clade']:
                 self.lineage = { 'lineage': str(ivarlin['lineage'] + " (FB: %s)" %(fblin['lineage'])),
                                 'pangolin_ver': ivarlin['pangolin_ver'],
                                 'pangodata_ver': ivarlin['pangodata_ver'],
                                 'clade': ivarlin['clade'] }
            else:
                 self.lineage = { 'lineage': str(ivarlin['lineage'] + " (FB: %s)" %(fblin['lineage'])),
                                 'pangolin_ver': ivarlin['pangolin_ver'],
                                 'pangodata_ver': ivarlin['pangodata_ver'],
                                 'clade': str(ivarlin['clade'] + " (FB: %s)" %(fblin['clade'])) }
        else:
            self.lineage = ivarlin

    # Compare sample consensus and variant outputs if both iVar and FreeBayes are present
    # Only will report iVar but values that differ from FreeBayes will be flagged with an asterisks

    # QC metrics comparison
        param = ["qc_gfrac", "qc_indel"]
        status = ["FAIL", "WARN", "PASS"]
        for item in param:
        # If values don't match, compare whether FreeBayes shows improvement (i.e., no indels) by indexing status
            if self.quast[item] != self.quast_freebayes[item]:
                ivar = status.index(self.quast[item])
                fb = status.index(self.quast_freebayes[item])
                if fb > ivar:
                    self.quast[item] = str(self.quast[item]) + "*" # Ex. WARN*

    # Variant call comparison
        if len(self.freebayes['variants']) > 0:
            variants = []
            # Check if iVar variant found in FreeBayes: remove from FreeBayes if found, tag if not
            # Final output should be a list of unique FreeBayes variants and a list of appropriately tagged iVar variants
            for var in self.ivar['variants']:
                if var in self.freebayes['variants']:
                    variants.append(var)
                    self.freebayes['variants'].remove(var)
                else:
                    variants.append(var + "*")
            self.ivar = { 'variants': variants }
            #if len(self.freebayes['variants']) == 0: self.freebayes['variants'].append("None")

    # Compare consensus assembly (N counts at 5' and 3')
    # Run quick_align to highlight specific differences across consensus genomes
        for itemvar, itemfb in zip(self.consensus, self.consensus_freebayes):
            assert itemvar == itemfb
        # Check if # of N's is fewer than in the FreeBayes consensus (assumed less ambiguious)
            if (self.consensus_freebayes[itemfb] is not None) and (float(self.consensus_freebayes[itemfb]) < float(self.consensus[itemvar])):
                self.consensus[itemvar] = str(self.consensus[itemvar]) + "*"


class Pipeline:
    """Must be constructed from toplevel pipeline directory."""

    def __init__(self, sample_csv_filename):
        sample_csv = pd.read_csv(sample_csv_filename)
        sample_names = sorted(sample_csv['sample'].drop_duplicates().values)

        self.iv_lineage = parse_lineage(f"lineage_assignments.tsv", sample_names)
        self.fb_lineage = parse_lineage(f"freebayes_lineage_assignments.tsv", sample_names)

        self.samples = [ Sample(s, self.iv_lineage['samples'][s], self.fb_lineage['samples'][s]) for s in sample_names ]

        if len(self.samples) == 0:
            raise RuntimeError(f"{sample_csv_filename} contains zero samples, nothing to do!")

    def write_summary_plot1(self):
        """Writes toplevel summary plot: %SARS versus completeness."""

        plt.figure(figsize=(7.6,5.7))

        xvec = [ ]
        yvec = [ ]

        for s in self.samples:
            x = s.kraken2['sars_cov2_percentage']
            y = s.quast['genome_fraction']

            if (x is None) or (y is None):
                continue

            xvec.append(x)
            yvec.append(y)

            if (x < 90.) or (y < 90.):
                depth = s.coverage['mean_coverage']
                label = f'{s.name}'
                if depth is not None:
                    label += f', {int(depth)}x'
                plt.annotate(label, (x,y), xytext=(x,y-5))

        plt.scatter(xvec, yvec, marker='o', facecolor='blue')
        plt.xlabel(r'SARS-COV-2 in Trimmed FASTQ (%)')
        plt.ylabel(r'Genome Fraction (%)')
        plt.xlim(0, 100)
        plt.ylim(0, 100)

        print('Writing summary_ncov2_in_reads_v_genome_fraction.png')
        plt.savefig('summary_ncov2_in_reads_v_genome_fraction.png')
        plt.clf()


    def write_summary_plot2(self):
        """Writes toplevel summary plot: depth versus completeness."""

        plt.figure(figsize=(7.6,5.7))

        xvec = [ ]
        yvec = [ ]

        for s in self.samples:
            x = s.coverage['mean_coverage']
            y = s.quast['genome_fraction']

            if (x is None) or (y is None):
                continue

            xvec.append(x)
            yvec.append(y)

            if True:
                scov2 = s.kraken2['sars_cov2_percentage']
                label = f'{s.name}'
                if scov2 is not None:
                    label += f', {scov2:.1f}%'
                plt.annotate(label, (x,y), xytext=(x,y-5))

        plt.scatter(xvec, yvec, marker='o', facecolor='red')
        plt.xlabel(r'Average Depth of Coverage')
        plt.ylabel(r'Genome Fraction (%)')
        plt.ylim(0, 100)

        print('Writing summary_average_depth_v_genome_fraction.png')
        plt.savefig('summary_average_depth_v_genome_fraction.png')
        plt.clf()

    def write_summary_plot3(self):
        """Writes toplevel summary plot: highly covered % versus completeness."""

        plt.figure(figsize=(7.6,5.7))

        xvec = [ ]
        yvec = [ ]

        for s in self.samples:
            x = s.coverage['cov100']
            y = s.quast['genome_fraction']

            if (x is None) or (y is None):
                continue

            xvec.append(x)
            yvec.append(y)

            if True:
                scov2 = s.kraken2['sars_cov2_percentage']
                label = f'{s.name}'
                if scov2 is not None:
                    label += f', {scov2:.1f}%'
                plt.annotate(label, (x,y), xytext=(x,y-5))

        plt.scatter(xvec, yvec, marker='o', facecolor='red')
        plt.xlabel(r'Fraction of Genome with >100x Coverage')
        plt.ylabel(r'Genome Fraction (%)')
        plt.xlim(0, 1)
        plt.ylim(0, 100)

        print('Writing summary_highly_covered_v_genome_fraction.png')
        plt.savefig('summary_highly_covered_v_genome_fraction.png')
        plt.clf()

    def write_reports(self):
        if len(self.samples) > 1:
            summary_writer = SummaryHTMLWriter('summary.html')
        else:
            summary_writer = SampleHTMLWriter('summary.html')

        for s in self.samples:
            if os.path.exists(s.name):
                w1 = SampleTextWriter(f'{s.name}/{s.name}_sample.txt')
                w2 = SampleHTMLWriter(f'{s.name}/{s.name}_sample.html')
                sample_writers = [ w1, w2 ]
            else:
                print(f"Warning: sample directory {s.name} does not exist")
                sample_writers = [ ]

            for w in [summary_writer] + sample_writers:
                w.write_sample(s)

            for w in sample_writers:
                w.close()

        summary_writer.close()


    def write_archive(self, debug=False):
        if not os.path.exists('summary.html'):
            print(f"write_archive: 'summary.html' does not exist, nothing to do")
            return

        a = Archive('summary.zip', debug)
        a.add_file('summary.html')
        a.add_file('summary_ncov2_in_reads_v_genome_fraction.png')
        a.add_file('summary_average_depth_v_genome_fraction.png')
        a.add_file('summary_highly_covered_v_genome_fraction.png')


        for sample in self.samples:
            s = sample.name
            a.add_glob(f'{s}/{s}_sample.txt')
            a.add_glob(f'{s}/{s}_sample.html')
            a.add_glob(f'{s}/coverage/{s}_coverage_plot.png')
            a.add_glob(f'{s}/adapter_trimmed/{s}_trim_galore.log')
            a.add_glob(f'{s}/adapter_trimmed/{s}_*_fastqc.html')
            a.add_glob(f'{s}/kraken2/{s}_kraken2.report')
            a.add_glob(f'{s}/host_removal/{s}_human_read_mapping.log')
            a.add_glob(f'{s}/quast/*.html')
            a.add_dir(f'{s}/quast/icarus_viewers')
            a.add_glob(f'{s}/breseq/{s}_breseq.log')
            a.add_dir(f'{s}/breseq/{s}_output')

        a.close()


####################################################################################################


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(f"Usage: c19_postprocess.py <sample_table.csv>", file=sys.stderr)
        print(f"Note: must be run from toplevel pipeline directory", file=sys.stderr)
        sys.exit(1)

    p = Pipeline(sys.argv[1])
    p.write_summary_plot1()
    p.write_summary_plot2()
    p.write_summary_plot3()
    p.write_reports()
    p.write_archive()

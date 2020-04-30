#!/usr/bin/env python3

import os
import re
import glob
import zipfile
import html.parser

import numpy as np


####################################################################################################


def file_is_missing(filename, allow_missing=True):
    if os.path.exists(filename):
        return False
    if not allow_missing:
        raise RuntimeError(f"File {filename} does not exist")
    print(f"Warning: file {filename} does not exist")
    return True


def read_file(filename, allow_missing=True, zname=None):
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

                    
class TextFileParser:
    def __init__(self):
        self._column_names = [ ]
        self._column_details = [ ]  # List of 4-tuples (regexp_pattern, regexp_group, dtype, required)
    
    def add_column(self, name, regexp_pattern, regexp_group=1, dtype=str, required=True):
        assert name not in self._column_names
        self._column_names.append(name)
        self._column_details.append((regexp_pattern, regexp_group, dtype, required))
    
    def parse_file(self, filename, allow_missing=True, zname=None):
        ret = { name: None for name in self._column_names }

        if file_is_missing(filename, allow_missing):
            return ret
        
        for line in read_file(filename, allow_missing, zname):
            for (name, (regexp_pattern, regexp_group, dtype, _)) in zip(self._column_names, self._column_details):
                m = re.match(regexp_pattern, line)
                if m is None:
                    continue
                if ret[name] is not None:
                    raise RuntimeError(f"{filename}: attempt to set field '{name}' twice")                    
                ret[name] = dtype(m.group(regexp_group))
        
        for (name, (_,_,_,required)) in zip(self._column_names, self._column_details):
            if required and ret[name] is None:
                raise RuntimeError(f"{filename}: failed to parse column '{name}'")
        
        return ret

    
def comma_separated_int(s):
    s = s.replace(',','')
    return int(s)


####################################################################################################


class SimpleHTMLTableParser(html.parser.HTMLParser):
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
    with open(html_filename) as f:
        p = SimpleHTMLTableParser()
        p.feed(f.read())
        return p.tables
    
def show_html_tables(html_tables):
    """Intended for debugging."""
        
    for (it,t) in enumerate(html_tables):
        print(f"Table {it}")
        for (ir,r) in enumerate(t):
            print(f"  Row {ir}")
            for (ic,c) in enumerate(r):
                print(f"    Col {ic}: {c}")


####################################################################################################


def binop(x, y, op):
    """Returns op(x,y), with None treated as zero."""

    if (x is None) and (y is None):
        return None

    x = x if (x is not None) else 0
    y = y if (y is not None) else 0
    return op(x,y)


def xround(x, ndigits):
    return round(x,ndigits) if (x is not None) else None


####################################################################################################


def parse_cutadapt_log(filename, allow_missing=True):
    t = TextFileParser()
    t.add_column('read_pairs_processed', r'Total read pairs processed:\s+([0-9,]+)', dtype=comma_separated_int)
    t.add_column('R1_with_adapter', r'\s*Read 1 with adapter:\s+([0-9,]+)\s+', dtype=comma_separated_int)
    t.add_column('R2_with_adapter', r'\s*Read 2 with adapter:\s+([0-9,]+)\s+', dtype=comma_separated_int)
    t.add_column('read_pairs_written', r'Pairs written \(passing filters\):\s+([0-9,]+)\s+', dtype=comma_separated_int)
    t.add_column('base_pairs_processed', r'Total basepairs processed:\s+([0-9,]+)\s+', dtype=comma_separated_int)
    t.add_column('base_pairs_written', r'Total written \(filtered\):\s+([0-9,]+)\s+', dtype=comma_separated_int)

    return t.parse_file(filename, allow_missing)


def parse_fastqc_output(zip_filename, allow_missing=True):
    assert zip_filename.endswith('_fastqc.zip')
    zname_data = f"{os.path.basename(zip_filename[:-4])}/fastqc_data.txt"
    zname_summ = f"{os.path.basename(zip_filename[:-4])}/summary.txt"
    
    t = TextFileParser()
    t.add_column('total_sequences', r'Total Sequences\s+(\d+)', dtype=int)
    t.add_column('flagged_sequences', r'Sequences flagged as poor quality\s+(\d+)', dtype=int)
    
    ret = t.parse_file(zip_filename, allow_missing, zname_data)
    ret['summary'] = { }   # dict (text -> flavor) pairs, where flavor is in ['PASS','WARN','FAIL']
    
    for line in read_file(zip_filename, allow_missing, zname_summ):
        (flavor, text, _) = line.split('\t')
        assert flavor in ['PASS','WARN','FAIL']
        assert text not in ret['summary']
        ret['summary'][text] = flavor

    return ret


def parse_fastqc_pair(zip_filename1, zip_filename2, allow_missing=True):
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
        
    return { 'total_sequences': seq_tot, 'flagged_sequences': flagged_tot, 'read_pairs': read_pairs, 'summary': summary }


def parse_kraken2_report(report_filename, allow_missing=True):
    t = TextFileParser()
    t.add_column('sars_cov2_percentage', r'\s*([\d\.]*)\s+.*Severe acute respiratory syndrome coronavirus 2', dtype=float)
    return t.parse_file(report_filename, allow_missing)


def parse_hostremove_hisat2_log(log_filename, allow_missing=True):
    t = TextFileParser()
    t.add_column('alignment_rate', r'([\d\.]*)%\s+overall alignment rate', dtype=float)
    return t.parse_file(log_filename, allow_missing)


def parse_quast_report(report_filename, allow_missing=True):
    t = TextFileParser()

    # Note that some fields are "required=False" here.
    # Sometimes, QUAST doesn't write every field, but I didn't investigate why.
    t.add_column("genome_length", r'Total length \(>= 0 bp\)\s+(\S+)', dtype=int)
    t.add_column("genome_fraction", r'Genome fraction \(%\)\s+(\S+)', dtype=float, required=False)   # Note: genome "fraction" is really a percentage
    t.add_column("genomic_features", r'# genomic features\s+(\S+)', required=False)
    t.add_column("Ns_per_100_kbp", r"# N's per 100 kbp\s+(\S+)", dtype=float)
    t.add_column("mismatches_per_100_kbp", r"# mismatches per 100 kbp\s+(\S+)", dtype=float, required=False)
    t.add_column("indels_per_100_kbp", r"# indels per 100 kbp\s+(\S+)", dtype=float, required=False)
    
    ret = t.parse_file(report_filename, allow_missing=True)
    
    gfrac = ret['genome_fraction']
    ret['qc_gfrac'] = "PASS" if ((gfrac is not None) and (gfrac >= 90)) else "FAIL"
    
    return ret


def parse_consensus_assembly(fasta_filename, allow_missing=True):
    if file_is_missing(fasta_filename, allow_missing):
        return { 'N5prime': None, 'N3prime': None }
    
    lines = open(fasta_filename).readlines()
    assert len(lines) == 2
    
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
        'qc_cov1000': 'FAIL'
    }
    
    if file_is_missing(depth_filename, allow_missing):
        return ret

    coverage = []
    for line in open(depth_filename):
        t = line.split('\t')
        assert len(t) == 3
        coverage.append(int(t[2]))

    coverage = np.array(coverage)
    bin_assignments = np.searchsorted(np.array(delims), coverage, side='left')
    bin_fractions = np.bincount(bin_assignments, minlength=nbins) / float(len(coverage))
    assert bin_fractions.shape == (nbins,)
    
    ret['bin_fractions'] = [ xround(f,3) for f in bin_fractions ]
    ret['mean_coverage'] = xround(np.mean(coverage), 1)
    ret['qc_meancov'] = "PASS" if (np.mean(coverage) >= 2000) else "FAIL"
    ret['qc_cov100'] = "PASS" if (np.mean(coverage >= 100) >= 0.9) else "FAIL"
    ret['qc_cov1000'] = "PASS" if (np.mean(coverage >= 1000) >= 0.9) else "WARN"
    
    return ret


def parse_lmat_output(lmat_dirname, allow_missing=True):
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
    nreads_cumul = 0

    for (nreads, score, rank, name) in reversed(sorted(taxa)):
        # Roundoff-tolerant version of (nreads_cumul >= 0.9 * nreads_tot)
        if 10*nreads_cumul >= 9*nreads_tot:
            break

        percentage = 100.*nreads/nreads_tot
        top_taxa.append(f"{name} ({rank}, {percentage:.1f}%)")
        nreads_cumul += nreads

    return { 'taxonomic_composition': top_taxa }


def parse_ivar_variants(tsv_filename, allow_missing=True):
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


def parse_breseq_output(html_filename, allow_missing=True):
    if file_is_missing(html_filename, allow_missing):
        return { 'variants': [], 'qc_varfreq': 'FAIL' }
    
    tables = parse_html_tables(html_filename)
    
    assert len(tables) >= 2
    assert tables[1][0] in [ ['Predicted mutation'], ['Predicted mutations'] ]
    assert tables[1][1] == [ 'evidence', 'position', 'mutation', 'freq', 'annotation', 'gene', 'description']

    for t in tables[2:]:
        assert t[0] in [ ['Unassigned missing coverage evidence'], ['Unassigned new junction evidence'] ]
    
    variants = [ ]
    qc_varfreq = 'PASS'

    for row in tables[1][2:]:
        assert len(row) == 7
        (evi, pos, mut, freq, ann, gene, desc) = row

        assert freq.endswith('%')
        freq = freq[:-1]
        
        variant = f"{evi} {pos} {mut} {desc} ({freq}% of reads) {ann}"
        variants.append(variant)
        
        if float(freq) < 90:
            qc_varfreq = 'WARN'
    
    return { 'variants': variants, 'qc_varfreq': qc_varfreq }


####################################################################################################


class WriterBase:
    def __init__(self, filename):
        self.filename = filename
        self.f = open(filename, 'w')

    def start_sample(self, sample_name, link_to=None):
        raise RuntimeError('To be overridden by subclass')

    def start_group(self, text):
        raise RuntimeError('To be overridden by subclass')
    
    def write(self, key, val=None, indent=0, optional=False, abbrev=None, qc=False):
        raise RuntimeError('To be overridden by subclass')

    def end_group(self):
        raise RuntimeError('To be overridden by subclass')        
    
    def end_sample(self):
        raise RuntimeError('To be overridden by subclass')
    
    def close(self):
        if self.f is not None:
            self.f.close()
            self.f = None
            print(f"Wrote {self.filename}")


class HTMLWriterBase(WriterBase):
    def __init__(self, filename):
        WriterBase.__init__(self, filename)
        
        print("<!DOCTYPE html>", file=self.f)
        print("<html>", file=self.f)
        print("<head></head>", file=self.f)
        print("<body>", file=self.f)
        
        url = "https://github.com/jaleezyy/covid-19-sequencing"
        print(f'<h3>SARS-CoV-2 Genome Sequencing Consensus & Variants</h3>', file=self.f)
        print(f'<a href="{url}">{url}</a>', file=self.f)

    def close(self):
        if self.f is not None:
            print("</body>", file=self.f)
            print("</html>", file=self.f)
        WriterBase.close(self)


class SampleTextWriter(WriterBase):
    def __init__(self, filename):
        WriterBase.__init__(self, filename)
        
        print("SARS-CoV-2 Genome Sequencing Consensus & Variants", file=self.f)
        print("https://github.com/jaleezyy/covid-19-sequencing", file=self.f)
        print("", file=self.f)
        
    def start_sample(self, sample_name, link_to=None):
        pass

    def start_group(self, text):
        print(text, file=self.f)

    def write(self, key, val=None, indent=0, optional=False, abbrev=None, qc=False):
        s = ' ' * indent
        v = val if (val is not None) else ''
        t = f'{v}  {key}' if qc else f'{key}: {v}'
        print(f"{s}{t}", file=self.f)

    def end_group(self):
        print(file=self.f)

    def end_sample(self):
        self.close()


class SampleHTMLWriter(HTMLWriterBase):
    def __init__(self, filename):
        HTMLWriterBase.__init__(self, filename)

    def start_sample(self, sample_name, link_to=None):
        print("<p><table>", file=self.f)
        
    def start_group(self, text):
        print(f'<tr><th colspan="2" style="text-align: left">{text}</th>', file=self.f)

    def write(self, key, val=None, indent=0, optional=False, abbrev=None, qc=False):
        v = val if (val is not None) else ''
        td1 = f'<td style="padding-left: {5*indent}px">{key}</td>'
        td2 = f'<td>{v}</td>'
        print(f"<tr>{td1}{td2}</tr>", file=self.f)

    def end_group(self):
        pass
    
    def end_sample(self):
        print("</table>", file=self.f)
        self.close()


class SummaryHTMLWriter(HTMLWriterBase):
    def __init__(self, filename):
        HTMLWriterBase.__init__(self, filename)

        self.first_row = ''
        self.second_row = ''
        self.current_row = ''
        self.first_rows_written = False

        self.current_group_text = None
        self.current_group_colspan = 0
        
        self.css_lborder = 'border-left: 1px solid black;'
        self.css_bborder = 'border-bottom: 1px solid black;'
        
        self.css_qc = {
            'PASS': 'background-color:darkseagreen;',
            'WARN': 'background-color: gold;',
            'FAIL': 'background-color: #ff3030;'
        }
        
        print('<p><table style="border-collapse: collapse; border: 1px solid black;">', file=self.f)

        
    def start_sample(self, sample_name, link_to):
        assert isinstance(link_to, SampleHTMLWriter)
        assert self.current_group_text is None

        url = os.path.relpath(link_to.filename, os.path.basename(self.filename))
        
        self.first_row += f'<th>Sample</th>\n'
        self.second_row += f'<th style="{self.css_bborder}"></th>\n'
        self.current_row += f'<td><a href="{url}">{sample_name}</a></td>'
        
        
    def start_group(self, text):
        assert self.current_group_text is None
        self.current_group_text = text
        self.current_group_colspan = 0

        
    def write(self, key, val=None, indent=0, optional=False, abbrev=None, qc=False):
        assert self.current_group_text is not None

        if optional:
            return
        if abbrev is not None:
            key = abbrev

        s1 = self.css_lborder if (self.current_group_colspan == 0) else ''
        s2 = self.css_qc[val] if qc else ''
        v = val if (val is not None) else ''
        
        self.second_row += f'<td style="{self.css_bborder} {s1}">{key}</td>\n'
        self.current_row += f'<td style="{s1} {s2}">{val}</td>\n'
        self.current_group_colspan += 1

        
    def end_group(self):
        assert self.current_group_text is not None
        self.first_row += f'<th colspan="{self.current_group_colspan}" style="{self.css_lborder}">{self.current_group_text}</th>\n'
        self.current_group_text = None
        self.current_group_colspan = 0

        
    def end_sample(self):
        assert self.current_group_text is None

        if not self.first_rows_written:
            print(f'<tr>{self.first_row}</tr>', file=self.f)
            print(f'<tr>{self.second_row}</tr>', file=self.f)
            self.first_rows_written = True

        print(f'<tr>{self.current_row}</tr>', file=self.f)
        self.first_row = ''
        self.second_row = ''
        self.current_row = ''

    def close(self):
        print('</table>', file=self.f)
        HTMLWriterBase.close(self)


####################################################################################################


class Sample:
    def __init__(self, dirname, sample_name):
        self.dirname = dirname
        self.sample_name = sample_name
        
        self.cutadapt = parse_cutadapt_log(f"{dirname}/fastq_primers_removed/cutadapt.log")
        self.post_trim_qc = parse_fastqc_pair(f"{dirname}/fastq_trimmed/R1_paired_fastqc.zip", f"{dirname}/fastq_trimmed/R2_paired_fastqc.zip")
        self.kraken2 = parse_kraken2_report(f"{dirname}/kraken2/report")
        self.hostremove = parse_hostremove_hisat2_log(f"{dirname}/host_removed/hisat2.log")
        self.quast = parse_quast_report(f"{dirname}/quast/report.txt")
        self.consensus = parse_consensus_assembly(f"{dirname}/consensus/virus.consensus.fa")
        self.coverage = parse_coverage(f"{dirname}/coverage/depth.txt")
        self.lmat = parse_lmat_output(f"{dirname}/lmat")
        self.ivar = parse_ivar_variants(f"{dirname}/ivar_variants/ivar_variants.tsv")
        self.breseq = parse_breseq_output(f"{dirname}/breseq/output/index.html")


    def write(self, w, link_to=None):
        w.start_sample(self.sample_name, link_to)
            
        w.start_group("Data Volume")
        w.write("Raw Data (read pairs)", self.cutadapt['read_pairs_processed'], indent=4)
        w.write("Raw Data (base pairs)", self.cutadapt['base_pairs_processed'], indent=4, optional=True)
        w.write("Post Primer Removal (read pairs)", self.cutadapt['read_pairs_written'], indent=4, optional=True)
        w.write("Post Primer Removal (base pairs)", self.cutadapt['base_pairs_written'], indent=4, optional=True)
        w.write("Post Trim (read pairs)", self.post_trim_qc['read_pairs'], indent=4)
        w.write("Post human purge (%)", self.hostremove['alignment_rate'], indent=4)
        w.end_group()
    
        w.start_group("Quality Control Flags")
        w.write("Genome Fraction greater than 90%", self.quast['qc_gfrac'], indent=4, abbrev='gfrac >90%', qc=True)
        w.write("Depth of coverage >= 2000x", self.coverage['qc_meancov'], indent=4, abbrev='depth >2000', qc=True)
        w.write("All variants with at least 90% frequency among reads", self.breseq['qc_varfreq'], indent=4, abbrev='vars >90%', qc=True)
        w.write("Reads per base sequence quality", self.post_trim_qc['summary'].get('Per base sequence quality','FAIL'), indent=4, abbrev='fastqc quality', qc=True)
        w.write("Sequencing adapter removed", self.post_trim_qc['summary'].get('Adapter Content','FAIL'), indent=4, abbrev='fastqc adapter', qc=True)
        w.write("At least 90% of positions have coverage >= 100x", self.coverage['qc_cov100'], indent=4, abbrev='cov >100 >90%', qc=True)
        w.write("At least 90% of positions have coverage >= 1000x", self.coverage['qc_cov1000'], indent=4, abbrev='cov >1000 >90%', qc=True)
        w.end_group()

        w.end_sample()
        return
        
        for flavor in [ 'FAIL', 'WARN' ]:
            messages = [ k for (k,v) in self.post_trim_qc['summary'].items() if v==flavor ]
            if len(messages) == 0:
                continue
            
            w.write("FASTQC {flavor}:")
            for k in messages:
                w.write("    {k}")
            w.write("")

        w.write("Reads SARS-CoV-2 (Kraken2)", self.kraken2['sars_cov2_percentage'], indent=4)
        w.write("")

        w.write("Genome Statistics (QUAST)")
        w.write("    Genome Length (bp)", self.quast['genome_length'], indent=4)
        w.write("    Genome Fraction (%)", self.quast['genome_fraction'], indent=4)
        w.write("    Genomic Features", self.quast['genomic_features'], indent=4)
        w.write("    N's per 100 kbp", self.quast['Ns_per_100_kbp'], indent=4)
        w.write("    Mismatches per 100 kbp", self.quast['mismatches_per_100_kbp'], indent=4)
        w.write("    Indels per 100 kbp", self.quast['indels_per_100_kbp'], indent=4)
        w.write("    Average Depth of Coverage", self.coverage['mean_coverage'], indent=4)

        for (l,f) in zip(self.coverage['bin_labels'], self.coverage['bin_fractions']):
            w.write(l, f, indent=8)

        w.write("5' Ns", self.consensus['N5prime'], indent=4)
        w.write("3' Ns", self.consensus['N3prime'], indent=4)
        w.write("")

        w.write("Taxonomic Composition of Assembly (LMAT)", self.lmat['taxonomic_composition'])
        w.write("")
        
        w.write("Variants in Consensus Genome (iVar)", self.ivar['variants'])
        w.write("")
        
        w.write("Variants in Read Alignment (BreSeq)", self.breseq['variants'])


class Pipeline:
    def __init__(self, dirname, sample_names):
        self.dirname = dirname
        self.sample_names = sample_names
        self.samples = [ Sample(os.path.join(dirname,s), s) for s in sample_names ]

    
    def write(self):
        w1 = SummaryHTMLWriter(os.path.join(self.dirname, 'summary.html'))

        for s in self.samples:
            w2 = SampleTextWriter(os.path.join(s.dirname, 'sample.txt'))
            w3 = SampleHTMLWriter(os.path.join(s.dirname, 'sample.html'))
        
            s.write(w1, link_to=w3)
            s.write(w2)
            s.write(w3)
            
            w2.close()
            w3.close()

        w1.close()
        

####################################################################################################


if __name__ == '__main__':
    p = Pipeline('Paper', [
        'ARTIC-NEG',
        'Iran1-ARTIC',
        'Iran1-LA',
        'Iran1-LN',
        'Iran1-RA',
        'Iran1-RN',
        'Wuhan1-ARTIC',
        'Wuhan1-LA',
        'Wuhan1-LN',
        'Wuhan1-RA',
        'Wuhan1-RN'
    ])

    p.write()


#!/usr/bin/env python3
# coding: utf-8
"""
@file: fasta_read_parse.py
@description: All-in-one parsing fasta format file
@author: Yanwei Qi
@email: qiyanwei1@genomics.cn
@last modified by: qiyanwei

change log:
    2023/09/21  first version released.
    2023/09/21  script created.
"""

import re, os, sys, time, math, gzip
import string
import json
import pandas as pd
import numpy as np
from pandas import Series, DataFrame
from mimetypes import guess_type
from functools import partial
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser

# basic stent function 
## check file
def file_checkin(infile, func):
    if os.path.exists(infile):
        print(func + ': check...' + infile + ' file existed.')
        return 1
    else:
        print(func + ': check...' + infile + ' file not found.')
        return 0
##/
## open fasta or fasta.gz or fa.gz
def open_checkin(infile, func):
    if file_checkin(infile,func):
        filetype = guess_type(infile)[1]
        openfile = partial(gzip.open, mode='rt') if filetype == 'gzip' else open
        return openfile
#/
def base_ATCG_GC_count(base):
    base = base.upper()
    A = base.count("A")
    T = base.count("T")
    C = base.count("C")
    G = base.count("G")
    B = [A, T, C, G]
    GC =  round( (C+G)/sum(B), 3)
    ATCG = {'A': A, 'T': T, 'C': C, 'G': G}
    return ATCG, GC
## Nx calculator
def Nx_reads_length(list_read_length):
    # Calculate the total length of the sequences
    total_length = sum(list_read_length)
    list_read_length.sort(reverse=True)
    # Calculate the N10, N20, N30, N40, N50, N60, N70, N80 and N90 values
    n_values = []
    for i in range(10, 100, 10):
        n_value = total_length * i / 100
        n_sum = 0
        for length in list_read_length:
            n_sum += length
            if n_sum >= n_value:
                 n_values.append(length)
                 break
    return n_values


##  LENGTH REPORTER
def read_length_statistics_reporter(list_read_length):
    dict_read_length = {}
    npLength = np.array(list_read_length)
    dict_read_length['total_read'] = len(list_read_length)
    dict_read_length['total_base'] = sum(npLength)
    dict_read_length['seq_length_min'] = min(npLength)
    dict_read_length['seq_length_max'] = max(npLength)
    dict_read_length['seq_length_avg'] = round(np.mean(npLength),3)
    Q1,Q2,Q3 = np.percentile(npLength, (25, 50, 75), method='midpoint')
    dict_read_length['seq_length_1st'] = Q1
    dict_read_length['seq_length_mid'] = Q2
    dict_read_length['seq_length_3rd'] = Q3
    # statistics the length list
    LEN_Freq_Gap = list(range(0, 10000, 1000)) + list(range(10000, 310000, 10000)) +  [400000, 500000, 600000, 1000000]
    LEN_Frequency_Distribution = pd.cut(npLength, LEN_Freq_Gap).value_counts().sort_index()
    LEN_Percentage_Distribution = LEN_Frequency_Distribution/len(list_read_length)*100
    LEN_GAP_STAT_TABLE = pd.concat([LEN_Frequency_Distribution, LEN_Percentage_Distribution], axis=1)
    LEN_GAP_STAT_TABLE = LEN_GAP_STAT_TABLE.rename(columns = {0:'Frequency',1:'Percentage'})
    print(LEN_GAP_STAT_TABLE)
    dict_read_length['seq_length_list'] = ','.join([ str(i) for i in list_read_length ])
    return dict_read_length

# fasta parse function
class fastaParse:
    'All-in-one parsing fasta format file'
    def __init__(self, args):
        self.args = args
    def fasta_stat(self):
        openfile = open_checkin(self.args.input, 'fasta_stat')
        with openfile(self.args.input)  as handle:
            list_read_length = []
            ATCG_DT = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
            for name, base in SimpleFastaParser(handle):
                read_length = len(base)
                list_read_length.append(read_length)
                atcg_dt, gc = base_ATCG_GC_count(base)
                for b in ['A', 'T', 'C', 'G']:
                    ATCG_DT[b] += atcg_dt[b]
            GC = round( (ATCG_DT['C'] + ATCG_DT['G']) / sum(ATCG_DT.values()), 3)
            Nx = Nx_reads_length(list_read_length)
            Read_LL = read_length_statistics_reporter(list_read_length)
            basic_stat_dict = {
                'GC': GC,
                'Read_length_summary': Read_LL,
                'N10-N90': Nx,
            }
        return basic_stat_dict
    def fasta_filtering(self):
        openfile = open_checkin(self.args.input, "fasta_filtering")
        length_threshold = self.args.qualified_length_threshold
        Ofasta = gzip.open(self.args.output + "_filtered_L" + str(length_threshold) +".fa.gz", 'wt')
        with openfile(self.args.input) as handle:
            list_read_length = []
            ATCG_DT = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
            for name, base in SimpleFastaParser(handle):
                read_length = len(base)
                if read_length >= int( length_threshold ):
                    list_read_length.append(read_length)
                    record = ">%s\n%s\n" % (name, base)
                    Ofasta.write( record.format('fasta'))
                    atcg_dt, gc = base_ATCG_GC_count(base)
                    for b in ['A', 'T', 'C', 'G']:
                        ATCG_DT[b] += atcg_dt[b]
                else:
                    pass
            GC = round( (ATCG_DT['C'] + ATCG_DT['G']) / sum(ATCG_DT.values()), 3)
            Nx = Nx_reads_length(list_read_length)
            Read_LL = read_length_statistics_reporter(list_read_length)
            basic_stat_dict = {
                'GC': GC,
                'N10-N90': Nx,
                'Read_length_summary': Read_LL,
            }
            return basic_stat_dict
    def fasta_cutting(self):
        cut_head_len = args.cut_head_length
        cut_tail_len = args.cut_tail_length
        Ofasta = gzip.open(self.args.output + "_cut_Head" + str(cut_head_len) + 'bp_Tail' + str(cut_tail_len) +"bp.fa.gz", 'wt')
        openfile = open_checkin(self.args.input, "fastq_cut")
        with openfile(self.args.input) as handle:
            list_read_length = []
            for name, base in SimpleFastaParser(handle):
                read_length = len(base)
                if ( cut_head_len + cut_tail_len ) <= read_length and cut_tail_len > 0:
                    seq1 = base[ cut_head_len: -cut_tail_len ]
                    list_read_length.append( len(seq1) )
                    record = ">%s\n%s\n" % (name, seq1)
                    Ofasta.write( record.format('fasta'))
                elif ( cut_head_len + cut_tail_len ) <= read_length and cut_tail_len == 0:
                    seq1 = base[ cut_head_len: ]
                    list_read_length.append( len(seq1) )
                    record = ">%s\n%s\n" % (name, seq1)
                    Ofasta.write( record.format('fasta'))
                else:
                    pass
        dict_read_length = read_length_statistics_reporter(list_read_length)
        return dict_read_length
    def fasta_trimming(self):
        head_homo_len = int(self.args.trim_head_homopolybase)
        def head_homopolybase_detector(base, head_homo_len):
            for b in ['A', 'G', 'C', 'T']:
                homopolybases = b*int(head_homo_len)
                if base.startswith(homopolybases):
                    patterns = re.compile(r"%s{%d,}" % (b, head_homo_len))
                    head_homopolybase_length = len(patterns.findall( base )[0])
                    return head_homopolybase_length
                else:
                    return 0
        tail_homo_len = int(self.args.trim_tail_homopolybase)
        def tail_homopolybase_detector(base, tail_homo_len):
            for b in ['A', 'G', 'C', 'T']:
                homopolybases = b*int(tail_homo_len)
                if base.endswith(homopolybases):
                    patterns = re.compile(r"%s{%d,}" % (b, tail_homo_len))
                    tail_homopolybase_length = len(patterns.findall( base )[-1])
                    return tail_homopolybase_length
                else:
                    return 0
        openfile = open_checkin(self.args.input, "fasta_trimming")
        Ofasta = gzip.open(self.args.output + "_trim_end_homobase_Head" + str(head_homo_len) + '_Tail' + str(tail_homo_len) +".fa.gz", 'wt')
        with openfile(self.args.input) as handle:
            list_read_length = []
            long_homo_ends_iter, long_homo_head_iter, long_homo_tail_iter = 0, 0, 0
            for name, base in SimpleFastaParser(handle):
                read_length = len(base)
                if (head_homo_len + tail_homo_len) <= read_length and tail_homo_len > 0:
                    head_homopolybase_length = head_homopolybase_detector(base, head_homo_len)
                    tail_homopolybase_length = tail_homopolybase_detector(base, tail_homo_len)
                    if tail_homopolybase_length > 0:
                        seq1 = base[ head_homopolybase_length : -tail_homopolybase_length]
                        list_read_length.append( len(seq1) )
                        record = ">%s\n%s\n" % (name, seq1)
                        Ofasta.write( record.format('fasta'))
                        long_homo_tail_iter += 1
                        if head_homopolybase_length > 0:
                            long_homo_ends_iter += 1
                            long_homo_head_iter += 1
                    else:
                        seq1 = base[ head_homopolybase_length: ]
                        list_read_length.append( len(seq1) )
                        record = ">%s\n%s\n" % (name, seq1)
                        Ofasta.write( record.format('fasta'))
                        if head_homopolybase_length > 0:
                            long_homo_head_iter += 1
                elif (head_homo_len + tail_homo_len) <= read_length and tail_homo_len == 0:
                    head_homopolybase_length = head_homopolybase_detector(base, head_homo_len)
                    seq1 = base[ head_homopolybase_length: ]
                    list_read_length.append( len(seq1) )
                    record = ">%s\n%s\n" % (name, seq1)
                    Ofasta.write( record.format('fasta'))
                    if head_homopolybase_length > 0:
                        long_homo_head_iter += 1
                else:
                    read_number_iter += 1
            dict_read_length = read_length_statistics_reporter(list_read_length)
            dict_read_length.update({ "homobase_in_longread_head_number": long_homo_head_iter, 
                                      "homobase_in_longread_tail_number": long_homo_tail_iter, 
                                      "homobase_in_longread_ends_number": long_homo_ends_iter})
            return dict_read_length
    def fasta_sampling(self):
        """
            Row = R5000 # row5000lines
            overLen = OL5000 # over length 5000bp
            lessLen = LL5000 # less length 5000bp
            genome_size, depth = 3000000000,30 # genome_size,depthx
        """
        dict_sampling_stat = {}
        if int(self.args.pickup_number_of_reads) > 0:
            list_read_length = []
            openfile = open_checkin(self.args.input, "fastq_sampling")
            with openfile(self.args.input)  as handle:
                pickup_lines_number_read = int(self.args.pickup_number_of_reads)
                Ofasta = gzip.open(self.args.output + '_pick_' + str(pickup_lines_number_read) + "row_reads.fa.gz", 'wt')
                read_line_iter = 0
                for name, base in SimpleFastaParser(handle):
                    read_line_iter += 1
                    if read_line_iter <= pickup_lines_number_read:
                        list_read_length.append( len(base) )
                        record = ">%s\n%s\n" % (name, base)
                        Ofasta.write( record.format('fasta'))
                    else:
                        Ofasta.close()
                        pass
            dict_read_length = read_length_statistics_reporter(list_read_length)
            dict_sampling_stat.update( { "sampling_rows_read_summary": dict_read_length})
        if int(self.args.pickup_over_length_of_reads) > 0:
            list_read_length = []
            openfile = open_checkin(self.args.input, "fastq_sampling")
            with openfile(self.args.input)  as handle:
                pickup_over_length_read = int(self.args.pickup_over_length_of_reads)
                Ofasta = gzip.open(self.args.output + '_over_' + str(pickup_over_length_read) + "bp_reads.fq.gz", 'wt')
                for name, base in SimpleFastaParser(handle):
                    read_length = len(base)
                    if read_length >= pickup_over_length_read:
                        list_read_length.append( len(base) )
                        record = ">%s\n%s\n" % (name, base)
                        Ofasta.write( record.format('fasta'))
                    else:
                        pass
                Ofasta.close()
            dict_read_length = read_length_statistics_reporter(list_read_length)
            dict_sampling_stat.update( {'sampling_long_read_summary':dict_read_length} )
        if int(self.args.pickup_less_length_of_reads) > 0:
            list_read_length = []
            openfile = open_checkin(self.args.input, "fastq_sampling")
            with openfile(self.args.input)  as handle:
                pickup_less_length_read = int(self.args.pickup_less_length_of_reads)
                Ofasta = gzip.open(self.args.output + '_less_' + str(pickup_less_length_read) + "bp_reads.fq.gz", 'wt')
                for name, base in SimpleFastaParser(handle):
                    read_length = len(base)
                    if read_length < pickup_less_length_read:
                        list_read_length.append( len(base) )
                        record = ">%s\n%s\n" % (name, base)
                        Ofasta.write( record.format('fasta'))
                    else:
                        pass
                Ofasta.close()
            dict_read_length = read_length_statistics_reporter(list_read_length)
            dict_sampling_stat.update( {'sampling_short_read_summary':dict_read_length} )
        if  float(self.args.downsample_genomesize) > 0 and float(self.args.downsample_depthx) > 0:
            list_read_length = []
            openfile = open_checkin(self.args.input, "fastq_sampling")
            with openfile(self.args.input)  as handle:
                genomesize, depthx = float(self.args.downsample_genomesize), float(self.args.downsample_depthx)
                Ofasta = gzip.open(self.args.output + '_depth_' + str( self.args.downsample_depthx ) + "x.fq.gz", 'wt')
                base_size = int(genomesize * depthx)
                base_size_iter = 0
                for name, base in SimpleFastaParser(handle):
                    read_length = len(base)
                    if base_size_iter <= base_size:
                        record = ">%s\n%s\n" % (name, base)
                        Ofasta.write( record.format('fasta'))
                        base_size_iter += read_length
                        list_read_length.append( len(base) )
                    else:
                        Ofasta.close()
                        break
            dict_read_length = read_length_statistics_reporter(list_read_length)
            dict_sampling_stat.update( {'downsampling_fold_read_summary':dict_read_length} )
        return dict_sampling_stat
    def fasta_kmer(self):
        kmersize = self.args.kmer_size
        output = self.args.output
        pwd_config_file = os.path.realpath(__file__)
        meryl = '/'.join(pwd_config_file.split('/')[:-1]) + '/tools/meryl'
        if self.args.input.endswith('gz'):
            os.system('{} count  k={} {} output {}.meryl'.format(meryl, str(kmersize), self.args.input,  output) )
            os.system('{} print {}.meryl |sort -k2nr > {}_kmer_{}_freq.txt'.format(meryl, output, output, str(kmersize)))
            os.system('rm -f {}.fasta'.format( output ))
            os.system('rm -rf {}.meryl'.format( output ))
        else:
            print('please determining the input file suffix is fastq or fq or fq.gz!')
    def fasta_homopolymer(self):
        from collections import Counter
        import seaborn as sns
        import matplotlib.pyplot as plt
        from matplotlib.ticker import FixedLocator, MaxNLocator
        def homopolymer_freq(base, homopolymer_len, homopolymer_dict):
            POLY = []
            ATCG = ['A', 'G', 'C', 'T']
            for b in ATCG:
                patterns = re.compile(r"%s{%d,}" % (b, homopolymer_len))
                POLY.extend( patterns.findall( base ))
            for k, v in dict(Counter(POLY)).items():
                if k in homopolymer_dict.keys():
                    homopolymer_dict[k] += v
                else:
                    homopolymer_dict[k] = v
            return homopolymer_dict
        hpmer_len_min = self.args.homopolymer_length_min # minimum length of homopolymer
        hpmer_len_max = self.args.homopolymer_length_max # maximum length of homopolymer
        homopolymer_dict = {}
        openfile = open_checkin(self.args.input, "fasta_homopolymer")
        with openfile(self.args.input)  as handle:
             for name, base in SimpleFastaParser(handle):
                 homopolymer_dict = homopolymer_freq(base, hpmer_len_min, homopolymer_dict)
        hpm = sorted(homopolymer_dict.keys())
        frq = {}
        for p in hpm:
            if len(p) < hpmer_len_max and p in homopolymer_dict.keys():
                 frq["Poly"+p[0]+'_'+str(len(p))] = homopolymer_dict[p]
            elif len(p) >= hpmer_len_max and p in homopolymer_dict.keys():
                 hpkey = "Poly"+p[0]+'_'+str(hpmer_len_max)
                 if hpkey in frq.keys():
                     frq[hpkey] += homopolymer_dict[p]
                 else:
                     frq[hpkey] = homopolymer_dict[p]
            else:
                 print(p)
        qrf = []
        out_prefix = self.args.output + '_homopolymer_range'+str(hpmer_len_min) + '-' + str(hpmer_len_max) + 'bp_freq.txt'
        homofile = open(out_prefix + '.txt', 'w')
        for p, n in frq.items():
            qrf.append( p.split("_") +[ n ])
            homofile.write('\t'.join(p.split("_") +[ str(n) ]) + '\n')
        HPmer = pd.DataFrame(data=qrf, columns=['Base_type', 'Homopolymer_len', 'Homopolymer_freq'])
        fig = plt.subplots(figsize=(20, 20))
        fig = sns.relplot(data=HPmer, x="Homopolymer_len", y="Homopolymer_freq", hue="Base_type", kind='line', )
        ax = plt.gca()
        ax.xaxis.set_major_locator(MaxNLocator(steps=[10]))
        ax.xaxis.set_minor_locator(FixedLocator(range(hpmer_len_min-5, hpmer_len_min+5)))
        fig.fig.set_dpi(300)
        plt.savefig( out_prefix + ".png")
        return { 'PolyBase': ','.join( frq.keys() ), 'HomoFreq': ','.join([ str(i) for i in frq.values() ]) }
 
         
import argparse

parser = argparse.ArgumentParser(description='Long reads sequencing QC tool')
parser.add_argument('--input', 
                    required='True',
                    help='raw data: fasta/fa.gz')
parser.add_argument('--stat',
                    default=1,
                    type=int,
                    help='statistic the length of reads in fasta file.')
parser.add_argument('--filtering',
                    default=0,
                    type=int,
                    help='filter short or low-quality read or not: 1 or 0.')
parser.add_argument('--qualified_length_threshold',
                    default=1000,
                    type=int,
                    help='filtered quality of read (default: 1000).')
parser.add_argument('--cut_read_ends',
                    default=0,
                    type=int,
                    help='cut the head or tail segment of read or not: 1 or 0.')
parser.add_argument('--cut_head_length',
                    default=100,
                    type=int,
                    help='cut the head segment of read (default: 100) under --cut_double_ends_length with 1, only cut head segment setting --cut_tail_length to 0.')
parser.add_argument('--cut_tail_length',
                    default=100,
                    type=int,
                    help='cut the head segment of read (default: 100) under --cut_double_ends_length with 1, only cut tail segment setting --cut_head_length to 0.')
parser.add_argument('--trim_ends_homopolybase',
                    default=0,
                    type=int,
                    help='trim the long homopoly bases segment of read or not: 1 or 0.')
parser.add_argument('--trim_head_homopolybase',
                    default=100,
                    type=int,
                    help='trim the head homopolybase of read (default: 100) under --trim_ends_homopolybase with 1, only cut head segment setting --trim_tail_homopolybase to 0.')
parser.add_argument('--trim_tail_homopolybase',
                    default=100,
                    type=int,
                    help='trim the head homopolybase of read (default: 100) under --trim_ends_homopolybase with 1, only cut tail segment setting --trim_head_homopolybase to 0.')
parser.add_argument('--downsampling',
                    default=0,
                    type=int,
                    help='downsampling or not: 1 or 0.')
parser.add_argument('--pickup_number_of_reads',
                    default=0,
                    type=int,
                    help='pick up the number of read, parameter under --downsampling with 1.')
parser.add_argument('--pickup_over_length_of_reads',
                    default=0,
                    type=int,
                    help='pick up the minimum length of read, parameter under --downsampling with 1.')
parser.add_argument('--pickup_less_length_of_reads',
                    default=0,
                    type=int,
                    help='pick up the maximum length of read, parameter under --downsampling with 1.')
parser.add_argument('--downsample_genomesize',
                    type=str,
                    help='the genome size of downsample data followed by sample species, binding with --downsample_depthx under --downsample with 1.')
parser.add_argument('--downsample_depthx',
                    type=str,
                    help='the depth (x) fold of downsample data custom-made by user, binding with --downsample_genomesize under --downsample with 1.')
parser.add_argument('--kmer',
                    default=0,
                    type=int,
                    help='homopolymer analysis or not: 1 or 0.')
parser.add_argument('--kmer_size',
                    default=6,
                    type=int,
                    help='observed the minimum length of kmer(default: 6).')
parser.add_argument('--homopolymer',
                    default=0,
                    type=int,
                    help='homopolymer analysis or not: 1 or 0.')
parser.add_argument('--homopolymer_length_min',
                    default=5,
                    type=int,
                    help='observed the minimum length of homopolymer(default: 5).')
parser.add_argument('--homopolymer_length_max',
                    default=50,
                    type=int,
                    help='observed the maximum length of homopolymer(default: 50).')
parser.add_argument('--output',
                    required='True',
                    default='Fasta_output',
                    help='the sign of output file name, following by sample name is better.')
parser.add_argument("--xlsx",
                    default='LRstat.xlsx', 
                    type=str,
                    help='output table file name, filetype=XLSX')
parser.add_argument('--json',
                    default='LRstat.json', 
                    type=str,
                    help='output table file name, filetype=JSON')
args = parser.parse_args()
print(args)

fp = fastaParse(args)
if args.stat:
    fa_basic_stat = fp.fasta_stat()
    fa_stat = {'raw_long_reads_stat': fa_basic_stat}
if args.cut_read_ends:
    fa_cutting_stat = fp.fasta_cutting()
    fa_stat.update({ 'cutted_long_reads_stat': fa_cutting_stat })
if args.trim_ends_homopolybase:
    fa_trim_stat = fp.fasta_trimming()
    fa_stat.update({ 'trimmed_long_reads_homopolymer_stat': fa_trim_stat })
if args.filtering:
    fa_filter_stat = fp.fasta_filtering()
    fa_stat.update({ 'filtered_long_reads_statistics': fa_filter_stat })
if args.downsampling:
    fa_sampling_stat = fp.fasta_sampling()
    fa_stat.update({ 'sampled_long_reads_stat': fa_sampling_stat })
if args.kmer:
    fp.fasta_kmer()
if args.homopolymer:
    fa_homo_stat = fp.fasta_homopolymer()
    fa_stat.update({'homopolymer_frequency_stat': fa_homo_stat })

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)

with open(args.output+'_'+args.json, 'w') as file:
    json.dump(fa_stat, file, indent='\t',
              separators=(', ', ': '), ensure_ascii=False,
              cls=NpEncoder)



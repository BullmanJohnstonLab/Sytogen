#!/usr/bin/env python

__author__ = 'Gianmarco Piccinno'
__version__ = '1.0.0'
__date__ = '26th February 2021'

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# SynPlART: Syngenic Plasmid Artificer
# Authors: Gianmarco Piccinno,
#          Serena Manara,
#          Nicola Segata,
#          Christopher Johnston
# Version: 1.0.0
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


import argparse as ap
import Bio
from Bio import Seq, SeqIO, SeqFeature
from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition
from Bio.SeqUtils import MeltingTemp as mt
from collections import OrderedDict
import copy
import configparser
import datetime
from dna_features_viewer import BiopythonTranslator, CircularGraphicRecord
import functools
from functools import reduce
import glob
import itertools
from itertools import product
import logging
import math
import matplotlib.backends.backend_pdf
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import multiprocessing
from multiprocessing import Pool
from multiprocessing.pool import ThreadPool
from multiprocessing import Process, Value, Array
import operator
from operator import mul
import os
import pandas as pd
import pickle
from pydna.assembly import Assembly
from pydna.design import assembly_fragments
from pydna.design import primer_design
from pydna.dseqrecord import Dseqrecord
from pydna.amplify import pcr
from pydna.tm import dbd_program
import random
import re
import sre_yield
import subprocess
import sys
import textwrap
import time
from types import SimpleNamespace
import warnings

if sys.version_info[0] < 3:
    sys.exit(
        f'SynToGen {__version__} requires Python 3, your current Python version is {sys.version_info[0]}.{sys.version_info[1]}.{sys.version_info[2]}')

camth = """rec_seq	Methyl-Modification introduced (5'-3')	meth_base	comp_meth_base	comp_meth_type	meth_type	availability
ATCGAT	ATCGm6AT	5	2	6	6	M.ClaI (Takara Bio); M.BseCI (Minotech Biotechnology)
AAGCTT	m6AAGCTT	1	6	6	6	M.HindIII (Takara Bio)
GGATCC	GGATm4CC	5	2	4	4	M.BamHI (Takara Bio); M.BamHI (New England Biolabs)
GAATTC	GAm6ATTC	3	4	6	6	M.EcoRI (Takara Bio); M.EcoRI (New England Biolabs); M.EcoRI (Nippon Gene); 
GGATG	GGm6ATG	3	3	6	6	M3.BstF5I (SibEnzyme)
GCNGC	Gm5CNGC	2	4	5	5	M.Fsp4HI (SibEnzyme)
AGCT	AGm5CT	3	2	5	5	M.AluI (Takara Bio); M.AluI (New England Biolabs)
TCGA	TCGm6A	4	1	6	6	M.TaqI (New England Biolabs)
GATC	Gm6ATC	2	3	6	6	M.EcoKDam (New England Biolabs)
GGCC	GGm5CC	3	2	5	5	M.HaeIII (Takara Bio); M.HaeIII (New England Biolabs)
GCGC	Gm5CGC	2	3	5	5	M.HspAI (SibEnzyme)
CCGG	m5CCGG	1	4	5	5	M.MspI (New England Biolabs)
CCGG	Cm5CGG	2	3	5	5	M.HpaII (Takara Bio); M.HpaII (New England Biolabs)
GC	Gm5C	2	1	5	5	M.CviPI (New England Biolabs)
CG	m5CG	1	2	5	5	M.SssI (New England Biolabs); M.SssI (ThermoFischer); M.SssI (Zymo Research)
A	m6A	1	1	6	6	M.EcoGII (New England Biolabs)"""

synpl_header = """+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+                     WELCOME TO SyToGen!                           +
+ You have run SynPlEST successfully!                               +
+ .synpl format specifications                                      +
+ If ORIGINAL|In_CDS:0 --> The slice is in a non-coding region      +
+ If ORIGINAL|In_CDS:1 --> The slice is in a coding region          +
+ Each block is made of                                             +
+ \--| Indicates the annotation for the CDS                         +
+ \--> Indicates the CDS in triplets                                +
+ \*-> Indicates the amino acid sequence of the CDS                 +
+ This block can be replicated n times, where n is the number of    +
+ overlapping CDS in the slice                                      +
+ IMPORTANT:                                                        +
+ ALL THE SEQUENCE ARE SHOWN AS ON THE + STRAND OF THE PLASMID!!!   +
+ The indication before codons (+/-) indicates on which strand      +
+ of the plasmid the CDS is placed                                  +
+ IMPORTANT:                                                        +
+ All the indexes start from 0.                                     +
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
_____________________________________________________________________
_____________________________________________________________________\n"""

synpl_map_header = """+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+                     WELCOME TO SyToGen!                           +
+ This is temporary .synpl output from SyToGen.                     +
+ This file is intended to provide access to the sequence of the    +
+ current slice.                                                    +
+ Change case (from lower to upper case) of bases in the            +
+ >SelectTargetPositions slice to indicate to SyToGen               +
+ which positions to generalize in the current slice.               +
+ Please do not edit the >Original slice.                           +
+ Modifications to the  >Original slice will be discarded and       +
+ will not be considered in the following procedure.                +
_____________________________________________________________________
_____________________________________________________________________\n
"""

report_html_ = """<!DOCTYPE html>
<html lang="en">
<head>
  <title>SyToGen Suite Report</title>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.4.1/css/bootstrap.min.css">
  <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.4.1/jquery.min.js"></script>
  <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.4.1/js/bootstrap.min.js"></script>
</head>
<body>

<div class="jumbotron text-center">
  <h1>SyToGen Suite Report</h1>
  <p>Piccinno G., Manara S., Johnston C., Segata N.</p> 
</div>
  
<div class="container" style="width: 75%;">
  <div class="row">
    <div class="col-lg-5 pre-scrollable" data-spy="scroll" data-offset="80" style="font-family:monospace; font-size:12px; max-height: 75vh; height: calc(100vh - 200px);
  overflow-y: auto;">
      <h3>Plasmid/minicircle information</h3>
      <p>Plasmid annotations:</p>
      <p>__ADD_ANNOTATIONS__</p>
      <p>Target motifs found on the original plasmid/minicircle:</p>
      <p>__ADD_OR_TARGET_MOTIFS</p>
      <p>Target motifs found on the proposed candidate syngenic plasmid plasmid/minicircle:</p>
      <p>__ADD_CAND_TARGET_MOTIFS</p>
      <p>Modified positions:</p>
      <p>__ADD_MODIFIED_POSITIONS</p>
    </div>
    <div class="col-lg-7 pre-scrollable" data-spy="scroll" data-offset="80" style="font-family:monospace; font-size:12px; max-height: 75vh; height: calc(100vh - 200px);
  overflow-y: auto;">
      
    <style>
            red {color:#FF4136;}
            green {color:#66bd63;}
            yellow {color:#FFDC00;}
            navy {color:#001f3f;}
            blue {color:#0074D9;}
            aqua {color:#7FDBFF;}
            teal {color:#39CCCC;}
            olive {color:#2ECC40;}
            lime {color:#01FF70;}
            orange {color:#f46d43;}
            maroon{color:#85144b;}
            fuchsia {color:#F012BE;}
            purple {color:#B10DC9;}
            black {color:#111111;}
            gray {color:#AAAAAA;}
            silver {color:#DDDDDD;}
            darkred {color:rgb(197,27,125)}
            darkgreen {color:rgb(0, 110, 0)}
            darkblue {color:rgb(77,146,33)}
            darkyellow {color:rgb(190, 170, 0)}
    </style>
      <h3>Plasmid alignment</h3>
      <p>__ADD_PLASMID_ALIGNMENT__</p>
    </div>

  </div>
</div>

</body>
</html>
"""

# SyToGen Sys

def info(m, init_new_line=True, exit=False, exit_value=0):
    if init_new_line:
        sys.stdout.write('\n')

    sys.stdout.write(f'{m}')
    sys.stdout.flush()

    if exit:
        sys.exit(exit_value)

def error(m, init_new_line=True, exit=True, exit_value=1):
    if init_new_line:
        sys.stderr.write('\n')

    sys.stderr.write(f'Error: {m}\n')
    sys.stderr.flush()

    if exit:
        sys.exit(exit_value)

def progress_status(progress, total, message='Current process'):

    barLength, status = 50, ""
    progress = float(progress) / float(total)
    if progress >= 1.:
        progress, status = 1, "\r\n"
    block = int(round(barLength * progress))
    text = "\r {}: [{}] {:.0f}% {}".format(
        message, "#" * block + "-" * (barLength - block), round(progress * 100, 0),
        status)
    sys.stdout.write(text)
    sys.stdout.flush()

# SyToGen Utils

def reverse_complement(input_sequence):

    complement_ = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    if not (set(input_sequence.upper()) - set(['A', 'T', 'C', 'G'])):
        return ''.join([complement_[el] for el in input_sequence.upper()[::-1]])
    else:
        sys.exit('The input sequence is not a genetic sequence [only canonical bases]!')

def translate(input_sequence, codon_usage_table):
    translated_sequence = ''.join([str(codon_usage_table[x]['Translation']) for x in re.split(
        r'(\w{3})', input_sequence) if x != ''])
    return translated_sequence

def irange(i, j, by=1, l=0):
    '''

    This function computes the range for a circular sequence.
    This function takes as input three integers:

    i: starting position of the range (included)
    j: ending position of the range (excluded)
    l: length of the sequence. The starting position of the sequence is supposed to be 0.

    The assumption is that the ranges are indicated according to c/python language convention:
    - numbering starts with 0
    - last position of a sequence is the actual length of the sequence reduced by one
    - ranges are made with first position included, second excluded

    The output is a list of lists with:
    - a single element (list) that corresponds to the actual range if the first position is lower than the second
    - two elements (lists) if the first position is higher than the second, indicating that we are considering
    the section of the circular sequence that goes from before the end of the sequence, across the junction
    to the starting portion of the sequence

    Example:

    >>> ss = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
    >>> intervals = [[10, 1], [0, 2], [8, 9], [10, 11]]

    >>> intervals = [irange(10, 2, len(ss)), irange(-2, 3, len(ss)), irange(5, 9, len(ss)), irange(8, 11, len(ss))]
    >>> info(intervals)
    [[[10, 11, 12, 13, 14, 15], [0, 1]], [[14, 15], [0, 1, 2]], [[5, 6, 7, 8]], [[8, 9, 10]]]

    '''

    i_ = i % l
    j_ = j % l

    if i_ <= j_:
        return list(range(i_, j_, by))
    else:
        return list(range(i_, l, by)) + list(range(0, j_, by))

def handle_ranges(list_ranges_):
    '''
    This function defines overlapping ranges and return a list of non overlapping ranges, in which the union of those overlapping is considered.
    In SyToGen, each irange must be elongate on the left and on the right by the length of the longest target motif. This must be done when defining
    the inputs of the irange function irange(i-m, j+m, l). This permits to retrieve from this function slices' positions that are at most contiguous
    and non-overlapping. Contigous slices preserve the correctness of the slice evaluation because of the two tails introduced by the irange.

    In this respect this function permits to define the positions of each bases of the slice on the plasmid.
    To get the irange take the first position and the last + 1 for each list in the list.

    >>> for el in list:
    >>>     info( irange(el[0][0], el[-1][-1]+1)) <-- This is the correct irange

    Each slice is made by:

    ---n [[], [], []] n--- :
    -two external n tails from the original sequence, where n is the length of the longest target motif
    -a list of lists of inner positions (from handle_ranges)

    The input is a list of iranges.

    Example:

    >>> ss = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
    >>> intervals = [[10, 1], [0, 2], [8, 9], [10, 11]]

    >>> intervals = [irange(10, 2, len(ss)), irange(-2, 3, len(ss)), irange(5, 9, len(ss)), irange(8, 11, len(ss))]
    >>> info(intervals)
    [[[10, 11, 12, 13, 14, 15], [0, 1]], [[14, 15], [0, 1, 2]], [[5, 6, 7, 8]], [[8, 9, 10]]]

    >>> out_intervals = handle_ranges(intervals)

    >>> info(out_intervals)
    [[[8, 9, 10, 11, 12, 13, 14, 15], [0, 1, 2]], [[5, 6, 7, 8]]]

    '''

    list_ranges = copy.deepcopy(list_ranges_)
    set_tot_overlap = [list_ranges[0]]

    for i in range(1, len(list_ranges)):
        overlap = False
        for j in range(len(set_tot_overlap)):

            for k in range(len(list_ranges[i])):
                for l in range(len(list_ranges[j])):
                    if set.intersection(set(list_ranges[i][k]), set(list_ranges[j][l])):
                        set_tot_overlap[j][l] = list(
                            set.union(set(list_ranges[i][k]), set(list_ranges[j][l])))
                        overlap = True

        if not overlap:
            set_tot_overlap.append(list_ranges[i])

    return set_tot_overlap

def findpatterns(sequence, patterns_, circular=True):
    '''
    This function takes as input the sequence of the plasmid, and the list of target motifs (a list of strings)
    and returns a dictionary with keys each target motifs.
    For each key in the output a list of the following elements is provided:

    -the sequence on the plasmid corresponding to the target motif
    -the starting position of the target (included)
    -the ending position of the target (excluded)
    -the strand on which it is found (1 indicates the forward, -1 indicates the reverse complement)
    -junction/inner if the sequence is found on the 0 position or inside the actual sequence of the plasmid
    '''

    patterns = []

    def plan_ambiguity(pattern):
        val = Bio.Data.IUPACData.ambiguous_dna_values
        re_pattern = ''
        for el in pattern:
            re_pattern = re_pattern + '[' + val[el] + ']'
        return re_pattern

    for pattern in patterns_:
        patterns.append(plan_ambiguity(pattern))
    patts = {}

    if circular:

        for i in range(len(patterns_)):

            tmp_fw = plan_ambiguity(patterns_[i])
            tmp_rv = plan_ambiguity(
                str(Seq(patterns_[i]).reverse_complement()))

            data = sequence + sequence[:len(patterns_[i])]

            patts[patterns_[i]] = [(data[j.start() % len(sequence):j.end() % len(sequence)], j.start(
            ) % len(sequence), j.end() % len(sequence), 1) for j in re.finditer(tmp_fw, data)]

            patts[patterns_[i]] += [(data[(len(data)-j.end()) % len(sequence):(len(data)-j.start()) % len(sequence)], (len(data)-j.end()) % len(
                sequence), (len(data)-j.start()) % len(sequence), -1) for j in re.finditer(tmp_fw, str(Seq(data).reverse_complement()))]
    else:
        for i in range(len(patterns)):

            tmp_fw = plan_ambiguity(patterns_[i])
            tmp_rv = plan_ambiguity(
                str(Seq(patterns_[i]).reverse_complement()))

            data = sequence

            patts[patterns_[i]] = [(data[j.start() % len(sequence):j.end() % len(sequence)], j.start(
            ) % len(sequence), j.end() % len(sequence), 1) for j in re.finditer(tmp_fw, data)]

            patts[patterns_[i]] += [(data[(len(data)-j.end()) % len(sequence):(len(data)-j.start()) % len(sequence)], (len(data)-j.end()) % len(
                sequence), (len(data)-j.start()) % len(sequence), -1) for j in re.finditer(tmp_fw, str(Seq(data).reverse_complement()))]

    return patts

def punctuate_targets(f_patts, l, max_len, is_pos_list = False):
    '''
    This function takes as input the target motifs mapped against the sequence (output from findpatterns) and a genetic sequence (better if it is a slice)
    and keep only the position on the sequence to which at least one non-ambiguous basis of the mapped target motifs is found

    '''

    iranges = []
    non_ambiguous_positions = []

    if not is_pos_list:
        n_poss = {}
        for key in f_patts.keys():
            for el in f_patts[key]:
                tmp = irange(el[1]-max_len, el[2]+max_len, l=l)
                iranges.append(tmp)
                non_ambiguous_positions.append([])
                for i in [el-max_len for el in range(max_len, len(tmp)-max_len)]:
                    if el[3] == -1:
                        tmp_base = key[::-1][i]
                    else:
                        tmp_base = key[i]
                    if tmp_base in ['A', 'T', 'C', 'G']:
                        non_ambiguous_positions[-1].append(tmp[i+max_len])

    else:
        for el in f_patts:
            tmp = irange(el-max_len, el+max_len+1, l=l)
            iranges.append(tmp)
            non_ambiguous_positions.append(el)

    print("\n--> HERE1: ", sorted(iranges))
            
    set_tot_overlap = [iranges[0]]
    set_tot_overlap_non_ambiguous = [[non_ambiguous_positions[0]]]

    for i in range(1, len(iranges)):
        overlap = False
        for j in range(len(set_tot_overlap)):
            if set.intersection(set(iranges[i]), set(set_tot_overlap[j])):
                set_tot_overlap[j] = sorted(list(
                    set.union(set(iranges[i]), set(set_tot_overlap[j]))))
                if not non_ambiguous_positions[i] in set_tot_overlap_non_ambiguous[j]:
                    set_tot_overlap_non_ambiguous[j].append(non_ambiguous_positions[i])
                overlap = True
        if not overlap:
            set_tot_overlap.append(iranges[i])
            set_tot_overlap_non_ambiguous.append([non_ambiguous_positions[i]])
            
    print("\n--> HERE2: ", sorted(set_tot_overlap), " ======" , sorted(set_tot_overlap_non_ambiguous))

    for i in range(len(set_tot_overlap)):
        
        if ((l-1 in set_tot_overlap[i]) & (0 in set_tot_overlap[i])):
            end = set_tot_overlap[i][0]
            start = set_tot_overlap[i][-1]
            for j in range(1, len(set_tot_overlap[i])):
                if (set_tot_overlap[i][j] - set_tot_overlap[i][j-1] > 1):
                    end = set_tot_overlap[i][j-1] + 1
                    start = set_tot_overlap[i][j]
            set_tot_overlap[i] = irange(start, end, l=l)
        else:
            set_tot_overlap[i] = sorted(set_tot_overlap[i])

    sorted_indexes = sorted(range(len(set_tot_overlap)),
                            key=lambda k: set_tot_overlap[k])
    set_tot_overlap = sorted(set_tot_overlap)
    set_tot_overlap_non_ambiguous = [
        sorted(set_tot_overlap_non_ambiguous[i]) for i in sorted_indexes]
    
    print("\n--> HERE3: ", sorted(set_tot_overlap), " ======" , sorted(set_tot_overlap_non_ambiguous))

    return (set_tot_overlap, set_tot_overlap_non_ambiguous)

def extract_subsequence(sequence, irange_):
    '''
    irange must be output form irange function
    '''

    subsequence = ""
    for index in irange_:
        subsequence += sequence[index]
    return subsequence

def combn(input_list, k):
    '''
    The input list must be a list of lists, where each sublist is a list of positions on the plasmid, e.g.

    a = [list(range(1,4)), list(range(5,7)), list(range(8,12))]

    '''

    plan_list = sum(input_list, [])
    tmp_combs = set(list(itertools.combinations(plan_list, k)))

    if len(input_list) >= 2:
        for el in input_list:
            tmp_combs = set.difference(
                tmp_combs, list(itertools.combinations(el, k)))

    for el1 in tmp_combs:
        tmp = []
        for el2 in input_list:
            if set.intersection(set(el2), set(el1)) == set():
                tmp.append(False)
            else:
                tmp.append(True)
        if set(tmp) != set([True]):
            tmp_combs = set.difference(tmp_combs, [el1])

    return sorted(list(tmp_combs))

def meth_tr_mapping(dict_rm_tr, slice_c, slice_range):

    meth_obj = {}
    for key in dict_rm_tr.keys():
        # 'rec_seq': 'ATCGAT', 'meth_base': 4, 'comp_meth_base': 1, 'meth_type': 5, 'comp_meth_type': 5
        tmp_ = (dict_rm_tr[key]['rec_seq'], str(dict_rm_tr[key]['meth_base']) + '+' + str(dict_rm_tr[key]['meth_type']), str(
            dict_rm_tr[key]['comp_meth_base']) + '-' + str(dict_rm_tr[key]['comp_meth_type']))
        ff = findpatterns(
            slice_c, [dict_rm_tr[key]['rec_seq']], circular=False)
        val = ff[list(ff.keys())[0]]

        keys = []
        if val != []:
            for el in val:
                if el[-1] == -1:
                    tmp_key = (str(el[2] - dict_rm_tr[key]['meth_base'] - 1) + '-' + str(dict_rm_tr[key]['meth_type']), str(el[2] - dict_rm_tr[key]
                                                                                                                            ['comp_meth_base'] - 1) + '+' + str(dict_rm_tr[key]['comp_meth_type']))
                else:
                    tmp_key = (str(dict_rm_tr[key]['meth_base']+el[1]) + '+' + str(dict_rm_tr[key]['meth_type']), str(dict_rm_tr[key]['comp_meth_base'] +
                                                                                                                      el[1]) + '-' + str(dict_rm_tr[key]['comp_meth_type']))
                if not set(tmp_key) in keys:
                    meth_obj[tuple(sorted(list(tmp_key)))] = [tmp_]
                    keys.append(set(tmp_key))
                else:
                    if not tmp_ in meth_obj[tuple(sorted(list(tmp_key)))]:
                        meth_obj[tuple(sorted(list(tmp_key)))].append(tmp_)

    return meth_obj

def cds_representation(sequence, annotations, codon_usage):

    cds_reprs = {}

    for annotation in annotations:
        tmp = []
        indeces = set(list(range(len(sequence))))
        internal_representation = {"Internal_index": []}
        i = int(annotation[0])
        j = int(annotation[1])
        tt = irange(i, j, l=len(sequence))
        indeces = set.difference(indeces, set(tt))
        tmp += sum([[el, el, el]
                    for el in irange(i, j, by=3, l=len(sequence))], [])
        
        internal_representation = dict((el, {'original': sequence[el] + sequence[(el+1)%len(sequence)] + sequence[(el+2)%len(sequence)],
                                             'translation': translate(reverse_complement(sequence[el] + sequence[(el+1)%len(sequence)] + sequence[(el+2)%len(sequence)]), codon_usage),
                                             'codon_usage': codon_usage[reverse_complement(sequence[el] + sequence[(el+1)%len(sequence)] + sequence[(el+2)%len(sequence)])]["Percentage"]
                                             })
                                       if annotation[-1] == "-"
                                       else
                                       (el, {'original': sequence[el] + sequence[(el+1)%len(sequence)] + sequence[(el+2)%len(sequence)],
                                             'translation': translate(sequence[el] + sequence[(el+1)%len(sequence)] + sequence[(el+2)%len(sequence)], codon_usage),
                                             'codon_usage': codon_usage[sequence[el] + sequence[(el+1)%len(sequence)] + sequence[(el+2)%len(sequence)]]["Percentage"]}) for el in irange(i, j, by=3, l=len(sequence)))

        tmp += list(indeces)
        tmp.sort()
        internal_representation["Internal_index"] = tmp

        internal_representation.update(
            (index, {'original': sequence[index], 'codon_usage': '-', 'translation': ' '}) for index in list(indeces))

        cds_reprs[tuple(annotation)] = internal_representation

    return cds_reprs

def single_slice_mapping(current_slice, slice_range, current_input_genetic_tool_sequence, current_input_genetic_tool_annotations, codon_usage_table, commercial_mtases, target_motifs, input_rebase_output, candidate_alternative=False):

    # Non-Synonymous slice [0 if the slice is synonymous, 1 if the slice is non-synonymous]
    # Number of target motifs found in the slice
    # Number of modifications introduced in the slice

    synonymous_counter = []
    number_of_targets = 0
    number_of_modifications = 0
    codon_bias_ranking = 0.
    cds_counter = 0

    tmp_pl = ''.join([current_input_genetic_tool_sequence[i] if i not in slice_range else current_slice[slice_range.index(
        i)] for i in range(len(current_input_genetic_tool_sequence))])

    in_cds = False

    cds_codon_usage = {}

    tmp = findpatterns(current_slice,
                       target_motifs, circular=False)

    for patt in tmp.keys():
        if tmp[patt] != []:
            for el in tmp[patt]:
                number_of_targets += 1

    out1 = meth_tr_mapping(
        input_rebase_output[0], current_slice, slice_range)
    out2 = meth_tr_mapping(
        commercial_mtases, current_slice, slice_range)

    mtases = []
    if (len(out1) > 0) & (len(out2) > 0):
        intersection_ = set(out1.keys()) & set(out2.keys())

        n_removed = 0

        if len(intersection_) > 0:
            tars = []
            for el in intersection_:
                tars += [val[0] for val in out1[el]]
            n_tars = list(set(tars))

            for tar in tars:
                if not tar == str(Seq(tar).reverse_complement()):
                    n_removed += 1
                else:
                    n_removed += 2

        number_of_targets = max(
            [0, number_of_targets - n_removed])

        for el in intersection_:
            mtases.append(
                'RM_target_motif:' + repr(out1[el]) + '-->CommercialMTase:' + repr(out2[el]))

    tmp_generated_slice = current_slice

    if not candidate_alternative:
        tmp_generated_slice = tmp_generated_slice.lower()
    else:
        tmp_generated_slice = tmp_generated_slice.upper()

    tmp_found_motifs_slice = findpatterns(
        current_slice, target_motifs, circular=False)
    representation = ''
    if not candidate_alternative:
        first_name = 'SelectTargetPositions'
    else:
        first_name = 'Original'
    representation += ''.join([f'>{first_name}|Slice:(' + repr(slice_range[0]) + ',' + repr(slice_range[-1]+1) + ')']) + \
        "|Number_of_mapped_targets:" + str(number_of_targets) + \
        ''.join(['|Normalized_codonBias_Rank:' + str(round(codon_bias_ranking, 2)) if (cds_counter > 0) else '']) + \
        ''.join(['|' + ''.join(mtases)
                 if ''.join(mtases) != '' else '']) + \
        '|Set_of_Types_of_RMSystems:' + \
        repr(input_rebase_output[-1]) + '\n'

    for tmp_found in tmp_found_motifs_slice.keys():
        if tmp_found_motifs_slice[tmp_found] != []:
            for el in tmp_found_motifs_slice[tmp_found]:
                representation += r"\\\\\\|" + " " * \
                    el[1] + str(tmp_found) + " " + str(el[3]) + "\n"
    representation += "|-5' - " + \
        tmp_generated_slice + " -3'- (+ strand)\n"

    if sum(synonymous_counter) <= 0:
        if (cds_counter > 0):
            cc_representation = cds_representation(
                tmp_pl, current_input_genetic_tool_annotations, codon_usage_table)
            for annotation in current_input_genetic_tool_annotations:
                if set.intersection(set(irange(slice_range[0], slice_range[-1]+1, l=len(current_input_genetic_tool_sequence))), set(irange(int(annotation[0]), int(annotation[1]), l=len(current_input_genetic_tool_sequence)))):
                    codons = []
                    codon_usage = []
                    amino_acids = []
                    representation += r"\--|CDS:" + str(annotation) + "\n"
                    representation += r"\--> "
                    for i in irange(slice_range[0], slice_range[-1]+1, l=len(current_input_genetic_tool_sequence)):
                        if i == slice_range[0]:
                            prev = cc_representation[tuple(
                                annotation)]["Internal_index"][i]

                            codons.append(str(
                                cc_representation[tuple(annotation)][cc_representation[tuple(annotation)]["Internal_index"][i]]["original"]))
                            codon_usage.append(cc_representation[tuple(annotation)][cc_representation[tuple(
                                annotation)]["Internal_index"][i]]["codon_usage"])
                            amino_acids.append(str(
                                cc_representation[tuple(annotation)][cc_representation[tuple(annotation)]["Internal_index"][i]]["translation"]))

                        else:
                            if cc_representation[tuple(annotation)]["Internal_index"][i] != prev:

                                codons.append(str(
                                    cc_representation[tuple(annotation)][cc_representation[tuple(annotation)]["Internal_index"][i]]["original"]))
                                codon_usage.append(cc_representation[tuple(annotation)][cc_representation[tuple(
                                    annotation)]["Internal_index"][i]]["codon_usage"])
                                amino_acids.append(str(
                                    cc_representation[tuple(annotation)][cc_representation[tuple(annotation)]["Internal_index"][i]]["translation"]))
                                prev = cc_representation[tuple(
                                    annotation)]["Internal_index"][i]
                    representation += "   ".join(codons)
                    representation += "\n"
                    representation += r"\+->"
                    representation += ''.join([codon_usage[i] + ' ' if (len(codon_usage[i]) > 1) else codon_usage[i] + '   ' if (len(codon_usage[i+1]) == 1) else codon_usage[i] + '   ' for i in range(len(codon_usage)-1)]) + codon_usage[-1]
                    representation += "\n"
                    representation += r"\*->  "
                    representation += "     ".join(amino_acids)
                    representation += "\n"

    return representation

def tmp_eval_function(inp):
    
    
    or_slice, candidate_alternative, slice_range, current_input_genetic_tool_sequence, current_input_genetic_tool_annotations, codon_usage_table, commercial_mtases, target_motifs, input_rebase_output = inp
    

    synonymous_counter = []
    number_of_targets = 0
    number_of_modifications = 0
    codon_bias_ranking = 0.
    cds_counter = 0

    tmp_pl = ''.join([current_input_genetic_tool_sequence[i] if i not in slice_range else candidate_alternative[slice_range.index(
        i)] for i in range(len(current_input_genetic_tool_sequence))])

    in_cds = False

    cds_codon_usage = {}

    if (set(slice_range) & set(sum([irange(el[0], el[1], l=len(current_input_genetic_tool_sequence))
                                    for el in current_input_genetic_tool_annotations], []))):
        for annotation in current_input_genetic_tool_annotations:
            cds_codon_usage[repr(annotation)] = {}

            tmp = irange(annotation[0], annotation[1],
                         l=len(current_input_genetic_tool_sequence))

            if set(slice_range) & set(irange(annotation[0], annotation[1], l=len(current_input_genetic_tool_sequence))):
                cds_counter += 1
                if annotation[-1] == "-":
                    tmp_candidate_alternative_seq = reverse_complement(
                        extract_subsequence(tmp_pl, tmp))
                    tmp_original_seq = reverse_complement(
                        extract_subsequence(current_input_genetic_tool_sequence, tmp))
                    if (translate(tmp_original_seq, codon_usage_table) == translate(tmp_candidate_alternative_seq, codon_usage_table)):
                        # The slice is synonymous
                        synonymous_counter.append(0)
                        tmp_codons_ranks = sum([codon_usage_table[x]['Inverse_ranking'] for x in re.split(
                            r'(\w{3})', str(tmp_candidate_alternative_seq)) if x != ''])
                        codon_bias_ranking += tmp_codons_ranks
                    else:
                        # The slice is non-synonymous
                        synonymous_counter.append(1)
                else:
                    tmp_candidate_alternative_seq = extract_subsequence(
                        tmp_pl, tmp)
                    tmp_original_seq = extract_subsequence(
                        current_input_genetic_tool_sequence, tmp)
                    if (translate(tmp_original_seq, codon_usage_table) == translate(tmp_candidate_alternative_seq, codon_usage_table)):
                        # The slice is synonymous
                        synonymous_counter.append(0)
                        tmp_codons_ranks = sum([codon_usage_table[x]['Inverse_ranking'] for x in re.split(
                            r'(\w{3})', str(tmp_candidate_alternative_seq)) if x != ''])
                        codon_bias_ranking += tmp_codons_ranks
                    else:
                        # The slice is non-synonymous
                        synonymous_counter.append(1)

    else:
        synonymous_counter = [0]

    tmp = findpatterns(candidate_alternative,
                       target_motifs, circular=False)
    for patt in tmp.keys():
        if tmp[patt] != []:
            for el in tmp[patt]:
                number_of_targets += 1

    out1 = meth_tr_mapping(
        input_rebase_output[0], candidate_alternative, slice_range)
    out2 = meth_tr_mapping(
        commercial_mtases, candidate_alternative, slice_range)

    mtases = []
    if (len(out1) > 0) & (len(out2) > 0):
        intersection_ = set(out1.keys()) & set(out2.keys())
        n_removed = 0
        if len(intersection_) > 0:
            tars = []
            for el in intersection_:
                tars += [val[0] for val in out1[el]]
            n_tars = list(set(tars))

            for tar in tars:
                if not (tar == str(Seq(tar).reverse_complement())):
                    n_removed += 1
                else:
                    n_removed += 1

        number_of_targets = max(
            [0, number_of_targets - n_removed])

        for el in intersection_:
            mtases.append(
                'RM_target_motif:' + repr(out1[el]) + '-->CommercialMTase:' + repr(out2[el]))

    tmp_generated_slice = ''
    for pos in range(len(candidate_alternative)):
        if candidate_alternative[pos] != or_slice[pos]:
            number_of_modifications += 1
            tmp_generated_slice += candidate_alternative[pos].upper()
        else:
            tmp_generated_slice += candidate_alternative[pos].lower()

    if (candidate_alternative == or_slice):
            synonymous_counter = [-1]
            tmp_generated_slice = tmp_generated_slice.upper()

    tmp_found_motifs_slice = findpatterns(
        candidate_alternative, target_motifs, circular=False)
    representation = ''
    representation += ''.join(['>SynonymousChange:Yes' if sum(synonymous_counter) == 0 else '>SynonymousChange:No' if sum(synonymous_counter) > 0 else '>Original|Slice:(' + repr(slice_range[0]) + ',' + repr(slice_range[-1]+1) + ')']) + \
        "|Number_of_left_targets:" + str(number_of_targets) + \
        '|Number_of_introduced_modifications:' + \
        str(number_of_modifications) + \
        ''.join(['|Normalized_codonBias_Rank:' + str(round(codon_bias_ranking, 2)) if (cds_counter > 0) else '']) + \
        ''.join(['|' + ''.join(mtases)
                 if ''.join(mtases) != '' else '']) + \
        '|Set_of_Types_of_RMSystems:' + \
        repr(input_rebase_output[-1]) + '\n'
    for tmp_found in tmp_found_motifs_slice.keys():
        if tmp_found_motifs_slice[tmp_found] != []:
            for el in tmp_found_motifs_slice[tmp_found]:
                representation += r"\\\\\\|" + " " * \
                    el[1] + str(tmp_found) + " " + str(el[3]) + "\n"
    representation += "|-5' - " + \
        tmp_generated_slice + " -3'- (+ strand)\n"

    if sum(synonymous_counter) <= 0:
        if (cds_counter > 0):
            cc_representation = cds_representation(
                tmp_pl, current_input_genetic_tool_annotations, codon_usage_table)
            for annotation in current_input_genetic_tool_annotations:
                if set.intersection(set(irange(slice_range[0], slice_range[-1]+1, l=len(current_input_genetic_tool_sequence))), set(irange(int(annotation[0]), int(annotation[1]), l=len(current_input_genetic_tool_sequence)))):
                    codons = []
                    codon_usage = []
                    amino_acids = []
                    representation += r"\--|CDS:" + str(annotation) + "\n"
                    representation += r"\--> "
                    for i in irange(slice_range[0], slice_range[-1]+1, l=len(current_input_genetic_tool_sequence)):
                        if i == slice_range[0]:
                            prev = cc_representation[tuple(
                                annotation)]["Internal_index"][i]

                            codons.append(str(
                                cc_representation[tuple(annotation)][cc_representation[tuple(annotation)]["Internal_index"][i]]["original"]))
                            codon_usage.append(cc_representation[tuple(annotation)][cc_representation[tuple(
                                annotation)]["Internal_index"][i]]["codon_usage"])
                            amino_acids.append(str(
                                cc_representation[tuple(annotation)][cc_representation[tuple(annotation)]["Internal_index"][i]]["translation"]))

                        else:
                            if cc_representation[tuple(annotation)]["Internal_index"][i] != prev:

                                codons.append(str(
                                    cc_representation[tuple(annotation)][cc_representation[tuple(annotation)]["Internal_index"][i]]["original"]))
                                codon_usage.append(cc_representation[tuple(annotation)][cc_representation[tuple(
                                    annotation)]["Internal_index"][i]]["codon_usage"])
                                amino_acids.append(str(
                                    cc_representation[tuple(annotation)][cc_representation[tuple(annotation)]["Internal_index"][i]]["translation"]))
                                prev = cc_representation[tuple(
                                    annotation)]["Internal_index"][i]
                    representation += "   ".join(codons)
                    representation += "\n"
                    representation += r"\+->"
                    representation += ''.join([codon_usage[i] + ' ' if (len(codon_usage[i]) > 1) else codon_usage[i] + '   ' if (
                        len(codon_usage[i+1]) == 1) else codon_usage[i] + '   ' for i in range(len(codon_usage)-1)]) + codon_usage[-1]
                    representation += "\n"
                    representation += r"\*->  "
                    representation += "     ".join(amino_acids)
                    representation += "\n"        
    return (sum(synonymous_counter), number_of_targets, number_of_modifications, codon_bias_ranking, tmp_generated_slice, representation)

def evaluate_slices(or_slice='', candidate_alternatives=[], slice_range=[], current_input_genetic_tool_sequence='', current_input_genetic_tool_annotations=[], codon_usage_table='', commercial_mtases={}, target_motifs=[], input_rebase_output={}, ncores=1):
    
    slices = []
    # Non-Synonymous slice [0 if the slice is synonymous, 1 if the slice is non-synonymous]
    # Number of target motifs found in the slice
    # Number of modifications introduced in the slice
    
    with Pool(ncores) as p:
        slices = list(p.map(
            tmp_eval_function, 
            [
                (
                    or_slice, 
                    candidate_alternative,
                    slice_range,
                    current_input_genetic_tool_sequence,
                    current_input_genetic_tool_annotations,
                    codon_usage_table,
                    commercial_mtases,
                    target_motifs,
                    input_rebase_output
                )
                
                for candidate_alternative in candidate_alternatives]))
    
#     for candidate_alternative in candidate_alternatives:
#         synonymous_counter = []
#         number_of_targets = 0
#         number_of_modifications = 0
#         codon_bias_ranking = 0.
#         cds_counter = 0

#         tmp_pl = ''.join([current_input_genetic_tool_sequence[i] if i not in slice_range else candidate_alternative[slice_range.index(
#             i)] for i in range(len(current_input_genetic_tool_sequence))])

#         in_cds = False

#         cds_codon_usage = {}

#         if (set(slice_range) & set(sum([irange(el[0], el[1], l=len(current_input_genetic_tool_sequence))
#                                         for el in current_input_genetic_tool_annotations], []))):
#             for annotation in current_input_genetic_tool_annotations:
#                 cds_codon_usage[repr(annotation)] = {}

#                 tmp = irange(annotation[0], annotation[1],
#                              l=len(current_input_genetic_tool_sequence))

#                 if set(slice_range) & set(irange(annotation[0], annotation[1], l=len(current_input_genetic_tool_sequence))):
#                     cds_counter += 1
#                     if annotation[-1] == "-":
#                         tmp_candidate_alternative_seq = reverse_complement(
#                             extract_subsequence(tmp_pl, tmp))
#                         tmp_original_seq = reverse_complement(
#                             extract_subsequence(current_input_genetic_tool_sequence, tmp))
#                         if (translate(tmp_original_seq, codon_usage_table) == translate(tmp_candidate_alternative_seq, codon_usage_table)):
#                             # The slice is synonymous
#                             synonymous_counter.append(0)
#                             tmp_codons_ranks = sum([codon_usage_table[x]['Inverse_ranking'] for x in re.split(
#                                 r'(\w{3})', str(tmp_candidate_alternative_seq)) if x != ''])
#                             codon_bias_ranking += tmp_codons_ranks
#                         else:
#                             # The slice is non-synonymous
#                             synonymous_counter.append(1)
#                     else:
#                         tmp_candidate_alternative_seq = extract_subsequence(
#                             tmp_pl, tmp)
#                         tmp_original_seq = extract_subsequence(
#                             current_input_genetic_tool_sequence, tmp)
#                         if (translate(tmp_original_seq, codon_usage_table) == translate(tmp_candidate_alternative_seq, codon_usage_table)):
#                             # The slice is synonymous
#                             synonymous_counter.append(0)
#                             tmp_codons_ranks = sum([codon_usage_table[x]['Inverse_ranking'] for x in re.split(
#                                 r'(\w{3})', str(tmp_candidate_alternative_seq)) if x != ''])
#                             codon_bias_ranking += tmp_codons_ranks
#                         else:
#                             # The slice is non-synonymous
#                             synonymous_counter.append(1)

#         else:
#             synonymous_counter = [0]

#         tmp = findpatterns(candidate_alternative,
#                            target_motifs, circular=False)
#         for patt in tmp.keys():
#             if tmp[patt] != []:
#                 for el in tmp[patt]:
#                     number_of_targets += 1

#         out1 = meth_tr_mapping(
#             input_rebase_output[0], candidate_alternative, slice_range)
#         out2 = meth_tr_mapping(
#             commercial_mtases, candidate_alternative, slice_range)

#         mtases = []
#         if (len(out1) > 0) & (len(out2) > 0):
#             intersection_ = set(out1.keys()) & set(out2.keys())
#             n_removed = 0
#             if len(intersection_) > 0:
#                 tars = []
#                 for el in intersection_:
#                     tars += [val[0] for val in out1[el]]
#                 n_tars = list(set(tars))

#                 for tar in tars:
#                     if not (tar == str(Seq(tar).reverse_complement())):
#                         n_removed += 1
#                     else:
#                         n_removed += 1

#             number_of_targets = max(
#                 [0, number_of_targets - n_removed])

#             for el in intersection_:
#                 mtases.append(
#                     'RM_target_motif:' + repr(out1[el]) + '-->CommercialMTase:' + repr(out2[el]))

#         tmp_generated_slice = ''
#         for pos in range(len(candidate_alternative)):
#             if candidate_alternative[pos] != or_slice[pos]:
#                 number_of_modifications += 1
#                 tmp_generated_slice += candidate_alternative[pos].upper()
#             else:
#                 tmp_generated_slice += candidate_alternative[pos].lower()

#         if (candidate_alternative == or_slice):
#                 synonymous_counter = [-1]
#                 tmp_generated_slice = tmp_generated_slice.upper()

#         tmp_found_motifs_slice = findpatterns(
#             candidate_alternative, target_motifs, circular=False)
#         representation = ''
#         representation += ''.join(['>SynonymousChange:Yes' if sum(synonymous_counter) == 0 else '>SynonymousChange:No' if sum(synonymous_counter) > 0 else '>Original|Slice:(' + repr(slice_range[0]) + ',' + repr(slice_range[-1]+1) + ')']) + \
#             "|Number_of_left_targets:" + str(number_of_targets) + \
#             '|Number_of_introduced_modifications:' + \
#             str(number_of_modifications) + \
#             ''.join(['|Normalized_codonBias_Rank:' + str(round(codon_bias_ranking, 2)) if (cds_counter > 0) else '']) + \
#             ''.join(['|' + ''.join(mtases)
#                      if ''.join(mtases) != '' else '']) + \
#             '|Set_of_Types_of_RMSystems:' + \
#             repr(input_rebase_output[-1]) + '\n'
#         for tmp_found in tmp_found_motifs_slice.keys():
#             if tmp_found_motifs_slice[tmp_found] != []:
#                 for el in tmp_found_motifs_slice[tmp_found]:
#                     representation += r"\\\\\\|" + " " * \
#                         el[1] + str(tmp_found) + " " + str(el[3]) + "\n"
#         representation += "|-5' - " + \
#             tmp_generated_slice + " -3'- (+ strand)\n"

#         if sum(synonymous_counter) <= 0:
#             if (cds_counter > 0):
#                 cc_representation = cds_representation(
#                     tmp_pl, current_input_genetic_tool_annotations, codon_usage_table)
#                 for annotation in current_input_genetic_tool_annotations:
#                     if set.intersection(set(irange(slice_range[0], slice_range[-1]+1, l=len(current_input_genetic_tool_sequence))), set(irange(int(annotation[0]), int(annotation[1]), l=len(current_input_genetic_tool_sequence)))):
#                         codons = []
#                         codon_usage = []
#                         amino_acids = []
#                         representation += r"\--|CDS:" + str(annotation) + "\n"
#                         representation += r"\--> "
#                         for i in irange(slice_range[0], slice_range[-1]+1, l=len(current_input_genetic_tool_sequence)):
#                             if i == slice_range[0]:
#                                 prev = cc_representation[tuple(
#                                     annotation)]["Internal_index"][i]

#                                 codons.append(str(
#                                     cc_representation[tuple(annotation)][cc_representation[tuple(annotation)]["Internal_index"][i]]["original"]))
#                                 codon_usage.append(cc_representation[tuple(annotation)][cc_representation[tuple(
#                                     annotation)]["Internal_index"][i]]["codon_usage"])
#                                 amino_acids.append(str(
#                                     cc_representation[tuple(annotation)][cc_representation[tuple(annotation)]["Internal_index"][i]]["translation"]))

#                             else:
#                                 if cc_representation[tuple(annotation)]["Internal_index"][i] != prev:

#                                     codons.append(str(
#                                         cc_representation[tuple(annotation)][cc_representation[tuple(annotation)]["Internal_index"][i]]["original"]))
#                                     codon_usage.append(cc_representation[tuple(annotation)][cc_representation[tuple(
#                                         annotation)]["Internal_index"][i]]["codon_usage"])
#                                     amino_acids.append(str(
#                                         cc_representation[tuple(annotation)][cc_representation[tuple(annotation)]["Internal_index"][i]]["translation"]))
#                                     prev = cc_representation[tuple(
#                                         annotation)]["Internal_index"][i]
#                         representation += "   ".join(codons)
#                         representation += "\n"
#                         representation += r"\+->"
#                         representation += ''.join([codon_usage[i] + ' ' if (len(codon_usage[i]) > 1) else codon_usage[i] + '   ' if (
#                             len(codon_usage[i+1]) == 1) else codon_usage[i] + '   ' for i in range(len(codon_usage)-1)]) + codon_usage[-1]
#                         representation += "\n"
#                         representation += r"\*->  "
#                         representation += "     ".join(amino_acids)
#                         representation += "\n"
#         slices.append(
#                 (sum(synonymous_counter), number_of_targets, number_of_modifications, codon_bias_ranking, tmp_generated_slice, representation))
    slices = list(set(slices))
    slices.sort()

    return slices

# Read inputs

def read_input_seq(input_file):
    '''
    The input_file must be in gbk format!

    This function returns a tuple of elements:
        (correct, annotations, total_sequence, lines)
    
    correct: True/False --> indicates if the input genetic sequence is correctly annotated and all the bases are canonical
    annotations: a dictionary for each sequence in gbk with all the CDS and the corresponding sequence. 
                CDS are indicated as a list of three elements [i, j, '-/+'],
                where, 
                    i: start position of the cds (included)
                    j: ending position of the cds (not included)
                    '-/+': indicates the orientation of the CDS
                        '-': indicates 3'-5' orientation --> to translate reverse complement is required
                        '+': indicates 5'-3' orientation
    
    total_sequence: is the concatenation of all the sequences in the .gbk file. 
                The concatenation is performed only if the annotations are correct,
                otherwise an empty string is returned
    
    lines: is the list of lines that compose the .gbk file. 
            This element is useful for retrieving the output syngenic sequence in .gbk format (not only fasta)
            without relying on Biopython --> change only the sequence

    The annotations are re-encoded according to c/python language convention:
    - numbering starts with 0
    - last position of a sequence is the actual length of the sequence reduced by one
    - ranges are made with first position included, second excluded

    '''

    out_config_file = {
        'correct': True,
        'messages': [],
        'output': [],
        'total_cds_concats': '',
        'annotations': {}
    }
    

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        
        try:
            sequences = list(SeqIO.parse(input_file, "genbank"))
        except ValueError as e:
            sys.exit(e)

        if len(w) > 0:
            out_config_file['correct'] = out_config_file['correct'] & False
            out_config_file['messages'].append('\n'.join([str(el.message) for el in w]))
            return out_config_file

    if not len(sequences) >= 1:
        out_config_file['correct'] = out_config_file['correct'] & False
        out_config_file['messages'].append('No input sequence found.')
        return out_config_file

    j = 0

    for seq in sequences:
        cds = [el for el in seq.features if el.type == 'CDS']
        out_config_file['annotations'][j] = {}
        out_config_file['annotations'][j]['cds'] = []
        out_config_file['annotations'][j]['cds_sequences'] = []
        out_config_file['annotations'][j]['sequence'] = str(seq.seq)
        
        for feat in cds:
            if feat.location == None:
                out_config_file['correct'] = False
            else:
                range_feat = list(feat.location)

                if feat.strand == 1:
                    start = range_feat[0]
                    end = range_feat[-1] + 1
                else:
                    start = range_feat[-1]
                    end = range_feat[0] + 1
                strand = '-' if feat.strand == -1 else '+'
                
                out_config_file['annotations'][j]['cds'].append([start, end, strand])
                
                tmp_str = str(feat.extract(seq).seq)

                out_config_file['annotations'][j]['cds_sequences'].append(tmp_str)
                if not len(tmp_str) % 3 == 0:
                    out_config_file['correct'] = False
                    out_config_file['messages'].append(f'Error in annotation {feat}: len(sequence) not a multiple of three: {tmp_str}!')
                if not ''.join(sorted(list(set(tmp_str.upper())))) == 'ACGT':
                    out_config_file['correct'] = False
                    out_config_file['messages'].append(f'Error in annotation {feat}: non-canonical bases found!')

                out_config_file['total_cds_concats'] += tmp_str
        j += 1

    return out_config_file

def read_commercial(input_commercial):
    '''
    This function reads the file CAM_list_formatted.xlsx, which contains the commercial MTases.
    the output is a dictionary in which for each MTase:

    -meth_base: position (numbering starts from 0) of the methylated base in 5'-3' orientation
    -comp_meth_base: position (numbering starts from 0) of the methylated base in 3'-5' orientation
    -meth_type: type of methylation on the base indicated as meth_base in the 5'-3' orientation
    -comp_meth_type: type of methylation on the base indicated as meth_base in the 3'-5' orientation

    **Methylation type does not correspond to the type of the RM system**
    **Methylation type indicates the position on the G/T molecule on which the modification is performed**

    '''

    def standardize_commercial(meth_commercial):
        n_dict = meth_commercial
        for key1 in n_dict.keys():
            for key2 in {'meth_base', 'comp_meth_base'}:
                n_dict[key1][key2] -= 1

        return n_dict

    tmp_lines = [el.split('\t') for el in input_commercial.split('\n')]
    header = tmp_lines[0]
    tmp_lines = tmp_lines[1:]

    d = {}

    for i in range(len(header)):
        d[header[i]] = []
        for el in tmp_lines:

            if i in [2, 3, 4, 5]:
                d[header[i]].append(int(el[i]))
            else:
                d[header[i]].append(el[i])

    d = pd.DataFrame(d)

    d = d.drop(["Methyl-Modification introduced (5'-3')",
                'availability'], axis=1)
    meth_commercial = d.T.to_dict(orient='dict')

    meth_commercial = standardize_commercial(meth_commercial)

    return meth_commercial

def read_synpl(file_, l, mapping=False):

    out_config_file = {
        'correct': True,
        'messages': [],
        'output': []
    }

    if not mapping:
        curr_header = synpl_header
    else:
        curr_header = synpl_map_header


    i = -1
    out_config_file['output'] = {}
    with open(file_, 'r') as handle:

        tmp_lines = handle.readlines()
        joined_lines = ''.join(tmp_lines)

        header = joined_lines[:len(curr_header)]

        if not header == curr_header:
            out_config_file['messages'].append(f'Invalid header during parsing {file_} .synpl file!')
            out_config_file['correct'] = out_config_file['correct'] & False
            return out_config_file

        oth_lines = joined_lines[len(curr_header):]

        num_separators = sum([1 if el.strip() == '_'*69 else 0 for el in tmp_lines])/2
        num_sequences = sum([1 if el.startswith('>') else 0 for el in tmp_lines])

        if not (num_separators == num_sequences):
            out_config_file['messages'].append(f'Invalid number of separators during parsing {file_} .synpl file!')
            out_config_file['correct'] = out_config_file['correct'] & False
            return out_config_file

        for line in tmp_lines:
            if ((line.startswith("+")) | (line.startswith("_"))):
                continue
            elif (line.startswith(">")):
                i += 1
                out_config_file['output'][i] = [line]
            elif ((line.startswith("|")) | (line.startswith("\\"))):
                out_config_file['output'][i].append(line)

        tmp_name = file_.split('/')[-1].split('.')[0]

        if not 'slice' in tmp_name:
            out_config_file['messages'].append(f'Incorrect .synpl file name {file_}.')
            out_config_file['correct'] = out_config_file['correct'] & False
            return out_config_file

        start, end = [int(el) for el in re.sub('\D', ' ', tmp_name).split(' ') if el != '']
        
        if not f'>Original|Slice:({start},{end})' in out_config_file['output'][0][0]:
            out_config_file['messages'].append(f'Incorrect range of the slice!')
            out_config_file['correct'] = out_config_file['correct'] & False
            return out_config_file

        for key in out_config_file['output'].keys():

            tmp_match = [el for el in out_config_file['output'][0] if re.match(r"^\|\-5' \- [ATCGatcg]+ \-3'\- \(\+ strand\)$", el.strip())]

            if not len(tmp_match) == 1:
                out_config_file['messages'].append(f'No slice sequence found while partsing {file_}')
                out_config_file['correct'] = out_config_file['correct'] & False
                return out_config_file
            
            tmp_search = re.search(r"[ATCGatcg]+", tmp_match[0].strip())

            if not (len(tmp_search.group()) == len(irange(start, end, l=l))):
                out_config_file['messages'].append(f'Found incorrect slice length while parsing {file_}')
                out_config_file['correct'] = out_config_file['correct'] & False
                return out_config_file
            
    return out_config_file

# SyToGen Functionalities Utils

def sequence_preprocess(input_sequence, backbone):

    out_config_file = {
        'correct': True,
        'messages': [],
        'output': []
    }

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        
        try:
            record = list(SeqIO.parse(input_sequence, "genbank"))
        except ValueError as e:
            sys.exit(e)

        if len(w) > 0:
            out_config_file['correct'] = out_config_file['correct'] & False
            out_config_file['messages'].append('\n'.join([str(el.message) for el in w]))
            return out_config_file

        if not len(record) == 1:
            out_config_file['correct'] = out_config_file['correct'] & False
            out_config_file['messages'].append('Number of input sequence != 1.')
            return out_config_file
        else:
            record = record[0]

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")

        try:
            backbone = list(SeqIO.parse(backbone, "fasta"))
        except ValueError as e:
            out_config_file['correct'] = out_config_file['correct'] & False
            out_config_file['messages'].append('\n'.join([str(el.message) for el in e]))
            return out_config_file
        
        if len(w) > 0:
            out_config_file['correct'] = out_config_file['correct'] & False
            out_config_file['messages'].append('\n'.join([str(el.message) for el in w]))
            return out_config_file

        if not len(backbone) == 1:
            out_config_file['correct'] = out_config_file['correct'] & False
            out_config_file['messages'].append('Number of backbone sequence != 1.')
            return out_config_file
        else:
            backbone = str(backbone[0].seq).upper()

    tmp = str(record.seq + record.seq[:len(backbone)-1]).upper()

    if not (len(set(tmp) - set(['A', 'C', 'T', 'G'])) == 0):
        out_config_file['correct'] = out_config_file['correct'] & False
        out_config_file['messages'].append(f'Input sequence contains non-canonical bases: {set(tmp) - set(["A", "C", "T", "G"])}')
        return out_config_file
    
    if not (len(set(backbone) - set(['A', 'C', 'T', 'G'])) == 0):
        out_config_file['correct'] = out_config_file['correct'] & False
        out_config_file['messages'].append(f'Input backbone sequence contains non-canonical bases: {set(backbone) - set(["A", "C", "T", "G"])}')
        return out_config_file

    if tmp.count(backbone) > 1:
        out_config_file['correct'] = out_config_file['correct'] & False
        out_config_file['messages'].append(f'Backbone found more than 1 time in the input genetic tool!')
        return out_config_file
    elif tmp.count(backbone) == 0:
        out_config_file['correct'] = out_config_file['correct'] & False
        out_config_file['messages'].append(f'Backbone not found in the input genetic tool!')
        return out_config_file
    else:
        start = tmp.find(backbone.upper()) % len(record.seq)
        end = (tmp.find(backbone.upper()) + len(backbone)) % len(record.seq)

        if start < end:
            out_record = record[:start] + record[end:]
        else:
            out_record = record[end:start]
        
        out_record.annotations["molecule_type"] = "DNA"
        out_config_file['output'] = out_record

    return out_config_file

# SyToGen Functionalities

# Functionality

def gibson_assembly(input_sequence):  

    out_config_file = {
        'correct': True,
        'messages': [],
        'representation': '',
        'primers': '',
        'assembled_sequence': '',
        'out_gibson': ''
    }

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        
        try:
            record = list(SeqIO.parse(input_sequence, "fasta"))
        except ValueError as e:
            sys.exit(e)

        if len(w) > 0:
            out_config_file['correct'] = out_config_file['correct'] & False
            out_config_file['messages'].append('\n'.join([str(el.message) for el in w]))
            return out_config_file

        if not len(record) > 1:
            out_config_file['correct'] = out_config_file['correct'] & False
            out_config_file['messages'].append('Number of input sequence not > 1.')
            return out_config_file


    tmp_fa = assembly_fragments([primer_design(Dseqrecord(str(el.seq))) for el in record] + [primer_design(Dseqrecord(str(record[0].seq)))])

    fa = pcr(tmp_fa[-1].forward_primer, tmp_fa[0].reverse_primer, primer_design(Dseqrecord(str(record[0].seq))))

    nn_l = [fa] + tmp_fa[1:-1]

    assemblyobj = Assembly(nn_l)

    nn = assemblyobj.assemble_circular()

    primers = [el.strip() for el in nn[0].detailed_figure().split('\n') if el.strip().isupper()]

    out_file = nn[0].detailed_figure()

    out_config_file['primers'] = primers
    out_config_file['representation'] = out_file
    
    out_config_file['out_gibson'] = 'Primers for assembly:\n\t' + '\n\t'.join(primers) + '\n\nAssembly representation:\n' + out_file

    features = []
    tmp_out_seq = Seq('')

    for tmp_record in record:
        if 'CDS' in tmp_record.id:
            if '-' in tmp_record.id:
                features.append(SeqFeature(FeatureLocation(len(tmp_out_seq), len(tmp_out_seq)+len(tmp_record.seq)), type="CDS", strand=-1))
            else:
                features.append(SeqFeature(FeatureLocation(len(tmp_out_seq), len(tmp_out_seq)+len(tmp_record.seq)), type="CDS", strand=+1))
        
        tmp_out_seq = tmp_out_seq + tmp_record.seq
    
    out_seq = SeqRecord(tmp_out_seq, features = features, name='')

    out_config_file['assembled_sequence'] = out_seq

    return out_config_file

def make_assembly(args):

    out_config_file = gibson_assembly(args.input_sequence)
        
    if not out_config_file['correct']:
        error('\n'.join(out_config_file['messages']))
    else:
        with open(os.path.join(args.output_folder, 'multifasta_example_out.fna'), 'w') as handle:
            handle.write(out_config_file['out_gibson'])

    return

# Functionality

def make_preprocess(args):

    out_config_file = sequence_preprocess(args.input_sequence, args.input_backbone)

    if not out_config_file['correct']:
        error('\n'.join(out_config_file['messages']))
    else:
        tmp_name = os.path.basename(args.input_sequence).split('.')[0]
        SeqIO.write(out_config_file['output'], os.path.join(args.output_folder, 'input_sequence.txt'), "genbank")

    return

# Functionality

def candidate_builder(args):

    synpl_files = [el for el in glob.glob(os.path.join(args.synpl_folder, "*.synpl")) if not 'temporary_mapping' in el]

    if not len(synpl_files) > 0:
        error(f'No input .synpl files in the provided directory {args.synpl_folder}.')

    if not os.path.exists(args.output_folder):
        error(f'Indicated output folder {args.output_folder} does not exist!')

    if not os.path.isdir(args.output_folder):
        error(f'Indicated output folder {args.output_folder} is not a directory!')

    tmp_auto_values = ['rank', 'auto', 'star']
    tmp_name = args.input_sequence.split('/')[-1].split('.')[0]
    log_file = os.path.join(args.output_folder, f'sytogen_candidate_syngeneic__{tmp_auto_values[args.auto]}__.log')

    if os.path.exists(log_file):
        os.remove(log_file)

    logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', filename=log_file, level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')

    # Read the original plasmid
    # Read the different input configuration files
    # Introduce the modifications
    # Produce the statistics for the report

    logging.info("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    logging.info("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    logging.info("Starting candidate syngenic assembler.")
    logging.info("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    logging.info("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

    if args.verbose:
        info("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        info("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        info("Starting candidate syngenic assembler.")
        info("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        info("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

    # Read the input plasmid

    input_plasmid = read_input_seq(args.input_sequence)

    if not input_plasmid['correct']:
        error('\n'.join(input_plasmid['messages']))
    else:
        if len(input_plasmid['annotations'].keys()) == 1:
            if len(input_plasmid['annotations'][0]['sequence']) > 0:
                plasmid_sequence = input_plasmid['annotations'][0]['sequence']
                plasmid_annotations = input_plasmid['annotations'][0]['cds']
                l = len(plasmid_sequence)
            else:
                error('More sequences in the input plasmid file OR the sequence of the input file is empty!')
        else:
            error('Incorrect number of inputs for the plasmid sequence!')

    # Read the target motifs

    input_rebase_output = read_rebase_output(args.input_rm_systems)
    input_patterns = [input_rebase_output[0][el]['rec_seq']
                        for el in input_rebase_output[0].keys()]
    max_len = max([len(inp_pat) for inp_pat in input_patterns])

    if 4 in input_rebase_output[1]:
        TYPE4 = True
    else:
        TYPE4 = False

    logging.info("")
    logging.info("Type IV RM System status: " + str(TYPE4))

    if args.verbose:
        info("")
        info("Type IV RM System status: " + str(TYPE4))

    # Present RM-systems -- commercially available MTases
    mtases = input_rebase_output[-1]

    # Circular representation
    nn_seq = {}
    for i in range(len(plasmid_sequence)):
        nn_seq[i] = plasmid_sequence[i]

    modified_positions = []

    # Operate target base modification
    for file_ in synpl_files:

        logging.info("")
        logging.info("Considering slice: " + file_.split("/")[-1] + '; ' + str(synpl_files.index(file_)) + '/' + str(len(synpl_files)))

        if args.verbose:
            info("")
            info("Considering slice: " + file_.split("/")
                    [-1] + '; ' + str(synpl_files.index(file_)) + '/' + str(len(synpl_files)))

        outs = read_synpl(file_, len(plasmid_sequence))#; print(outs); raise Exception("Here")
        print("Here")
        
        if not outs['correct']:
            # error('\n'.join(outs['messages']))
            print('\n'.join(outs['messages']))
            continue
        else:
            outs = outs['output']

        if len(outs) > 1:
            original = outs[0]; print(args.auto)

            if args.auto == 0:

                config = configparser.ConfigParser()
                config.read(args.input_config_file)
                if not 'CANDIDATE_ALTERNATIVES' in config.sections():
                        error(f'Missing SYNPLEST in {args.input_config_file} sections.')

                if file_.split("/")[-1] in config['CANDIDATE_ALTERNATIVES'].keys():
                    outs = outs[int(config['CANDIDATE_ALTERNATIVES'][file_.split("/")[-1]])]
                else:
                    n_left_targets_slice_original = int(outs[0][0].split("Number_of_left_targets:")[1].split("|")[0])
                    n_left_targets_slice1 = int(outs[1][0].split("Number_of_left_targets:")[1].split("|")[0])

                    print(file_, n_left_targets_slice_original, n_left_targets_slice1)

                    if n_left_targets_slice1 < n_left_targets_slice_original:
                        outs = outs[1]
                    else:
                        outs = original
            elif args.auto == 1:
                print("HEREEE")
                n_left_targets_slice_original = int(outs[0][0].split("Number_of_left_targets:")[1].split("|")[0])
                n_left_targets_slice1 = int(outs[1][0].split("Number_of_left_targets:")[1].split("|")[0])
                print(outs)
                print(file_, n_left_targets_slice_original, n_left_targets_slice1)

                if n_left_targets_slice1 < n_left_targets_slice_original:
                    outs = outs[1]
                else:
                    outs = original
            elif args.auto == 2:
                count_proposed_changes = 0

                for key in outs.keys():
                    if outs[key][0].startswith('>*'):
                        count_proposed_changes += 1
                        outs = outs[key]
                        break

                if count_proposed_changes == 0:
                    outs = outs[0]
                    
            mtases = set.union(mtases, set(re.findall('[A-Z]+-[0-9]+-[0-9]+-[0-9]+-[0-9]+', outs[0])))
            

            logging.info("")
            logging.info("The considered file at the present iteration is: " + os.path.basename(file_))
            
            if args.verbose:
                info("")
                info("The considered file at the present iteration is: " + os.path.basename(file_))


            new_seq = [el.split("-")[2].strip() for el in outs if el.startswith("|")][0]
            
            if not len(new_seq) < len(plasmid_sequence):
                error('Current slice is too long to be considered as correct!')

            if not len(new_seq) > 2*max_len:
                error(f'Incorrect lenght of the current slice {file_.split("/")[-1]}.')
            
            ran = [re.findall(r'\d+', el) for el in original if el.startswith(">")][0]

            rr = irange(int(ran[0]), int(ran[1]), l=l)
            for i in range(len(rr)):
                if nn_seq[rr[i]].upper() != new_seq.upper()[i].upper():
                    nn_seq[rr[i]] = new_seq.upper()[i].upper()
                    modified_positions.append(rr[i])

    out_primer_files = ''

    # for mod_pos in sorted(modified_positions):

    #     if args.verbose:
    #         info(f'\nConsidering modification at position {mod_pos}:')
    #     try:
    #         curr_seq = Dseqrecord(''.join([nn_seq[i].lower() if i != mod_pos else nn_seq[i].upper() for i in irange(mod_pos-5, mod_pos+24, l=len(plasmid_sequence))]))
    #         curr_primer = primer_design(curr_seq)
    #         curr_primer_forward = str(curr_primer.forward_primer.seq).upper()
    #         out_primer_files += f'--------------------------------------------------------------------\nConsidering modification at position {mod_pos}:\n'
    #         out_primer_files += str(curr_primer.figure())
    #         out_primer_files += '\n'
    #         out_primer_files += str(dbd_program(curr_primer))
    #         out_primer_files += '\n'
    #         out_primer_files += 'Primer: ' + str(curr_primer.forward_primer.seq)
    #         out_primer_files += '\nTm_Wallace: ' + '%0.2f' % mt.Tm_Wallace(curr_primer_forward) + '\nTm_Wallace: ' + '%0.2f' % mt.Tm_Wallace(curr_primer_forward) + '\nTm_GC: ' + '%0.2f' % mt.Tm_GC(curr_primer_forward) + '' + '%0.2f' % mt.Tm_NN(curr_primer_forward) + '\n'
    #     except ValueError as e:
    #         info(e)
    #         info(f'Error when considering position {mod_pos}.')
    #     except TypeError as e:
    #         info(e)
    #         info(f'Error when considering position {mod_pos}.')

    #     if args.verbose:
    #         info(out_primer_files)

    # if out_primer_files != '':
    #     with open(os.path.join(args.output_folder,f'primers__{tmp_auto_values[args.auto]}__.txt'), 'w') as handle:
    #         handle.write(out_primer_files)

    out_def = ''.join([nn_seq[i] for i in sorted(list(nn_seq.keys()))])

    tmp_out_gbk_sequence = [' '.join([el for el in re.split(r'(\w{10})', out_def.lower()) if el != ''][i:i+6]) for i in range(0,
                                                                                                                                len([el for el in re.split(r'(\w{10})', out_def.lower()) if el != '']), 6)]
    tmp_out_seq = '\n'.join([' '*(9-len(str(i*60+1))) + str(i*60+1) + ' ' +
                                tmp_out_gbk_sequence[i] for i in range(len(tmp_out_gbk_sequence))]) + '\n//'

    with open(args.input_sequence, 'r') as input_handle:
        tmp_lines = ''.join(input_handle.readlines())

    new_out = tmp_lines.split('ORIGIN')[0] + 'ORIGIN\n' + tmp_out_seq

    logging.info("")
    logging.info("The length of the original plasmid/minicircle is: " + str(len(plasmid_sequence)))
    logging.info("The length of the proposed syngenic plasmid/minicircle is: " + str(len(out_def)))

    if (len(out_def) != len(plasmid_sequence)):

        logging.info("")
        logging.info("Error!! Different lengths between the original and the new sequence!!")

        if args.verbose:
            info("")
            info("Error!! Different lengths between the original and the new sequence!!")


    with open(os.path.join(args.output_folder,f'candidate_syngenic_sequence__{tmp_auto_values[args.auto]}__.fa'), 'w') as handle:
        handle.write('>Annotations:' + str(plasmid_annotations) + '|Mtases:' + str([el for el in mtases]) + '\n' + out_def)
    with open(os.path.join(args.output_folder, f'candidate_syngenic_sequence__{tmp_auto_values[args.auto]}__.gbk'), 'w') as handle:
        handle.write(new_out)

    f_patterns = findpatterns(plasmid_sequence, input_patterns, circular = False)
    ff_patterns = findpatterns(out_def, input_patterns, circular = False)

    logging.info("Number of modifications: " + str(len(modified_positions)))

    n_remained_target_motifs = 0
    remained_target_motifs = []

    for key in ff_patterns.keys():
        if len(ff_patterns[key]) > 0:
            n_remained_target_motifs += 1
            remained_target_motifs.append((key, ff_patterns[key]))

    other_text = '<p>Number of not-removed target motifs: ' + str(n_remained_target_motifs)+ '</p>' + \
    '<p>Not-removed target motifs are: ' + str(remained_target_motifs) + '</p>'

    annotations_html, or_target_motifs_html, cand_target_motifs_html, pl_alignment = report(plasmid_sequence, out_def, plasmid_annotations, input_patterns, f_patterns, type4RM = TYPE4, other_text = other_text)

    report_html = report_html_.replace("__ADD_ANNOTATIONS__", str(annotations_html))
    report_html = report_html.replace("__ADD_OR_TARGET_MOTIFS", str(or_target_motifs_html))
    report_html = report_html.replace(
        "__ADD_CAND_TARGET_MOTIFS", str(cand_target_motifs_html))
    report_html = report_html.replace(
        "__ADD_MODIFIED_POSITIONS", str(modified_positions))
    report_html = report_html.replace("__ADD_PLASMID_ALIGNMENT__", str(pl_alignment))

    with open(os.path.join(args.output_folder, f'report__{tmp_auto_values[args.auto]}__.html'), 'w') as handle:
        handle.write(report_html)

    return

# Functionality

def codon_usage(genome, codonTable):

    '''
    This function computes codon frequencies and normalization over synonimous codons --> percentages --> rank --> rank/highest rank

    The function takes as input a concatenation of all the CDS sequences in 5'-3' orientation

    The output is a table that contains for each codon the following scores:

    Count: is the total number of occurrences of the corresponding codon
    Proportion: is the normalized number of occurrences of each codon 
                by the total number of counts of the synonymous codons
                (those that encode for the same amino acid - included
                the one considered)
    Ranking: is the corresponding rank for each codon ( the lowest is assigned to the codon with the highest proportion value)
    Ranking_ratio: is the ratio between the ranking value for the corresponding codon and the rank of the most used 
                    codon between the synonymous codons (considering also the same codon). 
                    This value is relevant because it permits to consider also the number of synonymous codons for each amino
                    acid when considering the codon bias as an aggregated value.

    AA	Codon	Count	Proportion	Ranking	Ranking_ratio
    T	 ACT	16266	   0.46        2        0.5
    T	 ACC	1085	   0.03  	   3        0.75
    T	 ACA	17386	   0.5	       1        0.25
    T	 ACG	344	       0.01	       4        1
    '''

    codon_usage = {}

    tmp = [x for x in re.split(r'(\w{3})', str(genome)) if x != '']

    b_cod_table = CodonTable.unambiguous_dna_by_name['Bacterial'].forward_table

    for cod in CodonTable.unambiguous_dna_by_name['Bacterial'].stop_codons:
        b_cod_table[cod] = '_Stop'

    aas = set(b_cod_table.values())

    for aa in aas:
        codon_usage[aa] = {}
        for codon in b_cod_table.keys():
            if b_cod_table[codon] == aa:
                codon_usage[aa][codon] = tmp.count(codon.split(' ')[0])

    tups = {(outerKey, innerKey): values for outerKey, innerDict in codon_usage.items()
            for innerKey, values in innerDict.items()}

    codon_usage_ = pd.DataFrame(pd.Series(tups), columns=['Count'])
    codon_usage_.index = codon_usage_.index.set_names(['AA', 'Codon'])
    codon_usage_['Proportion'] = codon_usage_.groupby(level=0).transform(lambda x: (x / x.sum()).round(2))
    codon_usage_.reset_index(inplace=True)
    codon_usage_['Ranking'] = codon_usage_['Proportion']
    codon_usage_['Ranking'] = codon_usage_.groupby("AA")["Count"].rank("dense", ascending=False)
    codon_usage_['Ranking_ratio'] = codon_usage_.groupby("AA")['Ranking'].transform(lambda x: (x / len(x)).round(2))

    return {'Dictionary': codon_usage, 'Tuples': tups, 'Table': codon_usage_}

def read_rebase_output(input_file):
    '''
    This function reads the rebase output.
    The output is a dictionary in which:

    -meth_base: position (numbering starts from 0) of the methylated base in 5'-3' orientation
    -comp_meth_base: position (numbering starts from 0) of the methylated base in 3'-5' orientation
    -meth_type: type of methylation on the base indicated as meth_base in the 5'-3' orientation
    -comp_meth_type: type of methylation on the base indicated as meth_base in the 3'-5' orientation

    **Methylation type does not correspond to the type of the RM system**
    **Methylation type indicates the position on the G/T molecule on which the modification is performed**
    '''

    lines = {0: {}}
    enz_types = []
    counter = 0
    with open(input_file, 'r') as handle:
        for line in handle.readlines():
            if '<>' in line:
                counter += 1
                lines[counter] = {}
            if 'enz_type' in line:
                enz_types.append(int(line.strip().split('>')[1]))
            if ((len(line) > 1) & ("*" not in line) & ("<>" not in line) & ("org" not in line) & ("genome" not in line) & (("rec_seq" in line) | ("meth" in line))):
                lines[counter][line.strip().split('>')[0][1:]] = line.strip().split('>')[1]

    n_dict = {}

    for key in lines.keys():
        if ((lines[key] != {}) & ('rec_seq' in lines[key].keys())):
            if sum([1 for el in lines[key]['rec_seq'] if el in ['A', 'C', 'G', 'T']]) >= 2:
                n_dict[key] = lines[key]

    for key1 in n_dict.keys():
        for key2 in {'meth_base', 'comp_meth_base'} & set(n_dict[key1].keys()):
            n_dict[key1][key2] = int(n_dict[key1][key2])-1

    for key1 in n_dict.keys():
        if (not 'meth_base' in n_dict[key1].keys()):
            n_dict[key1]['meth_base'] = -math.inf
        if (not 'meth_type' in n_dict[key1].keys()):
            n_dict[key1]['meth_type'] = -math.inf
        if (not 'comp_meth_base' in n_dict[key1].keys()):
            n_dict[key1]['comp_meth_base'] = -math.inf
        if (not 'comp_meth_type' in n_dict[key1].keys()):
            n_dict[key1]['comp_meth_type'] = -math.inf
            
    return n_dict, set(enz_types)

def basic_rm_check(file_, from_str=False):

        if not from_str:
            with open(file_, 'r') as handle:
                lines = ''.join(handle.readlines())
        else:
            lines = file_

        out_lines = lines.split('<>')

        intereseting_lines = ['enz_type', 'sub_type', 'enz_name', 'rec_seq', 'meth_base', 'meth_type', 'comp_meth_base', 'comp_meth_type']

        out_dict = {}
        i = 0

        for enz in out_lines:

            if "PacBio records assigned to enzymes" in enz:
                assigned = "PacBio records assigned to enzymes"
            elif "PacBio records not yet assigned to enzymes" in enz:
                assigned = 'PacBio records not yet assigned to enzymes'
            elif "Enzymes not assigned to PacBio records" in enz:
                assigned = 'Enzymes not assigned to PacBio records'

            tmp_lines = enz.split('\n')
            
            out_dict[i] = {}

            if "assigned" in locals():
                out_dict[i]["assigned"] = assigned

            for line in tmp_lines:
                for inter in intereseting_lines:
                    if inter in line:
                        out_repr = [el for el in '<'.join(line.split('>')).split('<') if el != '']
                        out_dict[i][out_repr[0]] = out_repr[1]

            curr_repr = {}

            if 'rec_seq' in out_dict[i].keys():

                for j in range(len(out_dict[i]['rec_seq'].strip())):
                    curr_repr[j] = [out_dict[i]['rec_seq'].strip()[j]]

                if 'meth_base' in out_dict[i].keys():
                    if out_dict[i]['meth_base'] != '':
                        curr_repr[int(out_dict[i]['meth_base']) - 1].append('+')

                    if 'meth_type' in out_dict[i].keys():

                        if out_dict[i]['meth_type'] != '':
                            curr_repr[int(out_dict[i]['meth_base']) - 1].append('m' + out_dict[i]['meth_type'])

                if 'comp_meth_base' in out_dict[i].keys():
                    if out_dict[i]['comp_meth_base'] != '':
                        curr_repr[int(out_dict[i]['comp_meth_base']) - 1].append('-')
                
                    if 'comp_meth_type' in out_dict[i].keys():

                        if out_dict[i]['comp_meth_type'] != '':
                            curr_repr[int(out_dict[i]['comp_meth_base']) - 1].append('m' + out_dict[i]['comp_meth_type'])

                out_dict[i]['curr_repr'] = ' '.join([''.join(curr_repr[key][::-1]) for key in curr_repr.keys()])

            i += 1

        out_file = {}

        for key1 in out_dict.keys():
            for key2 in out_dict[key1].keys():
                if not key2 in out_file.keys():
                    out_file[key2] = []

        for key1 in out_dict.keys():

            for key2 in out_dict[key1].keys():
                out_file[key2].append(out_dict[key1][key2])

            not_considered = set(out_file.keys()) - set(out_dict[key1].keys())

            for key2 in not_considered:
                out_file[key2].append('')

        out_file = pd.DataFrame.from_dict(out_file)

        return out_file

def check_rm_information(input_, type="output_rebase"):

    # Read the target motifs

    out_rm_information = {
        'correct': True,
        'messages': [],
        'out_rm_information': ''
    }

    if not os.path.exists(input_):
        out_rm_information['correct'] = out_rm_information['correct'] & False
        out_rm_information['messages'].append('Indicated input genome {input_} does not exist!')
        return out_rm_information

    if os.path.isfile(input_):
        out_rebase = input_
    else:
        out_rebase = glob.glob(input_ + '/*')

        if not len(out_rebase) == 1:
            out_rm_information['correct'] = out_rm_information['correct'] & False
            out_rm_information['messages'].append('Number of REBASE files different from 1!')
            return out_rm_information
        else:
            out_rebase = out_rebase[0]

    out_rm_information['out_rm_information'] = out_rebase

    tmp_check = basic_rm_check(out_rebase)

    if tmp_check.shape[0] == 0:
        out_rm_information['correct'] = out_rm_information['correct'] & False
        out_rm_information['messages'].append('Incorrect input RM file!')
        return out_rm_information
    
    input_rebase_output = read_rebase_output(out_rebase)

    if not len(input_rebase_output[0]) > 0:
        out_rm_information['correct'] = out_rm_information['correct'] & False
        out_rm_information['messages'].append(f'Error in {out_rebase}. Incorrect input: no target motif has been found.')
        return out_rm_information

    if not len(input_rebase_output[1]) > 0:
        out_rm_information['correct'] = out_rm_information['correct'] & False
        out_rm_information['messages'].append(f'Error in {out_rebase}. Incorrect input: no RM system type detected.')
        return out_rm_information

    return out_rm_information

def report(sequence, new_sequence, annotations, target_motifs, f_patterns, type4RM = False, other_text = ''):

    index = list(range(len(sequence)))

    representation = {}

    for i in range(len(sequence)):
        representation[i] = sequence[i]

    representation_ = {}

    for i in range(len(new_sequence)):
        if new_sequence[i] != sequence[i]:
            representation_[i] = "<b><orange>{}</orange></b>".format(new_sequence[i])
        else:
            representation_[i] = new_sequence[i]

    annotations_html = ""
    for feat in annotations:
        annotations_html += '<p>&nbsp;&nbsp;&nbsp;&nbsp;' + repr(feat) + '</p>\n'
    
    or_target_motifs_html = ""

    for f_pattern in f_patterns.keys():
        or_target_motifs_html += '<p>&nbsp;&nbsp;&nbsp;&nbsp;' + f_pattern + ':</p>\n'
        for val in f_patterns[f_pattern]:
            or_target_motifs_html += '<p>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;' + \
                str(val) + ':</p>\n'


    if type4RM:

        or_target_motifs_html += '<p>DETECTED Type IV RM System. Do not use MTases and pass the plasmid/minicircle through demethylating bacterial strain!</p>\n'
        or_target_motifs_html += "</span>\n<br>\n<span class='comparison'>\n"

    cand_target_motifs_html = ""

    cand_target_motifs_html += other_text
    cand_target_motifs_html += "</span>\n<br>\n<span class='comparison'>\n"    


    pl_alignment = ""

    cds = {}
    for ann in annotations:
        cds[repr(ann)] = {}
        direction = ann[-2]
        range_ = [int(ann[0]), int(ann[1])]
        for i in range(range_[0]):
            cds[repr(ann)][i] = '&nbsp;'
        if direction == "+":
            j = 0
            for i in range(range_[0], range_[1], 3):
                if j == 0:
                    cds[repr(ann)][i] = "<b><orange>></orange></b>"
                    cds[repr(ann)][i+1] = "<b><orange>></orange></b>"
                    cds[repr(ann)][i+2] = "<b><orange>></orange></b>"
                    j = 1
                else:
                    cds[repr(ann)][i] = "<b><green>></green></b>"
                    cds[repr(ann)][i+1] = "<b><green>></green></b>"
                    cds[repr(ann)][i+2] = "<b><green>></green></b>"
                    j = 0
        else:
            j = 0
            for i in range(range_[0], range_[1], 3):
                if j == 0:
                    cds[repr(ann)][i] = "<b><orange><</orange></b>"
                    cds[repr(ann)][i+1] = "<b><orange><</orange></b>"
                    cds[repr(ann)][i+2] = "<b><orange><</orange></b>"
                    j = 1
                else:
                    cds[repr(ann)][i] = "<b><green><</green></b>"
                    cds[repr(ann)][i+1] = "<b><green><</green></b>"
                    cds[repr(ann)][i+2] = "<b><green><</green></b>"
                    j = 0
        for i in range(range_[1], len(sequence)):
            cds[repr(ann)][i] = '&nbsp;'

    tmp1 = 'Original_seq:&nbsp;&nbsp;&nbsp;'
    tmp1_ ='CandSyng_seq:&nbsp;&nbsp;&nbsp;'
    tmp2 = 'PlasmidIndex:&nbsp;&nbsp;&nbsp;'
    tmp3 = '_'*len(tmp2)
    repr_cds = {}
    for a in cds.keys():
        repr_cds[a] = 'CodingSequen:&nbsp;&nbsp;&nbsp;'
    for pos in sorted(representation.keys()):
        tmp1 += '&nbsp;'*(len(str(pos))-1) + representation[pos] + '&nbsp;'
        tmp1_ += '&nbsp;'*(len(str(pos))-1) + representation_[pos] + '&nbsp;'
        tmp3 += '_'*((len(str(pos))-1)+len(representation[pos])+1)
        for a in repr_cds.keys():
            repr_cds[a] += '&nbsp;'*(len(str(pos))-1) + \
                cds[a][pos] + '&nbsp;'
        tmp2 += str(pos) + '&nbsp;'
        if (pos > 0) & (pos % 18 == 0):
            pl_alignment += '<br><br>'
            tmp1 += '<p>'
            tmp1_ += '<p>'
            tmp2 += '<p>'
            tmp3 += '<p>'
            for a in repr_cds.keys():
                repr_cds[a] += '<p>'
            for a in repr_cds.keys():
                pl_alignment += '<p>' + repr_cds[a]
            pl_alignment += '<p>' + tmp1 + '<p>' + tmp1_ + '<p>' + tmp2 + '<p>' + tmp3
            tmp1 = 'Original_seq:&nbsp;&nbsp;&nbsp;'
            tmp1_ ='CandSyng_seq:&nbsp;&nbsp;&nbsp;'
            tmp2 = 'PlasmidIndex:&nbsp;&nbsp;&nbsp;'
            tmp3 = '_'*len(tmp2)
            for a in repr_cds.keys():
                repr_cds[a] = 'CodingSequen:&nbsp;&nbsp;&nbsp;'
    tmp1 += '<p>'
    tmp1_ += '<p>'
    tmp2 += '<p>'
    for a in repr_cds.keys():
        repr_cds[a] += '<p>'
    for a in repr_cds.keys():
        pl_alignment += '<p>' + repr_cds[a]
    pl_alignment += tmp1 + '<p>' + tmp2

    pl_alignment += '</span>\n<br>\n</body>\n</html>'

    return annotations_html, or_target_motifs_html, cand_target_motifs_html, pl_alignment

def run_estimator(args, type_action='syngeneic_estimator'):

    type_actions = ['mapping_targets', 'partial_generalization', 'syngeneic_estimator']

    if not type_action in type_actions:
        info('Incorrect type of action!')
        return 1

    if type_action == 'partial_generalization':
        if not os.path.exists(args.synpl_folder):
            error(f'Input folder for temporary mapping {args.synpl_folder} does not exists!')

        existing_temporary_mapping_files = glob.glob(os.path.join(args.synpl_folder, "*_temporary_mapping.synpl"))

        if not len(existing_temporary_mapping_files) > 0:
            error('No temporary mapping files has been found!')

    start_time = time.time()

    log_file = os.path.join(args.output_folder, f'sytogen_{type_action}.log')

    if os.path.exists(log_file):
        os.remove(log_file)

    logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', filename=os.path.join(args.output_folder, f'sytogen_{type_action}.log'), level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')

    logging.info("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    logging.info("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    logging.info(f"Starting {type_action}")
    logging.info("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    logging.info("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
    logging.info("--- Start time: %s seconds ---" % (time.ctime()))
    if args.verbose:
        info("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        info("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        info(f"Starting {type_action}")
        info("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        info("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        info("--- Start time: %s seconds ---" % (time.ctime()))
       
    codonTable = args.codon_table

    # Read the input plasmid

    input_plasmid = read_input_seq(args.input_sequence)

    if not input_plasmid['correct']:
        info('The input plasmid is not correct! Please, check it.')
        error('\n'.join(input_plasmid['messages']))
    else:
        if len(input_plasmid['annotations'].keys()) == 1:
            if len(input_plasmid['annotations'][0]['sequence']) > 0:
                plasmid_sequence = input_plasmid['annotations'][0]['sequence']
                plasmid_annotations = input_plasmid['annotations'][0]['cds']
                l = len(plasmid_sequence)
            else:
                error('More sequences in the input plasmid file OR the sequence of the input file is empty!')
        else:
            error('Incorrect number of inputs for the plasmid sequence!')

    # Read the input genome

    input_genome = read_input_seq(args.input_strain_genome)

    config = configparser.ConfigParser()

    if not input_genome['correct']:
        info('The input genome is not correct! Please, find more info in the log file.')
        logging.info("")
        logging.info(
            "+------------------------------------------------------------------------------------------------------------")
        logging.info('\n'.join(input_genome['messages']))
        logging.info(
            "+------------------------------------------------------------------------------------------------------------")
        logging.info("")

        if not len(input_genome['total_cds_concats']) > 1000:
            error('\n'.join(input_genome['messages']))
        
    
    genome_codon_usage = codon_usage(input_genome['total_cds_concats'], codonTable)
        
    with open(os.path.join(args.output_folder, "codon_usage_table.csv"), 'w') as outf:
        genome_codon_usage['Table'].to_csv(outf, sep=',', index=False)

    genome_codon_usage_internal_repr = {
        el: {'Percentage': str(round(float(genome_codon_usage['Table'][genome_codon_usage['Table']['Codon'] == el]['Proportion'])*100, 1))+"%",
                'Inverse_ranking': float(genome_codon_usage['Table'][genome_codon_usage['Table']['Codon'] == el]['Ranking_ratio']),
                'Translation': [genome_codon_usage['Table'][genome_codon_usage['Table']['Codon'] == el]['AA'].values[0] if genome_codon_usage['Table'][genome_codon_usage['Table']['Codon'] == el]['AA'].values[0] != '_Stop' else '*'][0]} for el in genome_codon_usage['Table']['Codon']}
    logging.info(genome_codon_usage_internal_repr)
    
    input_rebase_output = read_rebase_output(args.input_rm_systems)
    input_patterns = [input_rebase_output[0][el]['rec_seq']
                      for el in input_rebase_output[0].keys()]
    max_len = max([len(inp_pat) for inp_pat in input_patterns])

    if 4 in input_rebase_output[1]:

        logging.info("")
        logging.info(
            "+------------------------------------------------------------------------------------------------------------")
        logging.info(
            "**TYPE IV RM Systems are found in the target motifs of the species!**")
        logging.info(
            "Demethylate the plasmid before use and do not use further methylation!")
        logging.info(
            "Produce the syngenic plasmid only through syngenic modifications!")
        logging.info(
            "+------------------------------------------------------------------------------------------------------------")
        logging.info("")

        if args.verbose:

            info("")
            info("+------------------------------------------------------------------------------------------------------------")
            info("**TYPE IV RM Systems are found in the target motifs of the species!**")
            info(
                "Demethylate the plasmid before use and do not use further methylation!")
            info("Produce the syngenic plasmid only through syngenic modifications!")
            info("+------------------------------------------------------------------------------------------------------------")
            info("")
    else:
        
        logging.info("")
        logging.info("+------------------------------------------------------------------------------------------------------------")
        logging.info("**TYPE IV RM Systems not found in the target motifs of the species!**")
        logging.info("+------------------------------------------------------------------------------------------------------------")
        logging.info("")

        if args.verbose:

            info("")
            info("+------------------------------------------------------------------------------------------------------------")
            info("**TYPE IV RM Systems not found in the target motifs of the species!**")
            info("+------------------------------------------------------------------------------------------------------------")
            info("")

    # Read commercial MTases

    commercial_mtases = read_commercial(camth)

    logging.info("")
    logging.info("The found RM systems are: " + str(input_patterns))

    logging.info("")
    logging.info("The found input patterns are: " + str(input_rebase_output))

    logging.info("")
    logging.info("The plasmid length is: " + str(l))

    logging.info("")
    logging.info("The longest target motif is: " + str(max_len))

    if args.verbose:
        info("The input plasmid is at location: " + args.input_sequence)
        info("The output from REBASE is at location: " +
              args.input_rm_systems)
        info("The species genome is at location: " + args.input_strain_genome)
        info("The indicated codon table is: " + str(args.codon_table))
        info("The indicated output folder is: " + args.output_folder)

        info("")
        info("The found RM systems are: " + str(input_patterns))

        info("")
        info("The found input patterns are: " + str(input_rebase_output))

        info("")
        info("The plasmid length is: " + str(l))

        info("")
        info("The longest target motif is: " + str(max_len))

    patterns_found_plasmid = findpatterns(plasmid_sequence, input_patterns, circular=True)

    logging.info("")
    logging.info("The number of target motifs found on the plasmid is: ")
    logging.info(len(sum(list(patterns_found_plasmid.values()), [])))
    logging.info("The target motifs found on the plasmid are: ")
    logging.info(patterns_found_plasmid)

    if args.verbose:
            info("")
            info("The target motifs found on the plasmid are: ")
            info(patterns_found_plasmid)
    
    if not len(sum(list(patterns_found_plasmid.values()), [])) > 0:
        info("")
        info(f"The number of target motifs found on the plasmid is: {len(sum(list(patterns_found_plasmid.values()), []))}")
        return

    out_target_slices = punctuate_targets(patterns_found_plasmid, l, max_len)
#     print("OOOUT:", out_target_slices)
#     raise Exception("Here")
#     logging.info("LOOK HERE: ", str(out_target_slices))
#     raise Exception("Here")

    if type_action == 'partial_generalization':

        temporary_file_target_positions_to_consider_absolute = []

        candidate_alternatives_to_consider = {}
        existing_temporary_mapping_ranges = sorted([(int(re.sub(r'.*[/\\]slice|_temporary_mapping.synpl', '', el).split('_')[0]), int(re.sub(r'.*[/\\]slice|_temporary_mapping.synpl', '', el).split('_')[1])) for el in existing_temporary_mapping_files])
        newly_generated_mapping_ranges = sorted(list(set([(el[0], el[-1]+1) for el in out_target_slices[0]]))); print(sorted(list(newly_generated_mapping_ranges)))

        tmp_existing_intersection = set.intersection(set(existing_temporary_mapping_ranges), set(newly_generated_mapping_ranges)); print(sorted(list(tmp_existing_intersection)))
        missing_temporary_mapping_files = sorted(list(set.difference(
            set(newly_generated_mapping_ranges), set(existing_temporary_mapping_ranges))))
        
        print(" ------------")
        print(existing_temporary_mapping_ranges)
        print(newly_generated_mapping_ranges)
        print(tmp_existing_intersection)
        print(missing_temporary_mapping_files)
        
        # raise Exception("Here")
        
        if (
            (len(tmp_existing_intersection) != len(existing_temporary_mapping_ranges)) | \
            (len(tmp_existing_intersection) != len(newly_generated_mapping_ranges)) ):
            error(f'Missing temporary files {missing_temporary_mapping_files}. The temporary files must correspond to the inner generated slices ranges.')

        number_of_user_selection = 0

        for i in range(len(out_target_slices[0])):

            tmp_range = out_target_slices[0][i]
            tmp_slice_sequence = extract_subsequence(plasmid_sequence, tmp_range)
            tmp_slice_targets = findpatterns(tmp_slice_sequence, input_patterns, circular=False)
        
            slice_non_ambiguous_positions_plasmid = out_target_slices[1][i]
            print(tmp_range, slice_non_ambiguous_positions_plasmid)
            slice_non_ambiguous_positions_slice = [[tmp_range.index(el1) for el1 in el0] for el0 in slice_non_ambiguous_positions_plasmid]
            out_temporary_mapping_file = read_synpl(os.path.join(args.synpl_folder, 'slice' + str(tmp_range[0]) + '_' + str(tmp_range[-1]+1) + '_temporary_mapping.synpl'), len(plasmid_sequence), mapping=True)

            # Perform check of the temporary mapping

            tmp_tmp_file = os.path.join(args.synpl_folder, 'slice' + str(tmp_range[0]) + '_' + str(tmp_range[-1]+1) + '_temporary_mapping.synpl')

            if not out_temporary_mapping_file['correct']:
                info(out_temporary_mapping_file['messages'])
                error(f'Incorrect temporary mapping (.synpl) file {tmp_tmp_file}.')

            if not 1 in out_temporary_mapping_file['output'].keys():
                info(out_temporary_mapping_file['messages'])
                error(f'Incorrect temporary mapping (.synpl) file {tmp_tmp_file}.')

            else:
                
                if not len(out_temporary_mapping_file['output'][1]):
                    info(out_temporary_mapping_file['messages'])
                    error(f'Incorrect temporary mapping (.synpl) file {tmp_tmp_file}.')

                tmp_split_check = re.sub(r'\|-5\' - | -3\'- \(\+ strand\)', '', out_temporary_mapping_file['output'][1][-1]).strip()

                temporary_file_target_positions_to_consider_absolute.append([(tmp_range[0] + i)%l for i in range(len(tmp_split_check)) if tmp_split_check[i].isupper()])

        temporary_file_target_positions_to_consider_absolute = sum(temporary_file_target_positions_to_consider_absolute, [])
        info(temporary_file_target_positions_to_consider_absolute)
        # out_target_slices = punctuate_targets(temporary_file_target_positions_to_consider_absolute, l, max_len, is_pos_list = True)
        print(' ================================ ')
        print(tmp_split_check)
        
        # raise Exception("Here")

    logging.info('')
    logging.info('The number of generated slices is: ' + str(len(out_target_slices[0])))
    logging.info('The ranges of generated slices are: ' + str([(el[0], el[-1] + 1) for el in out_target_slices[0]]))

    if args.verbose:
        info('')
        info('The number of generated slices is: ' +
              str(len(out_target_slices[0])))
        info('The ranges of generated slices are: ' + str([(el[0], el[-1] + 1) for el in out_target_slices[0]]))

    minima = []
    maxima = []

    slices_evaluation = {}

    num_slices = 0
    
    for i in range(len(out_target_slices[0])):

        if num_slices == 0:
                logging.info("")
                logging.info(
                    "+------------------------------------------------------------------------------------------------------------")
                logging.info(
                    "Starting slice generalization, alternative generation and alternative selection")
                logging.info(
                    "+------------------------------------------------------------------------------------------------------------")
                logging.info("")

                if args.verbose:
                    info("")
                    info("+------------------------------------------------------------------------------------------------------------")
                    info("Starting slice generalization, alternative generation and alternative selection")
                    info("+------------------------------------------------------------------------------------------------------------")
                    info("")
        
        num_slices += 1

        tmp_range = out_target_slices[0][i]
        tmp_slice_sequence = extract_subsequence(plasmid_sequence, tmp_range)
        tmp_slice_targets = findpatterns(tmp_slice_sequence, input_patterns, circular=False)
    
        slice_non_ambiguous_positions_plasmid = out_target_slices[1][i]
        print(out_target_slices[1])

        if type_action == 'mapping_targets':
            mapping_targets_or = single_slice_mapping(tmp_slice_sequence, tmp_range, plasmid_sequence, plasmid_annotations,
                            genome_codon_usage_internal_repr, commercial_mtases, input_patterns, input_rebase_output, candidate_alternative=True)
            mapping_targets_mod = single_slice_mapping(tmp_slice_sequence, tmp_range, plasmid_sequence, plasmid_annotations,
                                                    genome_codon_usage_internal_repr, commercial_mtases, input_patterns, input_rebase_output, candidate_alternative=False)

            tmp_str = '_____________________________________________________________________\n_____________________________________________________________________\n'
            with open(os.path.join(args.output_folder, 'slice' + str(tmp_range[0]) + '_' + str(tmp_range[-1]+1) + '_temporary_mapping.synpl'), 'w') as handle:
                handle.write(synpl_map_header)
                handle.write(mapping_targets_or + tmp_str + mapping_targets_mod)

            continue

        if type_action == 'partial_generalization':
            print(tmp_range)
            print(slice_non_ambiguous_positions_plasmid)
            slice_non_ambiguous_positions_slice = [[tmp_range.index(el) for el in sum(slice_non_ambiguous_positions_plasmid, [])]]
        
        else:
            slice_non_ambiguous_positions_slice = [[tmp_range.index(el1) for el1 in el0] for el0 in slice_non_ambiguous_positions_plasmid]


        logging.info('')
        logging.info("+------------------------------------------------------------------------------------------------------------")
        
        logging.info('Considering slice: [' + str(tmp_range[0]) + ', ' + str(
            tmp_range[-1]+1) + ']; ' + str(i+1) + '/' + str(len(out_target_slices[0])) +'; len: ' + str(len(tmp_range)))
        logging.info('Non-ambiguous target positions on the plasmid: ' +
                     repr(slice_non_ambiguous_positions_plasmid))
                     
        logging.info('Non-ambiguous target positions on the slice: ' + repr(
            slice_non_ambiguous_positions_plasmid))
        logging.info('Target plasmid sequence: ' + repr(tmp_slice_sequence))
        logging.info('Patterns found on the considered slice: ' +
                     repr({el:tmp_slice_targets[el] for el in tmp_slice_targets.keys() if tmp_slice_targets[el] != []}))

        logging.info('Cmin = ' + str((4**len(slice_non_ambiguous_positions_slice) *
                                      reduce(mul, [len(el) for el in slice_non_ambiguous_positions_slice], 1))))
        logging.info(
            'Cmax = ' + str((4**sum([len(el) for el in slice_non_ambiguous_positions_slice]))))

        if args.verbose:
            info('')
            info("+------------------------------------------------------------------------------------------------------------")
            info('Considering slice: [' + str(tmp_range[0]) + ', ' + str(tmp_range[-1]+1) + ']; ' + str(
                i+1) + '/' + str(len(out_target_slices[0])) + '; len: ' + str(len(tmp_range)))
            info('Non-ambiguous target positions on the plasmid: ' + repr(slice_non_ambiguous_positions_plasmid))
            info('Non-ambiguous target positions on the slice: ' + repr(slice_non_ambiguous_positions_plasmid))
            info('Target plasmid sequence: ' + repr(tmp_slice_sequence))
            info('Patterns found on the considered slice: ' + repr({el: tmp_slice_targets[el] for el in tmp_slice_targets.keys() if tmp_slice_targets[el] != []}))
            info('Cmin = ' + str((4**len(slice_non_ambiguous_positions_slice) * reduce(mul, [len(el) for el in slice_non_ambiguous_positions_slice], 1))))
            info('Cmax = ' + str((4**sum([len(el) for el in slice_non_ambiguous_positions_slice]))))

        minima.append((4**len(slice_non_ambiguous_positions_slice) * reduce(mul, [len(el) for el in slice_non_ambiguous_positions_slice], 1)))
        maxima.append(
            (4**sum([len(el) for el in slice_non_ambiguous_positions_slice])))

        minimum_c = ((4**len(slice_non_ambiguous_positions_slice)) * reduce(mul, [len(el) for el in slice_non_ambiguous_positions_slice], 1))
        maximum_c = (4**sum([len(el) for el in slice_non_ambiguous_positions_slice]))

        generalized_slices = []

        if maximum_c <= 4**6:
            
            combs_positions_to_generalize = sum(slice_non_ambiguous_positions_slice, [])
            n_slice_ = {i: tmp_slice_sequence[i]
                        for i in range(len(tmp_slice_sequence))}
            
            for position in combs_positions_to_generalize:
                n_slice_[position] = '(A|T|C|G)'
            nn_slice_ = ''
            for index in sorted(n_slice_.keys()):
                nn_slice_ += n_slice_[index]
            generalized_slices.append(nn_slice_)

            logging.info("")
            logging.info("Cmax <= 4**6 --> Do complete generalization.")
            logging.info("The produced generalized slice is: " + nn_slice_)

            if args.verbose:
                info("")
                info("Cmax <= 4**6 --> Do complete generalization.")
                info("The produced generalized slice is: " + nn_slice_)

        elif (maximum_c > 4**6) & (minimum_c <= 4**6):

            combs_positions_to_generalize = combn(slice_non_ambiguous_positions_slice, len(slice_non_ambiguous_positions_slice))
            
            for comb in combs_positions_to_generalize:
                    n_slice_ = {i: tmp_slice_sequence[i]
                                for i in range(len(tmp_slice_sequence))}
                    for position in comb:
                        n_slice_[position] = '(A|T|C|G)'
                    nn_slice_ = ''
                    for index in sorted(n_slice_.keys()):
                        nn_slice_ += n_slice_[index]
                    generalized_slices.append(nn_slice_)

            logging.info("")
            logging.info("Cmax > 4**6 & Cmin <= 4**6 --> Do partial generalization (use of combn function).")
            logging.info("Produced positions to generalize are: " + str(combs_positions_to_generalize))
            logging.info("The produced generalized slice is: " + str(generalized_slices))

            if args.verbose:
                info("")
                info("Cmax > 4**6 & Cmin <= 4**6 --> Do partial generalization (use of combn function).")

            logging.info(combs_positions_to_generalize)

        elif (minimum_c > 4**6) & (len(slice_non_ambiguous_positions_slice) <= 7):
            combs_positions_to_generalize = [
                [random.choice(el) for el in slice_non_ambiguous_positions_slice]]
            for comb in combs_positions_to_generalize:
                n_slice_ = {i: tmp_slice_sequence[i]
                            for i in range(len(tmp_slice_sequence))}
                for position in comb:
                    n_slice_[position] = '(A|T|C|G)'
                nn_slice_ = ''
                for index in sorted(n_slice_.keys()):
                    nn_slice_ += n_slice_[index]
                generalized_slices.append(nn_slice_)

            logging.info("")
            logging.info(
                "Cmin > 4**6 & number of non-ambiguous positions <= 7 --> Do random generalization")
            logging.info("Produced positions to generalized are: " + \
                            str(combs_positions_to_generalize))
            logging.info("The produced generalized slice is: " + \
                            str(generalized_slices))

            if args.verbose:
                info("")
                info(
                    "Cmin > 4**6 & number of non-ambiguous positions <= 7 --> Do random generalization")

        else:
            logging.info("")
            logging.info("Cmin > 4**6 & number of non-ambiguous positions > 7 --> The slice is too complex to be considered.")

            if args.verbose:
                info("")
                info("Cmin > 4**6 & number of non-ambiguous positions > 7 --> The slice is too complex to be considered.")

        nn_slice_ = []
        if len(generalized_slices) > 0:
            for generalized_slice in generalized_slices:
                nn_slice_ += list(sre_yield.AllStrings(generalized_slice))

        logging.info("")
        logging.info("\+------------------------------------------------------------------------------------------------------------")
        logging.info("\Starting alternative generation")
        logging.info("\+------------------------------------------------------------------------------------------------------------")
        logging.info("")

        if args.verbose:
            info("")
            info("\+------------------------------------------------------------------------------------------------------------")
            info("\Starting alternative generation")
            info("\+------------------------------------------------------------------------------------------------------------")
            info("")

        logging.info("")
        logging.info(
            "The cardinality of the expanded set of slices is: " + str(len(nn_slice_)))

        if args.verbose:
            info("")
            info("The cardinality of the expanded set of slices is: " + str(len(nn_slice_)))

        logging.info("--- Partial time: %s seconds ---" % (time.ctime()))
        logging.info("--- %s hh:mm:ss employed ---" %
              (datetime.timedelta(seconds=round(time.time() - start_time, 0))))

        if args.verbose:
            info("--- Partial time: %s seconds ---" % (time.ctime()))
            info("--- %s hh:mm:ss employed ---" %
              (datetime.timedelta(seconds=round(time.time() - start_time, 0))))
        out_evaluation = evaluate_slices(
            or_slice=tmp_slice_sequence, 
            candidate_alternatives=nn_slice_, 
            slice_range=tmp_range, 
            current_input_genetic_tool_sequence=plasmid_sequence, 
            current_input_genetic_tool_annotations=plasmid_annotations, 
            codon_usage_table=genome_codon_usage_internal_repr, 
            commercial_mtases=commercial_mtases, 
            target_motifs=input_patterns,  
            input_rebase_output=input_rebase_output, ncores=args.NCORES)
        tmp_str = '_____________________________________________________________________\n_____________________________________________________________________\n'.join(
            [el[-1] for el in out_evaluation])
        with open(os.path.join(args.output_folder, 'slice' + str(tmp_range[0]) + '_' + str(tmp_range[-1]+1) + '.synpl'), 'w') as handle:
            handle.write(synpl_header)
            handle.write(tmp_str)
            
    logging.info("--- End time: %s seconds ---" % (time.ctime()))
    logging.info("--- %s hh:mm:ss employed ---" %
          (datetime.timedelta(seconds=round(time.time() - start_time, 0))))
    if args.verbose:
        info("--- End time: %s seconds ---" % (time.ctime()))
        info("--- %s hh:mm:ss employed ---" %
                 (datetime.timedelta(seconds=round(time.time() - start_time, 0))))
        info('\n')

    return

# Functionality

def run_codon_bias(args):

    input_genome = args.input_strain_genome
    codonTable = args.codon_table
    input_genome = read_input_seq(input_genome)

    if not input_genome['correct']:
        info('The input genome is not correct! Please, find more info in the log file.')
        logging.info("")
        logging.info(
            "+------------------------------------------------------------------------------------------------------------")
        logging.info('\n'.join(input_genome['messages']))
        logging.info(
            "+------------------------------------------------------------------------------------------------------------")
        logging.info("")

        if not len(input_genome['total_cds_concats']) > 1000:
            error('\n'.join(input_genome['messages']))

    genome_codon_usage = codon_usage(input_genome['total_cds_concats'], codonTable)
        
    with open(os.path.join(args.output_folder, "codon_usage_table.csv"), 'w') as outf:
        genome_codon_usage['Table'].to_csv(outf, sep=',', index=False)

    return

# Functionality

def run_input_check(args):

    input_file = args.input_sequence
    output_folder = args.output_folder


    with open(input_file, 'r') as input_handle:
        tmp_lines = ''.join(input_handle.readlines())

    sequences = [el.strip().split('ORIGIN')
                for el in re.split(r'\s+//', tmp_lines) if el.strip() != '']

    k = 1
    i = 1

    sequence_n = []

    feature_n = []

    strand = []

    original_start = []
    original_end = []

    new_start = []
    new_end = []

    feat_length = []
    feat_type = []

    is_cds = []
    number_of_codons = []
    start_codon = []
    stop_codon = []

    on_join = []

    
    for seq in sequences:
        features = [el.replace('\n', '') for el in re.split('\n\s{5,5}', ''.join(seq[0].split(
            'Location/Qualifiers')[1].split('BASE COUNT')[0].strip().replace(' '*21, '__')))]
        current_sequence = re.sub('[\d\s]', '', seq[1]).upper()
        for el in features[1:]:

            sequence_n.append(k)
            feat = el.split('__')[0]

            feature_n.append(i)

            if not 'join' in feat:
                indeces = [int(el) for el in re.sub(
                    r'\D', ' ', feat).strip().split(' ') if el != '']
                on_join.append(False)

                or_start = indeces[0]
                or_end = indeces[-1]

                nw_start = indeces[0] - 1
                nw_end = indeces[-1]

            else:

                on_join.append(True)

                indeces = [int(el) for el in re.sub(
                    r'\D', ' ', feat).strip().split(' ') if el != '']
                indeces = [indeces[0], indeces[-1]]

                or_start = indeces[0]
                or_end = indeces[-1]

                nw_start = indeces[0] - 1
                nw_end = indeces[-1]

            irange_ = irange(i=nw_start, j=nw_end, l=len(current_sequence))

            tmp_sequence = ''.join([current_sequence[i] for i in irange_])

            original_start.append(or_start)
            original_end.append(or_end)

            new_start.append(nw_start)
            new_end.append(nw_end)

            feat_length.append(nw_end - nw_start)

            if 'complement' in feat:
                strand.append('-1')
            else:
                strand.append('+1')

            feat_type.append(feat.split(' ')[0])

            if feat.split(' ')[0] == 'CDS':
                is_cds.append(True)
                number_of_codons.append((nw_end - nw_start)/3)

                if 'complement' in feat:
                    start_codon.append(reverse_complement(tmp_sequence)[:3])
                    stop_codon.append(reverse_complement(tmp_sequence)[-3:])
                else:
                    start_codon.append(tmp_sequence[:3])
                    stop_codon.append(tmp_sequence[-3:])
            else:
                is_cds.append(False)
                number_of_codons.append('')
                start_codon.append('')
                stop_codon.append('')

            i += 1

        k += 1

    stats = pd.DataFrame(data={'Sequence_number': sequence_n,
                            'Number_of_features': feature_n,
                            'Strand': strand,
                            'On_join': on_join,
                            'Or_Feature_Start[INCLUDED]': original_start,
                            'Or_Feature_End[INCLUDED]': original_end,
                            'New_Feature_Start[INCLUDED]': new_start,
                            'New_Feature_End[EXCLUDED]': new_end,
                            'Feature_length': feat_length,
                            'Feature_type': feat_type,
                            'Is_CDS': is_cds,
                            'Number_of_codons': number_of_codons,
                            'Start_codon': start_codon,
                            'Stop_codon': stop_codon})


    stats.to_excel(os.path.join(args.output_folder, 'input_gbk_check.xlsx'))

    return

# Functionality

def run_output_check(args):

    or_input_file = args.input_sequence
    nw_input_file = args.input_sequence_new
    
    input_rebase_output = read_rebase_output(args.input_rm_systems)
    input_patterns = [input_rebase_output[0][el]['rec_seq']
                      for el in input_rebase_output[0].keys()]

    codon_table = args.codon_table
    output_folder = args.output_folder

    with open(or_input_file, 'r') as input_handle:
        or_input_sequence = ''.join(input_handle.readlines())

    with open(nw_input_file, 'r') as input_handle:
        nw_input_sequence = ''.join(input_handle.readlines())

    or_header = or_input_sequence.split('source')[0]
    nw_header = nw_input_sequence.split('source')[0]

    or_annotations = [el.replace('\n', '') for el in re.split('\n\s{5,5}', ''.join(or_input_sequence.split(
        'Location/Qualifiers')[1].split('ORIGIN')[0].strip().replace(' '*21, '__')))]
    nw_annotations = [el.replace('\n', '') for el in re.split('\n\s{5,5}', ''.join(nw_input_sequence.split(
        'Location/Qualifiers')[1].split('ORIGIN')[0].strip().replace(' '*21, '__')))]

    or_sequence = re.sub(
        r'[\d\s//]', '', or_input_sequence.split('ORIGIN')[1]).upper()
    nw_sequence = re.sub(
        r'[\d\s//]', '', nw_input_sequence.split('ORIGIN')[1]).upper()
    
    patterns_found_plasmid = findpatterns(nw_sequence, input_patterns, circular=True)

    table = ''

    table_header = ['Original sequence', 'New_sequence']
    table_rows = ['Length of the sequence', 'Number']

    table += 'Original sequence lenght: ' + str(len(or_sequence)) + '\n'
    table += 'New sequence lenght: ' + str(len(nw_sequence)) + '\n'
    
    table += 'Number of left RM recognition motifs: ' + str(len(sum(list(patterns_found_plasmid.values()), []))) + '\n'

    differences_pos = [i+1 for i in range(len(or_sequence)) if or_sequence[i] != nw_sequence[i]]

    out_to_genbank = copy.copy(nw_input_sequence)
    for j in [i+1 for i in range(len(or_sequence)) if or_sequence[i] != nw_sequence[i]]:
        out_to_genbank = out_to_genbank.replace("ORIGIN", "     " + f"SyT             {j}\n                     /label=\"{or_sequence[j-1]} --> {nw_sequence[j-1]}\"\nORIGIN")

    table += 'Number of differing positions: ' + str(len(differences_pos)) + '\n'
    table += 'List of differing positions: ' + \
        str(differences_pos) + '\n'
    table += 'All annotations check: ' + str(or_annotations == nw_annotations) + '\n'
    table += 'All annotations [header] check: ' + str(or_header == nw_header) + '\n'

    for feat in or_annotations:
        if feat.startswith('CDS'):            

            indeces = [int(el) for el in re.sub(
                r'\D', ' ', feat.split('__')[0]).strip().split(' ') if el != '']
    
            tmp_start = indeces[0] - 1
            tmp_end = indeces[1]

            tmp_or_sequence = ''.join([or_sequence[i] for i in irange(i = tmp_start, j = tmp_end, l = len(or_sequence))])
            tmp_nw_sequence = ''.join([nw_sequence[i] for i in irange(i=tmp_start, j=tmp_end, l=len(nw_sequence))])

            if 'complement' in feat:
                tmp_or_sequence = reverse_complement(tmp_or_sequence)
                tmp_nw_sequence = reverse_complement(tmp_nw_sequence)

            or_tr_sequence = Seq(tmp_or_sequence).translate(table=codon_table)
            nw_tr_sequence = Seq(tmp_nw_sequence).translate(table=codon_table)

            table += 'Check same protein sequence feature \n\t`' + feat + '`: ' + \
                ' \n\tLength original: ' + str(len(tmp_or_sequence)) + \
                ' \n\tLength new: ' + str(len(tmp_nw_sequence)) + \
                ' \n\tSame proteic sequence: ' + \
                str(or_tr_sequence == nw_tr_sequence) + '\n'
    
    with open(os.path.join(args.output_folder, 'candidate_syngenic_check.txt'), 'w') as input_handle:
        input_handle.write(table)

    with open(os.path.join(args.output_folder, 'candidate_syngenic_check__sequence.gbk'), "w") as handle:
        handle.write(out_to_genbank)

    return

# Functionality

def run_tool_representation(args):

    input_file = args.input_sequence
    output_rebase = args.input_rm_systems
    output_folder = args.output_folder

    with open(input_file, 'r') as input_handle:
        tmp_lines = ''.join(input_handle.readlines())

    sequences = [el.strip().split('ORIGIN') for el in re.split(r'\s+//', tmp_lines) if el.strip() != ''][0]

    features = [el.replace('\n', '') for el in re.split('\n\s{5,5}', ''.join(sequences[0].split(
                'Location/Qualifiers')[1].split('BASE COUNT')[0].strip().replace(' '*21, '__')))]

    current_sequence = re.sub('[\d\s]', '', sequences[1]).upper()

    # Read the target motifs

    input_rebase_output = read_rebase_output(output_rebase)
    input_patterns = [input_rebase_output[0][el]['rec_seq']
                    for el in input_rebase_output[0].keys()]
    max_len = max([len(inp_pat) for inp_pat in input_patterns])

    patterns_found_plasmid = findpatterns(
        current_sequence, input_patterns, circular=True)

    positive_strand = current_sequence
    negative_strand = reverse_complement(current_sequence)[::-1]

    feature_type = []
    feature_on_complement = []
    feature_on_join = []

    or_start = []
    or_end = []

    new_start = []
    new_end = []

    data = {'Position_index': list(range(1, len(current_sequence) + 1)),
            "Positive_strand[5'-3']": [x for x in re.split(r'(\w{1})', positive_strand) if x != ''],
            "Negative_strand[3'-5']": [x for x in re.split(r'(\w{1})', negative_strand) if x != '']}
    
    for key1 in input_rebase_output[0].keys():
        occ = 1
        key2 = input_rebase_output[0][key1]['rec_seq']
        if key2 in patterns_found_plasmid.keys():
            for el in patterns_found_plasmid[key2]:
                data['rec_seq:' + input_rebase_output[0][key1]['rec_seq'] + '|meth_base:' + str(input_rebase_output[0][key1]['meth_base']) + '|meth_type:' + str(input_rebase_output[0]
                                                                                                                                                                [key1]['meth_type']) + '|comp_meth_base:' + str(input_rebase_output[0][key1]['comp_meth_base']) + '|comp_meth_type:' + str(input_rebase_output[0][key1]['comp_meth_type']) + '|occurrence:' + str(occ) + '|strand' + str(el[-1]) + '|range:' + str([el[1]+1, el[2]])] = ['' for i in range(el[1])] + [x for x in re.split(r'(\w{1})', key2) if x != ''] + ['' for i in range(el[2], len(current_sequence))]
                occ += 1

    for feat in features:
        if not feat.startswith('source'):
            tmp_feat = feat.split('__')[0]

            indeces = [int(el) for el in re.sub(
                r'\D', ' ', tmp_feat).strip().split(' ') if el != '']

            tmp_start = indeces[0] - 1
            tmp_end = indeces[-1]

            if not 'join' in tmp_feat:
                if 'complement' in tmp_feat:
                    data[feat] = ['' for i in range(
                        tmp_start)] + [reverse_complement(current_sequence)[::-1][i] for i in range(tmp_start, tmp_end)] + ['' for i in range(
                            tmp_end, len(current_sequence))]
                    
                    if 'CDS' in tmp_feat:
                        codons = ['' for i in range(tmp_start)] + sum([[x, x, x] for x in re.split(r'(\w{3})', reverse_complement(
                            current_sequence[tmp_start:tmp_end]))[::-1] if x != ''], []) + ['' for i in range(tmp_end, len(current_sequence))]
                        new_index = ['' for i in range(tmp_start)] + sum([[el, el, el]
                                                                        for el in range(tmp_start, tmp_end, 3)], []) + ['' for i in range(tmp_end, len(current_sequence))]

                        data[feat + '[Codons_representation]'] = codons
                        data[feat + '[Index_representation]'] = new_index
                else:
                    data[feat] = ['' for i in range(
                        tmp_start)] + [current_sequence[i] for i in range(tmp_start, tmp_end)] + ['' for i in range(
                            tmp_end, len(current_sequence))]

                    if 'CDS' in tmp_feat:
                        codons = ['' for i in range(tmp_start)] + sum([[x, x, x] for x in re.split(r'(\w{3})', current_sequence[tmp_start:tmp_end]) if x != ''], []) + ['' for i in range(tmp_end, len(current_sequence))]
                        new_index = ['' for i in range(tmp_start)] + sum([[el, el, el]
                                                                        for el in range(tmp_start, tmp_end, 3)], []) + ['' for i in range(tmp_end, len(current_sequence))]

                        data[feat + '[Codons_representation]'] = codons
                        data[feat + '[Index_representation]'] = new_index
            else:
                if 'complement' in tmp_feat:
                    data[feat] = [reverse_complement(current_sequence)[::-1][i] for i in range(
                        tmp_end)] + ['' for i in range(tmp_end, tmp_start)] + [current_sequence[i] for i in range(
                            tmp_start, len(current_sequence))]

                    if 'CDS' in tmp_feat:
                        tmp_codons = sum([[x, x, x] for x in re.split(r'(\w{3})', reverse_complement(''.join([current_sequence[i] for i in range(tmp_start, len(current_sequence))]) + ''.join([current_sequence[i] for i in range(0, tmp_end)]))[::-1]) if x != ''], [])
                        ttt_irange = irange(tmp_start, tmp_end, l=len(current_sequence))
                        ttt_codons = {}
                        for i in range(len(ttt_irange)):
                            ttt_codons[ttt_irange[i]] = tmp_codons[i]
                        codons = [ttt_codons[i] for i in range(tmp_end)] + ['' for i in range(tmp_end, tmp_start)] + [ttt_codons[i] for i in range(tmp_start, len(current_sequence))]
                        new_index = sum([[el, el, el]
                                        for el in irange(tmp_start, tmp_end, by=3, l=len(current_sequence))], [])
                        nnn_index = {}
                        for i in range(len(ttt_irange)):
                            nnn_index[ttt_irange[i]] = new_index[i]
                        new_index = [nnn_index[i] for i in range(
                            tmp_end)] + ['' for i in range(tmp_end, tmp_start)] + [nnn_index[i] for i in range(
                                tmp_start, len(current_sequence))]
                    
                        data[feat + '[Codons_representation]'] = codons
                        data[feat + '[Index_representation]'] = new_index
                    
                else:
                    data[feat] = [current_sequence[i] for i in range(tmp_end)] + ['' for i in range(tmp_end, tmp_start)] + [current_sequence[i] for i in range(tmp_start, len(current_sequence))]
                    
                    if 'CDS' in tmp_feat:
                        tmp_codons = sum(
                            [
                                [x, x, x] 
                                for x in 
                                re.split(
                                    r'(\w{3})', ''.join([current_sequence[i] for i in range(tmp_start, len(current_sequence))]) + ''.join([current_sequence[i] for i in range(0, tmp_end)])
                                    )
                                if x != ''
                                ], [])
                        ttt_irange = irange(tmp_start, tmp_end, l=len(current_sequence))
                        ttt_codons = {}
                        for i in range(len(ttt_irange)):
                            ttt_codons[ttt_irange[i]] = tmp_codons[i]
                        codons = [ttt_codons[i] for i in range(tmp_end)] + ['' for i in range(tmp_end, tmp_start)] + [ttt_codons[i] for i in range(tmp_start, len(current_sequence))]
                        new_index = sum([[el, el, el]
                                        for el in irange(tmp_start, tmp_end, by=3, l=len(current_sequence))], [])
                        nnn_index = {}
                        for i in range(len(ttt_irange)):
                            nnn_index[ttt_irange[i]] = new_index[i]
                        new_index = [nnn_index[i] for i in range(
                            tmp_end)] + ['' for i in range(tmp_end, tmp_start)] + [nnn_index[i] for i in range(
                                tmp_start, len(current_sequence))]

                        data[feat + '[Codons_representation]'] = codons
                        data[feat + '[Index_representation]'] = new_index
                    

    stats = pd.DataFrame(data=data)

    stats.to_excel(os.path.join(output_folder, 'input_sequence_representation.xlsx'))

    return


def run_motifs_finder(args):

    # Params: input_tools [list], input_rm_files [list], verbose [bool], output_folder [str path]

    # SeqFeature(
    #   CompoundLocation(
    #       [
    #           FeatureLocation(ExactPosition(4373), ExactPosition(4574), strand=1), 
    #           FeatureLocation(ExactPosition(0), ExactPosition(1275), strand=1)
    #       ], 
    #       'join'
    #   ), 
    #   type='CDS', 
    #   location_operator='join'
    # )    

    color_map = {}

    # annotations = ["CDS", "RMM"]
    annotations = ["CDS"]

    class ExpressionUnitTranslator(BiopythonTranslator):
        def compute_feature_color(self, feature):            
            return color_map[feature.type]
        def compute_feature_label(self, feature):
            if feature.type not in annotations:
                return None
            else:
                return None

    input_tools = args.input_tools
    input_rm_files = args.input_rm_files

    name_of_the_project = args.project_name
    width_of_the_image = args.image_width

    verbose = args.verbose
    output_folder = args.output_folder

    pairs = list(itertools.product(input_tools, input_rm_files))

    if not len(pairs) >= 1:
        error('No input files provided!')
    
    out_dirs = [name_of_the_project + "-"+ "0"*(len(str(len(pairs)))-len(str(i)))+str(i) for i in range(1, len(pairs)+1)]

    identifiers = ''

    for pair in pairs:

        annotations = ["CDS"]

        curr_identifier = out_dirs.pop(0)

        tt_name1 = os.path.basename(pair[0]).split("__input_genetic_tool_file-")[-1]
        tt_name2 = os.path.basename(pair[1]).split("__input_rebase_file-")[-1]

        identifiers += f"Id: {curr_identifier}, \n\tGenetic tool: {tt_name1}, \n\tRM motifs: {tt_name2}\n"

        out_dir = os.path.join(output_folder, curr_identifier)

        try:
            os.mkdir(out_dir)
        except OSError:
            pass
        else:
            pass

        input_tool = pair[0]
        input_rm_file = pair[1]

        args = SimpleNamespace(input_sequence=input_tool, input_rm_systems=input_rm_file, output_folder=out_dir)
        run_tool_representation(args)

        out_tool_representation = pd.read_excel(os.path.join(out_dir, 'input_sequence_representation.xlsx'),engine='openpyxl')

        rec_mots = [el for el in out_tool_representation.columns if el.startswith('rec_seq:')]
        cds_s = [el for el in out_tool_representation.columns if ((el.startswith('CDS')) & (not el.endswith('[Index_representation]')) & (not el.endswith('[Codons_representation]')))]

        out_table = {
            'RM motifs': [],
            'non-CDS': []
        }

        # read basic

        input_plasmid = read_input_seq(input_tool)

        if not input_plasmid['correct']:
            info('The input plasmid is not correct! Please, check it.')
            error('\n'.join(input_plasmid['messages']))
        else:
            if len(input_plasmid['annotations'].keys()) == 1:
                if len(input_plasmid['annotations'][0]['sequence']) > 0:
                    plasmid_sequence = input_plasmid['annotations'][0]['sequence']
                    plasmid_annotations = input_plasmid['annotations'][0]['cds']
                    l = len(plasmid_sequence)
                else:
                    error('More sequences in the input plasmid file OR the sequence of the input file is empty!')
            else:
                error('Incorrect number of inputs for the plasmid sequence!')

        input_rebase_output = read_rebase_output(input_rm_file)
        input_patterns = [input_rebase_output[0][el]['rec_seq']
                        for el in input_rebase_output[0].keys()]
        max_len = max([len(inp_pat) for inp_pat in input_patterns])

        # end of read basic

        if ((len(rec_mots) > 0) & (len(cds_s) > 0)):
            for rec_mot in rec_mots:
                out_table['RM motifs'].append(rec_mot)
                
                subset_cols = [rec_mot] + cds_s
                tmp_out_tool_representation = copy.copy(out_tool_representation[subset_cols])
                nn_tmp_out_tool_representation = tmp_out_tool_representation[(tmp_out_tool_representation[rec_mot].notna())]

                for cds in cds_s:
                    
                    nn_tmp_out_tool_representation = nn_tmp_out_tool_representation[~(nn_tmp_out_tool_representation[cds].notna())]

                    if not cds in out_table.keys():
                        out_table[cds] = []
                    out_table[cds].append(out_tool_representation[[rec_mot, cds]].dropna().shape[0])

                out_table['non-CDS'].append(nn_tmp_out_tool_representation.shape[0])

            n_out_table = pd.DataFrame.from_dict(out_table)

            t_out_table = copy.copy(n_out_table)
            t_out_table.set_index('RM motifs', inplace=True)
            t_out_table[t_out_table > 0] = 1

            tt_out_table = copy.copy(t_out_table)

            tt_out_table['Recognition motif'] = ''
            tt_out_table['Methylation base (+ strand)'] = ''
            tt_out_table['Methylation type (+ strand)'] = ''
            tt_out_table['Methylation base (- strand)'] = ''
            tt_out_table['Methylation type (- strand)'] = ''
            tt_out_table['Occurrence range'] = ''
            tt_out_table['Number of occurrences in CDS'] = ''
            tt_out_table['Number of occurrences in non-CDS'] = ''

            for el in tt_out_table.index:
                ttt_ell = [s for s in el.split("|") if s != '']
                
                cccc_rec_seq = ttt_ell[0].replace('rec_seq:', '')
                cccc_meth_base = ttt_ell[1].replace('meth_base:', '') if ttt_ell[1].replace('meth_base:', '')!= '-inf' else ''
                cccc_meth_type = ttt_ell[2].replace('meth_type:', '') if ttt_ell[2].replace('meth_type:', '')!= '-inf' else ''
                cccc_comp_meth_base = ttt_ell[3].replace('comp_meth_base:', '') if ttt_ell[3].replace('comp_meth_base:', '')!= '-inf' else ''
                cccc_comp_meth_type = ttt_ell[4].replace('comp_meth_type:', '') if ttt_ell[4].replace('comp_meth_type:', '')!= '-inf' else ''

                tt_out_table.loc[el, 'Recognition motif'] = cccc_rec_seq

                if cccc_meth_base.isdigit():
                    if int(cccc_meth_base) in range(len(cccc_rec_seq)):
                        tt_out_table.loc[el, 'Methylation base (+ strand)'] = f"{cccc_rec_seq[:int(cccc_meth_base)]}*{cccc_rec_seq[int(cccc_meth_base)]}*{cccc_rec_seq[int(cccc_meth_base):]}"
                    else:
                        tt_out_table.loc[el, 'Methylation base (+ strand)'] = '-'
                else:
                    tt_out_table.loc[el, 'Methylation base (+ strand)'] = '-'

                if cccc_meth_type in ['6', '4', '5']:
                    tt_out_table.loc[el, 'Methylation type (+ strand)'] = 'm6A' if cccc_meth_type == '6' else 'm4C' if cccc_meth_type == '4' else 'm5C'
                else:
                    tt_out_table.loc[el, 'Methylation type (+ strand)'] = '-'

                if cccc_comp_meth_base.isdigit():
                    if int(cccc_comp_meth_base) in range(len(cccc_rec_seq)):
                        tt_out_table.loc[el, 'Methylation base (- strand)'] = f"{cccc_rec_seq[:int(cccc_comp_meth_base)]}*{cccc_rec_seq[int(cccc_comp_meth_base)]}*{cccc_rec_seq[int(cccc_comp_meth_base):]}"
                    else:
                        tt_out_table.loc[el, 'Methylation base (- strand)'] = '-'
                else:
                    tt_out_table.loc[el, 'Methylation base (- strand)'] = '-'

                if cccc_comp_meth_type in ['6', '4', '5']:
                    tt_out_table.loc[el, 'Methylation type (- strand)'] = 'm6A' if cccc_comp_meth_type == '6' else 'm4C' if cccc_comp_meth_type == '4' else 'm5C'
                else:
                    tt_out_table.loc[el, 'Methylation type (- strand)'] = '-'
                
                tt_out_table.loc[el, 'Occurrence range'] = ttt_ell[7].replace('range:', '')
            
            unmapped_rm_motifs = set(input_patterns) - set(tt_out_table['Recognition motif'].values)

            tmp_dat_to_add = {
                "Recognition motif": [],
                "Methylation base (+ strand)": [],
                "Methylation type (+ strand)": [],
                "Methylation base (- strand)": [],
                "Methylation type (- strand)": [],
                "Occurrence range": [],
                "Number of occurrences in CDS": [],
                "Number of occurrences in non-CDS": []
            }

            for key in input_rebase_output[0].keys():
                if input_rebase_output[0][key]['rec_seq'] in unmapped_rm_motifs:
                    tmp_dat_to_add["Recognition motif"].append(input_rebase_output[0][key]['rec_seq'])
                    tmp_dat_to_add["Methylation base (+ strand)"].append(input_rebase_output[0][key]['meth_base'])
                    tmp_dat_to_add["Methylation type (+ strand)"].append(input_rebase_output[0][key]['meth_type'])
                    tmp_dat_to_add["Methylation base (- strand)"].append(input_rebase_output[0][key]['comp_meth_base'])
                    tmp_dat_to_add["Methylation type (- strand)"].append(input_rebase_output[0][key]['comp_meth_type'])
                    tmp_dat_to_add["Occurrence range"].append("")
                    tmp_dat_to_add["Number of occurrences in CDS"].append(0)
                    tmp_dat_to_add["Number of occurrences in non-CDS"].append(0)


            tmp_dat_to_add_dataframe = pd.DataFrame.from_dict(tmp_dat_to_add)

            tt_out_table = tt_out_table.append(tmp_dat_to_add_dataframe, ignore_index=True)
            
            subset_cols = ['Recognition motif', 'Methylation base (+ strand)', 'Methylation type (+ strand)', 'Methylation base (- strand)', 'Methylation type (- strand)', 'Occurrence range']

            tt_out_table.drop_duplicates(subset=subset_cols, inplace=True)
            tt_out_table['Number of occurrences in CDS'] = tt_out_table[[el for el in tt_out_table.columns if el.startswith('CDS')]].sum(axis=1)
            cols_to_remove = [el for el in tt_out_table.columns if el.startswith('CDS')]
            tt_out_table.drop(cols_to_remove, axis=1, inplace=True)
            tt_out_table['Occurrence range'] = [1 if len(el)>0 else 0 for el in tt_out_table['Occurrence range'].values]
            subset_cols = ['Recognition motif', 'Methylation base (+ strand)', 'Methylation type (+ strand)', 'Methylation base (- strand)', 'Methylation type (- strand)']
            tt_out_table = tt_out_table.groupby(subset_cols)[['Number of occurrences in CDS', 'non-CDS', 'Occurrence range']].sum().reset_index()

            tt_out_table.rename(columns={'Occurrence range': 'Total # motifs found in the sequence', 'Number of occurrences in CDS': '# of motifs found in coding regions', 'non-CDS': '# of motifs found in non-coding regions'}, inplace=True)

            mm_out_table = copy.copy(tt_out_table)
            mm_out_table.drop(['Methylation base (+ strand)', 'Methylation type (+ strand)', 'Methylation base (- strand)', 'Methylation type (- strand)'], axis=1, inplace=True)

            with pd.ExcelWriter(os.path.join(out_dir, 'RM_statistics.xlsx')) as writer: 
                mm_out_table.to_excel(writer, sheet_name='Main', index=False) 
                tt_out_table.to_excel(writer, sheet_name='Main_spec', index=False)
                t_out_table.to_excel(writer, sheet_name='Aggregated')
                n_out_table.to_excel(writer, sheet_name='Overlapping_bases')

        patterns_found_plasmid = findpatterns(plasmid_sequence, input_patterns, circular=True)
        
        if not len(sum(list(patterns_found_plasmid.values()), [])) > 0:
            info("")
            info(f"The number of target motifs found on the plasmid is: {len(sum(list(patterns_found_plasmid.values()), []))}")
        else:
            sequence = list(SeqIO.parse(input_tool, "genbank"))

            if not len(sequence) == 1:
                error(f'Multiple sequence in the input genetic tool {input_tool}!')
            else:
                sequence = sequence[0]
            tmp_feat_dict = {}
            i = 0
            for key in patterns_found_plasmid.keys():
                targets = patterns_found_plasmid[key]

                tmp_feature_type = f"RMM-{i}"
                annotations.append(tmp_feature_type)
                tmp_feat_dict[tmp_feature_type] = key
                i += 1

                for tar in targets:
                    or_patt = key
                    found_seq = tar[0]
                    or_pos = (tar[1] - 10) % l
                    en_pos = (tar[2] + 10) % l
                    strand = tar[3]

                    if or_pos < en_pos:
                        tmp_start_pos = Bio.SeqFeature.ExactPosition(or_pos)
                        tmp_end_pos = Bio.SeqFeature.ExactPosition(en_pos)                
                        tmp_feature_location = FeatureLocation(tmp_start_pos, tmp_end_pos, strand = strand)
                        my_feature = SeqFeature(tmp_feature_location, type = tmp_feature_type, qualifiers={
                            'note': key#, 'label':key
                            })
                        sequence.features.append(my_feature)
                    else:
                        tmp_start_pos1 = Bio.SeqFeature.ExactPosition(or_pos)
                        tmp_end_pos1 = Bio.SeqFeature.ExactPosition(len(sequence.seq))
                        tmp_start_pos2 = Bio.SeqFeature.ExactPosition(0)
                        tmp_end_pos2 = Bio.SeqFeature.ExactPosition(en_pos)
                        l1 = Bio.SeqFeature.FeatureLocation(tmp_start_pos1, tmp_end_pos1, strand=strand)
                        l2 = Bio.SeqFeature.FeatureLocation(tmp_start_pos2, tmp_end_pos2, strand=strand)
                        
                        my_feature_l1 = SeqFeature(l1, type = tmp_feature_type, qualifiers={
                            'note': key#, 'label':key
                            })
                        sequence.features.append(my_feature_l1)
                        
                        my_feature_l2 = SeqFeature(l2, type = tmp_feature_type, qualifiers={
                            'note': key#, 'label':key
                            })
                        sequence.features.append(my_feature_l2)

            tmp_new_features = []

            for feat in sequence.features:
                if not "join" in str(feat.location):
                    tmp_new_features.append(feat)
                else:
                    ttt_feat1 = copy.copy(feat)
                    ttt_feat1.location = feat.location.parts[0]
                    ttt_feat2 = copy.copy(feat)
                    ttt_feat2.location = feat.location.parts[1]

                    tmp_new_features.append(ttt_feat1)
                    tmp_new_features.append(ttt_feat2)

            sequence.features = tmp_new_features

            with open(os.path.join(out_dir, f"{curr_identifier}_gen_tool.gb"), "w") as output_handle:
                count = SeqIO.write(sequence, output_handle, "genbank")
            
            feature_types = set([feat.type for feat in sequence.features])

            colors = set(["#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a"])

            while len(colors) < len(feature_types):
                colors |= set(["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])])
            
            colors = colors - set(["#000000", "#a6cee3"])

            for feat in feature_types:
                if not "source" in feat:
                    if "cds" in feat.lower():
                        color_map[feat] = "#a6cee3"
                    else:
                        color_map[feat] = colors.pop()
                else:
                    color_map[feat] = "#000000"

            translator = ExpressionUnitTranslator()
            graphic_record = translator.translate_record(os.path.join(out_dir, f"{curr_identifier}_gen_tool.gb"), record_class=CircularGraphicRecord)
            graphic_record.top_position = 4800

            try:
                ax1, _ = graphic_record.plot(figure_width=width_of_the_image)
                patches = []
                for key_nn in tmp_feat_dict.keys():
                    if key_nn in color_map.keys():
                        patches.append(mpatches.Patch(color=color_map[key_nn], label=tmp_feat_dict[key_nn]))
                patches.append(mpatches.Patch(color=color_map["CDS"], label="CDS"))
                ax1.legend(handles=patches)
                ax1.set_title(f"{tt_name1}\n{tt_name2}", fontsize=14)
                ax1.figure.savefig(os.path.join(out_dir, f"{curr_identifier}_circ_tool_repr.png"), bbox_inches="tight", dpi=400)
                ax1.figure.savefig(os.path.join(out_dir, f"{curr_identifier}_circ_tool_repr.svg"), bbox_inches="tight", format="svg")
            except:          
                info(f'Error during plot! Have a look to the annotated genetic tool at {os.path.join(out_dir, f"{curr_identifier}_gen_tool.gb")}.')

            annotations = ["CDS"]

            class ExpressionUnitTranslator(BiopythonTranslator):
                def compute_feature_color(self, feature):            
                    return color_map[feature.type]
                def compute_feature_label(self, feature):
                    if feature.type not in annotations:
                        return None
                    else:
                        return BiopythonTranslator.compute_feature_label(self, feature)

            translator = ExpressionUnitTranslator()
            graphic_record = translator.translate_record(os.path.join(out_dir, f"{curr_identifier}_gen_tool.gb"))
            graphic_record.top_position = 4800

            try:
                ax2, _ = graphic_record.plot(figure_width=width_of_the_image)
                patches = []
                for key_nn in tmp_feat_dict.keys():
                    if key_nn in color_map.keys():
                        patches.append(mpatches.Patch(color=color_map[key_nn], label=tmp_feat_dict[key_nn]))
                patches.append(mpatches.Patch(color=color_map["CDS"], label="CDS"))
                ax2.legend(handles=patches, loc='upper center', bbox_to_anchor=(1.25, 0.8))
                ax2.set_title(f"{tt_name1}\n{tt_name2}", fontsize=14)
                ax2.figure.savefig(os.path.join(out_dir, f"{curr_identifier}_lin_tool_repr.png"), bbox_inches="tight", dpi=400)
                ax2.figure.savefig(os.path.join(out_dir, f"{curr_identifier}_lin_tool_repr.svg"), bbox_inches="tight", format="svg")
            except:
                info(f'Error during plot! Have a look to the annotated genetic tool at {os.path.join(out_dir, f"{curr_identifier}_gen_tool.gb")}.')

            if ((ax1 is not None) & (ax2 is not None)):

                pdf = matplotlib.backends.backend_pdf.PdfPages(os.path.join(out_dir, f"{curr_identifier}_combined_plot_tool_repr.pdf"))
                pdf.savefig(ax1.figure, bbox_inches="tight")
                pdf.savefig(ax2.figure, bbox_inches="tight")
                pdf.close()

    with open(os.path.join(output_folder, f"{name_of_the_project}_map.txt"), 'w') as handle:
        handle.write(identifiers.strip())

    return

def run_motifs_finder2(args):


    # Params: input_tools [list], input_rm_files [list], verbose [bool], output_folder [str path]

    # SeqFeature(
    #   CompoundLocation(
    #       [
    #           FeatureLocation(ExactPosition(4373), ExactPosition(4574), strand=1),
    #           FeatureLocation(ExactPosition(0), ExactPosition(1275), strand=1)
    #       ], 
    #       'join'
    #   ), 
    #   type='CDS', 
    #   location_operator='join'
    # )    

    color_map = {}

    # annotations = ["CDS", "RMM"]
    annotations = ["CDS"]

    class ExpressionUnitTranslator(BiopythonTranslator):
        def compute_feature_color(self, feature):            
            return color_map[feature.type]
        def compute_feature_label(self, feature):
            if feature.type not in annotations:
                return None
            else:
                return None

    input_tools = args.input_tools
    input_rm_files = args.input_rm_files

    name_of_the_project = args.project_name
    width_of_the_image = args.image_width

    verbose = args.verbose
    output_folder = args.output_folder

    pairs = list(itertools.product(input_tools, input_rm_files))

    if not len(pairs) >= 1:
        error('No input files provided!')
    
    out_dirs = [name_of_the_project + "-"+ "0"*(len(str(len(pairs)))-len(str(i)))+str(i) for i in range(1, len(pairs)+1)]

    identifiers = ''

    for pair in pairs:

        annotations = ["CDS"]

        curr_identifier = out_dirs.pop(0)

        tt_name1 = os.path.basename(pair[0]).split("__input_genetic_tool_file-")[-1]
        tt_name2 = os.path.basename(pair[1]).split("__input_rebase_file-")[-1]

        identifiers += f"Id: {curr_identifier}, \n\tGenetic tool: {tt_name1}, \n\tRM motifs: {tt_name2}\n"

        out_dir = os.path.join(output_folder, curr_identifier)

        try:
            os.mkdir(out_dir)
        except OSError:
            pass
        else:
            pass

        input_tool = pair[0]
        input_rm_file = pair[1]

        # args = SimpleNamespace(input_sequence=input_tool, input_rm_systems=input_rm_file, output_folder=out_dir)
        # run_tool_representation(args)

        # out_tool_representation = pd.read_excel(os.path.join(out_dir, 'input_sequence_representation.xlsx'),engine='openpyxl')

        # rec_mots = [el for el in out_tool_representation.columns if el.startswith('rec_seq:')]
        # cds_s = [el for el in out_tool_representation.columns if ((el.startswith('CDS')) & (not el.endswith('[Index_representation]')) & (not el.endswith('[Codons_representation]')))]

        # out_table = {
        #     'RM motifs': [],
        #     'non-CDS': []
        # }

        # read basic

        input_plasmid = read_input_seq(input_tool)

        if not input_plasmid['correct']:
            info('The input plasmid is not correct! Please, check it.')
            error('\n'.join(input_plasmid['messages']))
        else:
            if len(input_plasmid['annotations'].keys()) == 1:
                if len(input_plasmid['annotations'][0]['sequence']) > 0:
                    plasmid_sequence = input_plasmid['annotations'][0]['sequence']
                    plasmid_annotations = input_plasmid['annotations'][0]['cds']
                    l = len(plasmid_sequence)
                else:
                    error('More sequences in the input plasmid file OR the sequence of the input file is empty!')
            else:
                error('Incorrect number of inputs for the plasmid sequence!')
        
        input_rebase_output = read_rebase_output(input_rm_file)
        input_patterns = [input_rebase_output[0][el]['rec_seq']
                        for el in input_rebase_output[0].keys()]

        if len(input_patterns) == 0:
            return
        max_len = max([len(inp_pat) for inp_pat in input_patterns])

        patterns_found_plasmid = findpatterns(plasmid_sequence, input_patterns, circular=True)

        out_table = {
            "Recognition motif": [],
            "# of motifs found in coding regions": [],
            "# of motifs found in non-coding regions": [],
            "Total # motifs found in the sequence": []
        }

        for rec_seq in patterns_found_plasmid.keys():
            n_occ = len(set([(el[1], el[2]) for el in patterns_found_plasmid[rec_seq]]))
            n_occ_in_cds = 0

            for cds in [set(irange(el[0], el[1], l=len(plasmid_sequence))) for el in plasmid_annotations]:
                n_occ_in_cds += sum([1 for el in set([(el1[1], el1[2]) for el1 in patterns_found_plasmid[rec_seq]]) if len(set(irange(el[0], el[1], l=len(plasmid_sequence))) & cds) > 0])

            n_occ_out_cds = sum([1 for el in patterns_found_plasmid[rec_seq] if not len(set(irange(el[1], el[2], l=len(plasmid_sequence))) & set.union(*[set(irange(el[0], el[1], l=len(plasmid_sequence))) for el in plasmid_annotations])) > 0])

            out_table["Recognition motif"].append(rec_seq)
            out_table["# of motifs found in coding regions"].append(n_occ_in_cds)
            out_table["# of motifs found in non-coding regions"].append(n_occ_out_cds)
            out_table["Total # motifs found in the sequence"].append(n_occ)

        
        out_table = pd.DataFrame.from_dict(out_table)

        s_out_table = {
            "Recognition motif": [],
            "Methylation base (+ strand)": [],
            "Methylation type (+ strand)": [],
            "Methylation base (- strand)": [],
            "Methylation type (- strand)": []
        }

        for el in input_rebase_output[0]:

            cccc_rec_seq = input_rebase_output[0][el]['rec_seq']
            cccc_meth_base = input_rebase_output[0][el]['meth_base']
            cccc_meth_type = input_rebase_output[0][el]['meth_type']
            cccc_comp_meth_base = input_rebase_output[0][el]['comp_meth_base']
            cccc_comp_meth_type = input_rebase_output[0][el]['comp_meth_type']
            
            s_out_table["Recognition motif"].append(cccc_rec_seq)

            if str(cccc_meth_base).isdigit():
                if int(cccc_meth_base) in range(len(cccc_rec_seq)):
                    s_out_table['Methylation base (+ strand)'].append(f"{cccc_rec_seq[:int(cccc_meth_base)]}*{cccc_rec_seq[int(cccc_meth_base)]}*{cccc_rec_seq[int(cccc_meth_base):]}")
                else:
                    s_out_table['Methylation base (+ strand)'].append('-')
            else:
                s_out_table['Methylation base (+ strand)'].append('-')

            if cccc_meth_type in ['6', '4', '5']:
                s_out_table['Methylation type (+ strand)'].append('m6A' if cccc_meth_type == '6' else 'm4C' if cccc_meth_type == '4' else 'm5C')
            else:
                s_out_table['Methylation type (+ strand)'].append('-')

            if str(cccc_comp_meth_base).isdigit():
                if int(cccc_comp_meth_base) in range(len(cccc_rec_seq)):
                    s_out_table['Methylation base (- strand)'].append(f"{cccc_rec_seq[:int(cccc_comp_meth_base)]}*{cccc_rec_seq[int(cccc_comp_meth_base)]}*{cccc_rec_seq[int(cccc_comp_meth_base):]}")
                else:
                    s_out_table['Methylation base (- strand)'].append('-')
            else:
                s_out_table['Methylation base (- strand)'].append('-')

            if cccc_comp_meth_type in ['6', '4', '5']:
                s_out_table['Methylation type (- strand)'].append('m6A' if cccc_comp_meth_type == '6' else 'm4C' if cccc_comp_meth_type == '4' else 'm5C')
            else:
                s_out_table['Methylation type (- strand)'].append('-')

        s_out_table = pd.DataFrame.from_dict(s_out_table)

        with pd.ExcelWriter(os.path.join(out_dir, 'RM_statistics.xlsx')) as writer: 
            out_table.to_excel(writer, sheet_name='Main', index=False) 
            s_out_table.to_excel(writer, sheet_name='Main_spec', index=False)
        
        if not len(sum(list(patterns_found_plasmid.values()), [])) > 0:
            info("")
            info(f"The number of target motifs found on the plasmid is: {len(sum(list(patterns_found_plasmid.values()), []))}")
            
            with open(os.path.join(out_dir, f"{curr_identifier}_n_motifs.txt"), "w") as output_handle:
                output_handle.write(f"The number of target motifs found on the plasmid is: {len(sum(list(patterns_found_plasmid.values()), []))}")
            
        else:
            sequence = list(SeqIO.parse(input_tool, "genbank"))
            sequence_repr = list(SeqIO.parse(input_tool, "genbank"))

            if not len(sequence) == 1:
                error(f'Multiple sequence in the input genetic tool {input_tool}!')
            else:
                sequence = sequence[0]
                sequence_repr = sequence_repr[0]
            tmp_feat_dict = {}
            i = 0

            for key in patterns_found_plasmid.keys():

                print("\n", key)
                print(sorted(patterns_found_plasmid[key]))

                new_targets = {}

                for tar in patterns_found_plasmid[key]:
                    if not (tar[0], tar[1], tar[2]) in new_targets.keys(): 
                        new_targets[(tar[0], tar[1], tar[2])] = 0
                    
                    new_targets[(tar[0], tar[1], tar[2])] += tar[3]

                print(new_targets)

                nn_targets = []

                for tar in new_targets.keys():
                    if new_targets[tar] < 0:
                        new_targets[tar] = -1
                    elif new_targets[tar] > 0:
                        new_targets[tar] = +1

                    nn_targets.append((tar[0], tar[1], tar[2], new_targets[tar]))
                
                # targets = patterns_found_plasmid[key]
                targets = sorted(nn_targets)

                tmp_feature_type = f"RMM-{i}"
                annotations.append(tmp_feature_type)
                tmp_feat_dict[tmp_feature_type] = key
                i += 1

                for tar in targets:
                    or_patt = key
                    found_seq = tar[0]
                    or_pos = tar[1]
                    en_pos = tar[2]
                    strand = tar[3]

                    if or_pos < en_pos:
                        tmp_start_pos = Bio.SeqFeature.ExactPosition(or_pos)
                        tmp_end_pos = Bio.SeqFeature.ExactPosition(en_pos)                
                        tmp_feature_location = FeatureLocation(tmp_start_pos, tmp_end_pos, strand = strand)
                        my_feature = SeqFeature(tmp_feature_location, type = tmp_feature_type, qualifiers={
                            'note': key#, 'label':key
                            })
                        sequence.features.append(my_feature)
                    else:
                        tmp_start_pos1 = Bio.SeqFeature.ExactPosition(or_pos)
                        tmp_end_pos1 = Bio.SeqFeature.ExactPosition(len(sequence.seq))
                        tmp_start_pos2 = Bio.SeqFeature.ExactPosition(0)
                        tmp_end_pos2 = Bio.SeqFeature.ExactPosition(en_pos)
                        l1 = Bio.SeqFeature.FeatureLocation(tmp_start_pos1, tmp_end_pos1, strand=strand)
                        l2 = Bio.SeqFeature.FeatureLocation(tmp_start_pos2, tmp_end_pos2, strand=strand)
                        
                        my_feature_l1 = SeqFeature(l1, type = tmp_feature_type, qualifiers={
                            'note': key#, 'label':key
                            })
                        sequence.features.append(my_feature_l1)
                        
                        my_feature_l2 = SeqFeature(l2, type = tmp_feature_type, qualifiers={
                            'note': key#, 'label':key
                            })
                        sequence.features.append(my_feature_l2)

                    # del or_pos
                    # del en_pos
                    # del tmp_start_pos
                    # del tmp_end_pos
                    # del tmp_feature_location
                    # del my_feature
                    # del tmp_start_pos1
                    # del tmp_end_pos1
                    # del tmp_start_pos2
                    # del tmp_end_pos2
                    # del l1
                    # del l2
                    # del my_feature_l1
                    # del my_feature_l2

                    or_pos = (tar[1] - 10) % l
                    en_pos = (tar[2] + 10) % l
                    
                    if or_pos < en_pos:
                        tmp_start_pos = Bio.SeqFeature.ExactPosition(or_pos)
                        tmp_end_pos = Bio.SeqFeature.ExactPosition(en_pos)                
                        tmp_feature_location = FeatureLocation(tmp_start_pos, tmp_end_pos, strand = strand)
                        my_feature = SeqFeature(tmp_feature_location, type = tmp_feature_type, qualifiers={
                            'note': key#, 'label':key
                            })
                        sequence_repr.features.append(my_feature)
                    else:
                        tmp_start_pos1 = Bio.SeqFeature.ExactPosition(or_pos)
                        tmp_end_pos1 = Bio.SeqFeature.ExactPosition(len(sequence_repr.seq))
                        tmp_start_pos2 = Bio.SeqFeature.ExactPosition(0)
                        tmp_end_pos2 = Bio.SeqFeature.ExactPosition(en_pos)
                        l1 = Bio.SeqFeature.FeatureLocation(tmp_start_pos1, tmp_end_pos1, strand=strand)
                        l2 = Bio.SeqFeature.FeatureLocation(tmp_start_pos2, tmp_end_pos2, strand=strand)
                        
                        my_feature_l1 = SeqFeature(l1, type = tmp_feature_type, qualifiers={
                            'note': key#, 'label':key
                            })
                        sequence_repr.features.append(my_feature_l1)
                        
                        my_feature_l2 = SeqFeature(l2, type = tmp_feature_type, qualifiers={
                            'note': key#, 'label':key
                            })
                        sequence_repr.features.append(my_feature_l2)

                    

            tmp_new_features = []

            for feat in sequence.features:
                if not "join" in str(feat.location):
                    tmp_new_features.append(feat)
                else:
                    ttt_feat1 = copy.copy(feat)
                    ttt_feat1.location = feat.location.parts[0]
                    ttt_feat2 = copy.copy(feat)
                    ttt_feat2.location = feat.location.parts[1]

                    tmp_new_features.append(ttt_feat1)
                    tmp_new_features.append(ttt_feat2)

            sequence.features = tmp_new_features

            tmp_new_features = []

            for feat in sequence_repr.features:
                if not "join" in str(feat.location):
                    tmp_new_features.append(feat)
                else:
                    ttt_feat1 = copy.copy(feat)
                    ttt_feat1.location = feat.location.parts[0]
                    ttt_feat2 = copy.copy(feat)
                    ttt_feat2.location = feat.location.parts[1]

                    tmp_new_features.append(ttt_feat1)
                    tmp_new_features.append(ttt_feat2)

            sequence_repr.features = tmp_new_features

            with open(os.path.join(out_dir, f"{curr_identifier}_gen_tool.gb"), "w") as output_handle:
                SeqIO.write(sequence, output_handle, "genbank")

            with open(os.path.join(out_dir, f"{curr_identifier}_gen_tool_repr.gb"), "w") as output_handle:
                SeqIO.write(sequence_repr, output_handle, "genbank")

            feature_types = set([feat.type for feat in sequence_repr.features])

            # colors = set(["#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a"])

            # while len(colors) < len(feature_types):
            #     colors |= set(["#"+''.join([random.choice('0123456789ABCDEF') for j in range(6)])])
            
            # colors = colors - set(["#000000", "#a6cee3"])

            # for feat in feature_types:
            #     if not "source" in feat:
            #         if "cds" in feat.lower():
            #             color_map[feat] = "#a6cee3"
            #         else:
            #             color_map[feat] = colors.pop()
            #     else:
            #         color_map[feat] = "#000000"

            # translator = ExpressionUnitTranslator()
            # graphic_record = translator.translate_record(os.path.join(out_dir, f"{curr_identifier}_gen_tool_repr.gb"), record_class=CircularGraphicRecord)
            # graphic_record.top_position = 4800

            # try:
            #     ax1, _ = graphic_record.plot(figure_width=width_of_the_image)
            #     patches = []
            #     for key_nn in tmp_feat_dict.keys():
            #         if key_nn in color_map.keys():
            #             patches.append(mpatches.Patch(color=color_map[key_nn], label=tmp_feat_dict[key_nn]))
            #     patches.append(mpatches.Patch(color=color_map["CDS"], label="CDS"))
            #     ax1.legend(handles=patches)
            #     ax1.set_title(f"{tt_name1}\n{tt_name2}", fontsize=14)
            #     ax1.figure.savefig(os.path.join(out_dir, f"{curr_identifier}_circ_tool_repr.png"), bbox_inches="tight", dpi=400)
            #     ax1.figure.savefig(os.path.join(out_dir, f"{curr_identifier}_circ_tool_repr.svg"), bbox_inches="tight", format="svg")
            # except:          
            #     info(f'Error during plot! Have a look to the annotated genetic tool at {os.path.join(out_dir, f"{curr_identifier}_gen_tool_repr.gb")}.')

            # annotations = ["CDS"]

            # class ExpressionUnitTranslator(BiopythonTranslator):
            #     def compute_feature_color(self, feature):            
            #         return color_map[feature.type]
            #     def compute_feature_label(self, feature):
            #         if feature.type not in annotations:
            #             return None
            #         else:
            #             return BiopythonTranslator.compute_feature_label(self, feature)

            # translator = ExpressionUnitTranslator()
            # graphic_record = translator.translate_record(os.path.join(out_dir, f"{curr_identifier}_gen_tool_repr.gb"))
            # graphic_record.top_position = 4800

            # try:
            #     ax2, _ = graphic_record.plot(figure_width=width_of_the_image)
            #     patches = []
            #     for key_nn in tmp_feat_dict.keys():
            #         if key_nn in color_map.keys():
            #             patches.append(mpatches.Patch(color=color_map[key_nn], label=tmp_feat_dict[key_nn]))
            #     patches.append(mpatches.Patch(color=color_map["CDS"], label="CDS"))
            #     ax2.legend(handles=patches, loc='upper center', bbox_to_anchor=(1.25, 0.8))
            #     ax2.set_title(f"{tt_name1}\n{tt_name2}", fontsize=14)
            #     ax2.figure.savefig(os.path.join(out_dir, f"{curr_identifier}_lin_tool_repr.png"), bbox_inches="tight", dpi=400)
            #     ax2.figure.savefig(os.path.join(out_dir, f"{curr_identifier}_lin_tool_repr.svg"), bbox_inches="tight", dpi="svg")
            # except:
            #     info(f'Error during plot! Have a look to the annotated genetic tool at {os.path.join(out_dir, f"{curr_identifier}_gen_tool_repr.gb")}.')

            # if ((ax1 is not None) & (ax2 is not None)):

            #     pdf = matplotlib.backends.backend_pdf.PdfPages(os.path.join(out_dir, f"{curr_identifier}_combined_plot_tool_repr.pdf"))
            #     pdf.savefig(ax1.figure, bbox_inches="tight")
            #     pdf.savefig(ax2.figure, bbox_inches="tight")
            #     pdf.close()

    with open(os.path.join(output_folder, f"{name_of_the_project}_map.txt"), 'w') as handle:
        handle.write(identifiers.strip())

    return

sytogen_description = 'SyToGen is a computational pipeline for genetic engineering. It includes a set of modules: make_assembly, make_preprocess, mapping_targets, partial_generalization, candidate_builder, run_estimator, run_codon_bias, run_input_check, run_output_check, run_tool_representation'
make_assembly_help = 'Procedure to produce a genetic tool sequence in genbank format. This module also estimates primers for Gibson assembly.'
make_preprocess_help = 'Procedure to remove the backbone sequence from the genetic tool sequence. The output is in genbank format.'
mapping_targets_help = 'Procedure to map target motifs on the sequence genetic tool.'
partial_generalization_help = 'Procedure to generalize previously indicated target positions of the sequence of the genetic tool.'
candidate_builder_help = 'Procedure to produce a candidate syngeneic sequence of the genetic tool.'
run_estimator_help = 'Procedure to run the complete pipeline of estimation of candidate syngeneic sequence for the input genetic tool.'
run_codon_bias_help = 'Procedure to estimate codon usage for the input genome sequence.'
run_input_check_help = 'Procedure to check format and annotations of the input genetic sequence.'
run_output_check_help = 'Procedure to check format and annotations of the cadidate synegenic sequence for the genetic tool.'
run_tool_representation_help = 'Procedure to produce a tabular representation of the input sequence (in genbank format).'
input_sequence_help = 'Input genbank tool sequence.'
output_sequence_help = 'Output genbank tool sequence.'
input_backbones_help = 'Input fasta backbone sequence.'
input_rm_systems_help = 'Input file containing RM systems information.'
input_strain_genome_help = 'Input strain genome in genbank format.'
codon_table_help = 'Code of the codon table.'
target_positions_help = 'Target positions to consider on the input tool sequence.'
input_config_file_help = 'Input configuration file.'
verbose_help = 'Make SyToGen verbose.'
output_folder_help = 'Output folder path where the outputs from SyToGen will be saved.'
motifs_finder_help = 'A tool that takes as input multiple rebase details and multiple genetic tools (in .gbk format) and generate a static circular representation of the genetic tool and flags where the RM target motifs are found, and a table of the RM target motifs and the number of times each exists in a coding sequence or in a non-coding region.'
input_tools_help = 'List of input_tools as space separated paths (sytogen --motifs_finder file1.gbk file2.gbk). ALso a wildcard indicating all the files with given pattern name can be used (sytogen --motifs_finder  ./path/*.gbk)'
input_rebase_files_help = 'List of input_rebase_files as space separated paths (sytogen --motifs_finder file1.txt file2.txt). ALso a wildcard indicating all the files with given pattern name can be used (sytogen --motifs_finder  ./path/*.txt)'
project_name_help = 'String indicating the name of the project (default: SYT)!'
image_width_help = 'Int indicating the width of the image in px (default: 9).'

available_codon_tables = sorted(list(Bio.Data.CodonTable.generic_by_id.keys()))
help_codon_tables = ''

for key in available_codon_tables:
    if key != 1:
        help_codon_tables += '\n'
    help_codon_tables += '\t' + str(key) + '. ' + ', '.join(
        [name for name in Bio.Data.CodonTable.generic_by_id[key].names if name != None])

def check_required_params(parser):

    if '--make_assembly' in sys.argv:
        required_params = ['input_sequence', 'output_folder']
    elif '--make_preprocess' in sys.argv:
        required_params = ['input_sequence', 'input_backbone', 'output_folder']
    elif '--mapping_targets' in sys.argv:
        required_params = ['input_sequence', 'input_rm_systems', 'input_strain_genome', 'codon_table', 'output_folder']
    elif '--partial_generalization' in sys.argv:
        required_params = ['input_sequence', 'input_rm_systems', 'input_strain_genome', 'codon_table', 'synpl_folder', 'output_folder', 'NCORES']
    elif '--candidate_builder' in sys.argv:
        required_params = ['input_sequence', 'input_rm_systems', 'auto', 'verbose', 'synpl_folder', 'output_folder']
    elif '--run_estimator' in sys.argv:
        required_params = ['input_sequence', 'input_rm_systems', 'input_strain_genome', 'codon_table', 'target_positions', 'output_folder', 'NCORES']
    elif '--run_codon_bias' in sys.argv:
        required_params = ['input_strain_genome', 'codon_table', 'output_folder']
    elif '--run_input_check' in sys.argv:
        required_params = ['input_sequence', 'output_folder']
    elif '--run_output_check' in sys.argv:
        required_params = ['input_sequence', 'input_sequence_new', 'input_rm_systems', 'codon_table', 'output_folder']
    elif '--run_tool_representation' in sys.argv:
        required_params = ['input_sequence', 'input_rm_systems', 'output_folder']
    elif '--motifs_finder' in sys.argv:
        required_params = ['input_tools', 'input_rm_files', 'project_name', 'image_width', 'output_folder']
    
    if 'input_sequence' in required_params:
        parser.add_argument('-i', '--input_sequence', type=str, required=True, help=input_sequence_help)
    
    if 'input_sequence_new' in required_params:
        parser.add_argument('-s', '--input_sequence_new', type=str, required=True, help=output_sequence_help)
    
    if 'input_backbone' in required_params:
        parser.add_argument('-b', '--input_backbone', type=str, required=True, help=input_backbones_help)
    
    if 'input_rm_systems' in required_params:
        parser.add_argument('-r', '--input_rm_systems', type=str, required=True, help=input_rm_systems_help)
    
    if 'auto' in required_params:
        parser.add_argument('-a', '--auto', type=int, required=False, help='', choices=[0, 1, 2], default=1)

        args, remaining = parser.parse_known_args()

        if args.auto == 0:
          parser.add_argument('-f', '--input_config_file', type=str, required=True, help=input_config_file_help)

    if 'synpl_folder' in required_params:
        parser.add_argument('-s', '--synpl_folder', type=str, required=True, help='')

    if 'input_strain_genome' in required_params:
        parser.add_argument('-g', '--input_strain_genome', type=str, required=True, help=input_strain_genome_help)
    
    if 'codon_table' in required_params:
        parser.add_argument('-c', '--codon_table', type=int, required=False, default=11, choices=available_codon_tables, help=codon_table_help) 
    
    if 'target_positions' in required_params:
        parser.add_argument('-t', '--target_positions', type=str, required=False, default=None, help=target_positions_help) 

    if 'output_folder' in required_params:
        parser.add_argument('-o', '--output_folder', type=str, required=True, help=output_folder_help)
        
    if 'NCORES' in required_params:
        parser.add_argument('-nc', '--NCORES', type=int, required=False, help='The number of cores to use', default=1)

    if 'input_tools' in required_params:
        parser.add_argument('-i', '--input_tools', type=str, nargs = "+", required=True, help=input_tools_help)
    if 'input_rm_files' in required_params:
        parser.add_argument('-r', '--input_rm_files', type=str, nargs = "+", required=True, help=input_rebase_files_help)
    if 'project_name' in required_params:
        parser.add_argument('-p', '--project_name', type=str, default='SYT', help=project_name_help)
    if 'image_width' in required_params:
        parser.add_argument('-w', '--image_width', type=int, default=9, help=image_width_help)

    parser.add_argument('--verbose', action='store_true', required=False, help=verbose_help)

    return required_params

def check_inputs(args, required_params):

    if 'input_sequence' in required_params:
        if not os.path.exists(args['input_sequence']):
            error(f"Input sequence file {args['input_sequence']} does not exists!")

        if not os.path.isfile(args['input_sequence']):
            error(f"Input sequence file {args['input_sequence']} is not a file!")            
    
    if 'input_sequence_new' in required_params:
        if not os.path.exists(args['input_sequence_new']):
            error(f"Input alternative sequence file {args['input_sequence_new']} does not exists!")

        if not os.path.isfile(args['input_sequence_new']):
            error(f"Input alternative sequence file {args['input_sequence_new']} is not a file!")   
    
    if 'input_backbone' in required_params:
        if not os.path.exists(args['input_backbone']):
            error(f"Input backbone sequence file {args['input_backbone']} does not exists!")

        if not os.path.isfile(args['input_backbone']):
            error(f"Input backbone sequence file {args['input_backbone']} is not a file!")
    
    if 'input_rm_systems' in required_params:
        if not os.path.exists(args['input_rm_systems']):
            error(f"Input RM file {args['input_rm_systems']} does not exists!")

        if not os.path.isfile(args['input_rm_systems']):
            error(f"Input RM file {args['input_rm_systems']} is not a file!")
        
        out_rm_information = check_rm_information(args['input_rm_systems'])

        if not out_rm_information['correct']:
            error('\n'.join(out_rm_information['messages']))
    
    if 'input_strain_genome' in required_params:
        if not os.path.exists(args['input_strain_genome']):
            error(f"Input genome file {args['input_strain_genome']} does not exists!")

        if not os.path.isfile(args['input_strain_genome']):
            error(f"Input genome file {args['input_strain_genome']} is not a file!")

    if 'input_config_file' in required_params:
        if not os.path.exists(args['input_strain_genome']):
            error(f"Input genome file {args['input_strain_genome']} does not exists!")

        if not os.path.isfile(args['input_strain_genome']):
            error(f"Input genome file {args['input_strain_genome']} is not a file!")
    
    if 'output_folder' in required_params:
        (os.path.exists(args['output_folder']) & os.path.isdir(args['output_folder']))

        if not os.path.exists(args['output_folder']):
            error(f"Output folder {args['output_folder']} does not exists!")

        if not os.path.isdir(args['output_folder']):
            error(f"Output folder {args['output_folder']} is not a file!")

    # Motifs finder module

    if 'input_tools' in required_params:

        input_tools = args['input_tools']

        for input_tool in input_tools:
            if not os.path.exists(input_tool):
                error(f"Input sequence file {input_tool} does not exists!")

            if not os.path.isfile(input_tool):
                error(f"Input sequence file {input_tool} is not a file!")
    
    if 'input_rm_files' in required_params:

        input_rm_files = args['input_rm_files']

        for input_rm_file in input_rm_files:
            if not os.path.exists(input_rm_file):
                error(f"Input sequence file {input_rm_file} does not exists!")

            if not os.path.isfile(input_rm_file):
                error(f"Input sequence file {input_rm_file} is not a file!")

    return

def sytogen_main():

    parser = ap.ArgumentParser(description=("SyToGen is an automatic pipeline for estimating candidate syngeneic DNA and mobile "
                                        "genetic elements that can be accepted as self by the target microbial strain. Refer "
                                        "to the online documentation {} for a detailed description of SyToGen functionalities. ".format("..")
                                        ), formatter_class=ap.RawTextHelpFormatter)
    
    group = parser.add_mutually_exclusive_group()

    group.add_argument('--make_assembly', action='store_true', required=False, help=make_assembly_help)
    group.add_argument('--make_preprocess', action='store_true', required=False, help=make_preprocess_help)
    
    group.add_argument('--mapping_targets', action='store_true', required=False, help=mapping_targets_help)
    group.add_argument('--partial_generalization', action='store_true', required=False, help=partial_generalization_help)
    group.add_argument('--candidate_builder', action='store_true', required=False, help=candidate_builder_help)
    
    group.add_argument('--run_estimator', action='store_true', required=False, help=run_estimator_help)

    group.add_argument('--run_codon_bias', action='store_true', required=False, help=run_codon_bias_help)
    group.add_argument('--run_input_check', action='store_true', required=False, help=run_input_check_help)
    group.add_argument('--run_output_check', action='store_true', required=False, help=run_output_check_help)
    group.add_argument('--run_tool_representation', action='store_true', required=False, help=run_tool_representation_help)

    group.add_argument('--motifs_finder', action='store_true', required=False, help=motifs_finder_help)

    params = [el for el in sys.argv[1:] if '-' in el]

    commands = [
        '--make_assembly',
        '--make_preprocess',
        '--mapping_targets',
        '--partial_generalization',
        '--candidate_builder',
        '--run_estimator',
        '--run_codon_bias',
        '--run_input_check',
        '--run_output_check',
        '--run_tool_representation',
        '--motifs_finder'
    ]

    if not len(params) >= 1:
        info('Incorrect input parameters!', init_new_line=False)
        error(f'Run python{sys.version_info[0]}.{sys.version_info[1]} sytogen.py -h to get basic usage!')

    if params[0] in ['-h', '--help']:
        info(f"basic usage is:\npython{sys.version_info[0]}.{sys.version_info[1]} sytogen.py [{'|'.join(commands)}]\n", init_new_line=False, exit=True)

    if not params[0] in commands:
        info(f'Incorrect input parameters: {params[0]} not in {commands}')
        error(f'Run python{sys.version_info[0]}.{sys.version_info[1]} sytogen.py -h to get basic usage!')

    required_params = check_required_params(parser)

    args = parser.parse_args()

    check_inputs(vars(args), required_params)

    if params[0] == '--make_assembly':

        make_assembly(args)

        info('--make_assembly run successfully!\n')

    elif params[0] == '--make_preprocess':

        make_preprocess(args)
        
        info('--make_preprocess run successfully!\n')

    elif params[0] == '--mapping_targets':
        
        run_estimator(args, type_action = 'mapping_targets')
        
        info('--mapping_targets run successfully!\n')

    elif params[0] == '--partial_generalization':

        run_estimator(args, type_action = 'partial_generalization')
        
        info('--partial_generalization run successfully!\n')

    elif params[0] == '--candidate_builder':

        candidate_builder(args)
        
        info('--candidate_builder run successfully!\n')

    elif params[0] == '--run_estimator':
        
        run_estimator(args)

        info('--run_estimator run successfully!\n')

    elif params[0] == '--run_codon_bias':

        run_codon_bias(args)
        
        info('--run_codon_bias run successfully!\n')

    elif params[0] == '--run_input_check':

        run_input_check(args)
        
        info('--run_input_check run successfully!\n')

    elif params[0] == '--run_output_check':

        run_output_check(args)
        
        info('--run_output_check run successfully!\n')

    elif params[0] == '--run_tool_representation':
        
        run_tool_representation(args)
        
        info('--run_tool_representation run successfully!\n')

    elif params[0]  == '--motifs_finder':

        # run_motifs_finder(args)
        run_motifs_finder2(args)

        info('--motifs_finder run successfully!\n')

    return

if __name__ == '__main__':
    t0 = time.time()
    sytogen_main()
    t1 = time.time()
    info(f'\nTotal elapsed time {int(t1 - t0)}s\n')
    sys.exit(0)
#!/usr/bin/env python

import re
import argparse as ap
import pandas as pd
import os

import Bio
from Bio import Seq, SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable

from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import CompoundLocation

import sys

import warnings

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
        out_config_file['messages'].append('Number of input sequence > 1.')
        return out_config_file

    with open(input_file, 'r') as input_handle:
        tmp_lines = ''.join(input_handle.readlines())

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
                if not ''.join(sorted(list(set(tmp_str.upper())))) == 'ACGT':
                    out_config_file['correct'] = False

                out_config_file['total_cds_concats'] += tmp_str
        j += 1

    return out_config_file

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

if __name__ == '__main__':

    parser = ap.ArgumentParser(description=
    'This script takes as input gff, gbk and single or multifasta files and \n'
    'compute the codon usage for a specified codon table.\n'
    'Usage:\n'
    'python codon_bias_estimator.py -i example.gbk -t genebank -o gbk_example -c 11\n'
    'python codon_bias_estimator.py -i example.ffn -t fasta -o fasta_example -c 11\n'
    'python codon_bias_estimator.py -i example.gff -t gff -o gff_example -c 11\n',
    formatter_class=ap.RawTextHelpFormatter)

    parser.add_argument('-i', '--input_genome',
                        help='The path to the input file', required=True)

    # parser.add_argument('-t','--type', help=
    # 'The format of the file [genebank, fasta, ]', required=True)

    parser.add_argument('-c','--codonTable', help=
        'The codon table code to be used [1, 2, 3 ...]\n'
        '1. The Standard Code,\n'
        '2. The Vertebrate Mitochondrial Code,\n'
        '3. The Yeast Mitochondrial Code,\n'
        '4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code,\n'
        '5. The Invertebrate Mitochondrial Code,\n'
        '6. The Ciliate, Dasycladacean and Hexamita Nuclear Code,\n'
        '9. The Echinoderm and Flatworm Mitochondrial Code,\n'
        '10. The Euplotid Nuclear Code,\n'
        '11. The Bacterial, Archaeal and Plant Plastid Code,\n'
        '12. The Alternative Yeast Nuclear Code,\n'
        '13. The Ascidian Mitochondrial Code,\n'
        '14. The Alternative Flatworm Mitochondrial Code,\n'
        '16. Chlorophycean Mitochondrial Code,\n'
        '21. Trematode Mitochondrial Code,\n'
        '22. Scenedesmus obliquus Mitochondrial Code,\n'
        '23. Thraustochytrium Mitochondrial Code,\n'
        '24. Rhabdopleuridae Mitochondrial Code,\n'
        '25. Candidate Division SR1 and Gracilibacteria Code,\n'
        '26. Pachysolen tannophilus Nuclear Code,\n'
        '27. Karyorelict Nuclear Code,\n'
        '28. Condylostoma Nuclear Code,\n'
        '29. Mesodinium Nuclear Code,\n'
        '30. Peritrich Nuclear Code,\n'
        '31. Blastocrithidia Nuclear Code,\n'
        '33. Cephalodiscidae Mitochondrial UAA-Tyr Code,\n', required=True)

    parser.add_argument(
        '-o', '--output_folder', help='output_dir', required=True)
                
    args = vars(parser.parse_args())

    codonTable = args['codonTable']
    input_genome = args['input_genome']

    input_genome = read_input_seq(input_genome)

    if input_genome['correct']:        
        genome_codon_usage = codon_usage(input_genome['total_cds_concats'], codonTable)
        
        with open(os.path.join(args['output_folder'], "codon_usage_table.csv"), 'w') as outf:
            genome_codon_usage['Table'].to_csv(outf, sep=',', index=False)

    else:
        info('\n'.join(input_genome['messages']))
        error('The input genome is misannotated')
#!/usr/bin/env python
'''
This tool extracts type sequences from the MLST-like schemes created by mlst_from_srst2.py.

It outputs a concatanated sequence of the constituent genes of one or more types. If multiple types
are requested, this tool can align them as well (requires Muscle to be installed).

Author: Ryan Wick
email: rrwick@gmail.com
'''

from __future__ import print_function
from __future__ import division
import argparse
from itertools import groupby
from mlst_from_srst2 import *



def main():
    args = get_arguments()
    gene_seqs = load_fasta_as_dict(args.gene_seqs)
    mlst = MlstScheme(scheme_table=args.scheme)
    seqs_by_type_and_gene = {} # key = ST value = list of gene sequences
    for st_num in get_types(args, mlst):
        seqs_by_type_and_gene[st_num] = []
        alleles = mlst.type_to_alleles[st_num]
        for allele in alleles:
            if allele not in gene_seqs:
                quit_with_error('Allele ' + allele + ' not in gene sequence FASTA')
            seqs_by_type_and_gene[st_num].append(gene_seqs[allele])
    if args.align:
        seqs_by_type_and_gene = align_seqs(seqs_by_type_and_gene)
    cat_seqs = sorted([(st, ''.join(seqs)) for st, seqs in seqs_by_type_and_gene.iteritems()])
    save_fasta(cat_seqs, args.out)


def get_arguments():
    '''
    Specify the command line arguments required by the script.
    '''
    parser = argparse.ArgumentParser(description='MLST from SRST2 - get type sequences')
    parser.add_argument('-s', '--scheme', type=str, required=True,
                        help='MLST scheme file')
    parser.add_argument('-g', '--gene_seqs', type=str, required=True,
                        help='FASTA file of gene sequences')
    parser.add_argument('-t', '--types', type=str, required=False,
                        help='Comma-delimited list of types to output')
    parser.add_argument('--all', action='store_true',
                        help='Output all types in the scheme')
    parser.add_argument('-o', '--out', type=str, required=True,
                        help='Filename for type sequence output')
    parser.add_argument('-a', '--align', action='store_true',
                        help='Align type sequences (only applicable when more than one type '
                             'sequences are outputted)')
    args = parser.parse_args()
    if not args.types and not args.all:
        parser.error('You must either provide one or more types (--types) or all types (--all) '
                     'for sequence output.')
    if args.types and args.all:
        parser.error('The --types and --all arguments cannot be both used together.')
    return parser.parse_args()

def get_types(args, mlst):
    '''
    This function returns a list of the sequence types to output based on the user's arguments.
    It also checks for problems with the chosen types.
    '''
    if args.all:
        types = mlst.type_to_alleles.keys()
        if not types:
            quit_with_error('The given MLST scheme has no sequence types.')
    else:
        types = args.types.split(',')
    types = sorted([string_to_int(x) for x in types])
    if None in types:
        quit_with_error('One or more sequence types is incorrectly formatted.')
    for st_num in types:
        if st_num not in mlst.type_to_alleles:
            quit_with_error('Sequence type ' + str(st_num) + ' is not in the given MLST scheme.')
    if len(types) != len(set(types)):
        quit_with_error('Not all sequence types are unique (there are duplicates).')
    return types

def align_seqs(seqs_by_type_and_gene):
    '''
    Uses Muscle to align the sequences. Alignments are performed independently for each gene (the
    sequences will be concatenated later).
    '''
    # TO DO
    # TO DO
    # TO DO
    # TO DO
    # TO DO
    # TO DO
    return aligned_seqs

def load_fasta_as_dict(fasta_filename):
    '''
    Returns the contents of the FASTA file as a dictionary where Key = header and Value = sequence.
    '''
    return dict(fasta_iter(fasta_filename))

def string_to_int(int_str):
    '''
    This function converts a string to an integer or to None, if the string can't be converted to
    an integer.
    '''
    try:
        integer = int(int_str)
        return integer
    except ValueError:
        return None

def fasta_iter(fasta_filename):
    '''
    https://www.biostars.org/p/710/
    '''
    fasta = open(fasta_filename)
    faiter = (x[1] for x in groupby(fasta, lambda line: line[0] == '>'))
    for header in faiter:
        header = header.next()[1:].strip()
        seq = ''.join(s.strip() for s in faiter.next())
        yield header, seq

def save_fasta(st_seqs, filename):
    '''
    Saves a FASTA file using the list of tuples where the first part is the header and the second
    is the sequence.
    '''
    fasta_file = open(filename, 'w')
    for st_num, seq in st_seqs:
        fasta_file.write('>ST' + str(st_num) + '\n')
        fasta_file.write(add_line_breaks_to_sequence(seq, 60))

def add_line_breaks_to_sequence(sequence, length):
    '''
    Wraps sequences to the defined length. All resulting sequences end in a line break.
    '''
    seq_with_breaks = ''
    while len(sequence) > length:
        seq_with_breaks += sequence[:length] + '\n'
        sequence = sequence[length:]
    if sequence:
        seq_with_breaks += sequence
        seq_with_breaks += '\n'
    return seq_with_breaks

if __name__ == '__main__':
    main()

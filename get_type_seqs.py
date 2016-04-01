#!/usr/bin/env python
'''
This tool extracts type sequences from the MLST-like schemes created by mlst_from_srst2.py.

It outputs a concatanated sequence of the constituent genes of one or more sequence types. If
multiple types are requested, this tool can align them as well (requires Muscle to be installed).

Author: Ryan Wick
email: rrwick@gmail.com
'''

from __future__ import print_function
from __future__ import division
import argparse
import subprocess
import sys
from mlst_from_srst2 import MlstScheme, quit_with_error

def main():
    '''
    Script execution starts here.
    '''
    args = get_arguments()
    if args.align and not find_program('muscle'):
        quit_with_error('Muscle must be installed to produce an aligned output.')
    gene_seqs = dict(load_fasta_file(args.gene_seqs))
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
        seqs_by_type_and_gene = align_seqs(seqs_by_type_and_gene, args.muscle_args)
    cat_seqs = sorted([(st, ''.join(seqs)) for st, seqs in seqs_by_type_and_gene.iteritems()])
    save_fasta(cat_seqs, args.out)

def get_arguments():
    '''
    Specify the command line arguments required by the script.
    '''
    fix_muscle_args()
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
    parser.add_argument('--muscle_args', type=str, required=False, default='',
                        help='Additional parameters (enclosed in quotes) to be passed to the '
                             'Muscle aligner')
    args = parser.parse_args()
    if not args.types and not args.all:
        parser.error('You must either provide one or more types (--types) or all types (--all) '
                     'for sequence output.')
    if args.types and args.all:
        parser.error('The --types and --all arguments cannot be both used together.')
    if args.muscle_args and not args.align:
        parser.error('--muscle_args requires the use of --align.')
    return parser.parse_args()

def fix_muscle_args():
    '''
    This function looks to see if the --muscle_args argument was used, and if so, it will possibly
    add a space before the following piece to prevent issues with argparse.
    '''
    if '--muscle_args' in sys.argv:
        i = sys.argv.index('--muscle_args') + 1
        if i < len(sys.argv) and sys.argv[i].startswith('-'):
            sys.argv[i] = ' ' + sys.argv[i]

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

def align_seqs(seqs_by_type, muscle_args):
    '''
    Uses Muscle to align the sequences. Alignments are performed independently for each gene (the
    sequences will be concatenated later).
    '''
    st_nums = seqs_by_type.keys()
    aligned_seqs_by_type = {st_num: [] for st_num in st_nums}
    gene_count = len(seqs_by_type.itervalues().next())
    command = ['muscle'] + muscle_args.split()
    for i in range(gene_count):
        gene_seqs = [(st_num, seqs[i]) for st_num, seqs in seqs_by_type.iteritems()]
        muscle_input = ''
        for gene_seq in gene_seqs:
            muscle_input += '>' + str(gene_seq[0]) + '\n'
            muscle_input += add_line_breaks_to_sequence(gene_seq[1], 60)
        process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        muscle_output, err = process.communicate(input=muscle_input)
        if '*** ERROR ***' in err:
            muscle_error = err.split('*** ERROR ***')[1].strip()
            quit_with_error('Muscle alignment failed\n' + muscle_error)
        if 'Invalid command line option' in err:
            muscle_error = err.split('\n')[0].strip()
            quit_with_error('Muscle alignment failed\n' + muscle_error)
        aligned_seqs = load_fasta_lines(muscle_output)
        for st_num, seq in aligned_seqs:
            aligned_seqs_by_type[int(st_num)].append(seq)
    return aligned_seqs_by_type

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

def load_fasta_file(filename): # type: (str) -> list[tuple[str, str]]
    '''
    Returns the names and sequences for the given fasta file.
    '''
    fasta_file = open(filename, 'r')
    return load_fasta_lines(fasta_file.read())

def load_fasta_lines(fasta_str): # type: (str) -> list[tuple[str, str]]
    '''
    Takes as input a single string for the whole FASTA file.
    '''
    fasta_seqs = []
    name = ''
    sequence = ''
    for line in line_iterator(fasta_str):
        line = line.strip()
        if not line:
            continue
        if line[0] == '>': # Header line = start of new contig
            if name:
                fasta_seqs.append((name.split()[0], sequence))
                name = ''
                sequence = ''
            name = line[1:]
        else:
            sequence += line
    if name:
        fasta_seqs.append((name.split()[0], sequence))
    return fasta_seqs

def line_iterator(string_with_line_breaks):
    '''
    Iterates over a string containing line breaks, one line at a time.
    '''
    prev_newline = -1
    while True:
        next_newline = string_with_line_breaks.find('\n', prev_newline + 1)
        if next_newline < 0:
            break
        yield string_with_line_breaks[prev_newline + 1:next_newline]
        prev_newline = next_newline

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

def find_program(name): # type: (str) -> bool
    '''
    Checks to see if a program exists.
    '''
    process = subprocess.Popen(['which', name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    return bool(out) and not bool(err)

if __name__ == '__main__':
    main()

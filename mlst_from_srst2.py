#!/usr/bin/env python
'''
MLST from SRST2

This is a tool to create or expand an MLST-like scheme using a table of
compiled gene allele results from SRST2.

Author: Ryan Wick
email: rrwick@gmail.com
'''

from __future__ import print_function
from __future__ import division
import argparse
import os

def main():
    args = get_arguments()

    # If an existing scheme was given, load it in.  If not, create an empty one.
    if (args.existing_scheme):
        mlst = MlstScheme(scheme_table = args.existing_scheme)
    else:
        mlst = MlstScheme(gene_list = args.mlst_genes)

    sample_types = open(args.out_types, 'w')
    sample_types.write('Sample\tST\talleles\n')
    existing_type_assignments = 0
    new_type_assignments = 0
    failed_type_assignments = 0

    with open(args.srst2_table , 'r') as table:

        # Read the first line of the SRST2 table and build a dictionary of
        # column number -> gene (cluster) name
        first_line = table.readline()
        table = open(args.srst2_table, 'r')
        first_line = table.readline()
        line_parts = first_line.strip().split('\t')
        gene_columns = {}
        for i in range(1, len(line_parts)):
            gene_columns[i] = line_parts[i]

        # Read each sample line in the SRST2 table and get its sequence type.
        for line in table:
            line_parts = line.strip().split('\t')
            assert len(line_parts) <= len(gene_columns) + 1
            sample = line_parts[0]
            sample_alleles = {}
            for i in range(1, len(line_parts)):
                gene = gene_columns[i]
                allele = line_parts[i]
                sample_alleles[gene] = allele
            sequence_type, status, allele_list = mlst.get_sequence_type(sample_alleles)
            sample_types.write(sample + '\t' + str(sequence_type) + '\t' + allele_list + '\n')
            if status == 'existing':
                existing_type_assignments += 1
            elif status == 'new':
                new_type_assignments += 1
            elif status == 'fail':
                failed_type_assignments += 1

    # Write the new MLST scheme to file.
    new_scheme = open(args.out_scheme, 'w')
    new_scheme.write(str(mlst))

    # Output results to user.
    successes = existing_type_assignments + new_type_assignments
    if successes:
        plural = ('' if successes == 1 else 's')
        print('Successfully assigned sequence type' + plural + ' to ' + str(successes) + ' sample' + plural)
    if failed_type_assignments:
        plural = ('' if failed_type_assignments == 1 else 's')
        print('Failed to assign sequence type' + plural + ' to ' + str(failed_type_assignments) + ' sample' + plural)
    if new_type_assignments:
        plural = ('' if new_type_assignments == 1 else 's')
        print('Created ' + str(new_type_assignments) + ' new sequence type' + plural)
    else:
        print('All assigned sequence types were previously known (MLST scheme is unchanged)')



def get_arguments():
    parser = argparse.ArgumentParser(description='MLST from SRST2')

    parser.add_argument('-s', '--srst2_table', type=str, required=True, help='Table of SRST2 results')
    parser.add_argument('-e', '--existing_scheme', type=str, required=False, help='Existing MLST scheme')
    parser.add_argument('-g', '--mlst_genes', type=str, required=False, help='Comma-delimited list of genes to use in new scheme')
    parser.add_argument('-t', '--out_types', type=str, required=True, help='Output file of sample type assignments')
    parser.add_argument('-o', '--out_scheme', type=str, required=True, help='Output file of updated MLST scheme')

    args = parser.parse_args()

    if not args.existing_scheme and not args.mlst_genes:
        parser.error("When an existing MLST scheme is not used, a list of MLST genes (--mlst_genes) is required.")
    if args.existing_scheme and args.mlst_genes:
        parser.error("When an existing MLST scheme is used, a list of MLST genes (--mlst_genes) cannot be used, as the script will use the genes in the existing scheme.")
    
    return parser.parse_args()




class MlstScheme:
    '''
    Class members:
    genes = a list of the gene clusters used in the scheme
    sequence_types = a dictionary where key = allele list (as a comma-delimited string) and value = ST number
    '''

    def __init__(self, scheme_table = None, gene_list = None):
        '''
        This constructor can be run in two ways: either with both
        scheme_table and seq_table, or with just gene_list.
        The first way will build the object using the existing scheme.
        The second way will create an empty scheme with the given genes.
        '''
        self.sequence_types = {}

        # Load in the existing scheme.
        if scheme_table:
            scheme_table_file = open(scheme_table, 'r')
            for line in scheme_table_file:
                stripped_line = line.strip()
                if not stripped_line:
                    continue
                line_parts = stripped_line.split('\t')

                # The first line is a header line.
                if line_parts[0] == 'ST':
                    self.genes = line_parts[1:]

                # Other lines are sequence types.
                else:
                    assert len(line_parts) == len(self.genes) + 1
                    st_num = int(line_parts[0])
                    allele_list = ','.join(line_parts[1:])
                    self.sequence_types[allele_list] = st_num

        # Create an empty scheme.
        elif gene_list:
            self.genes = gene_list.split(',')
            if len(self.genes) < 2:
                print('Error: at least two genes needed to make an MLST scheme')
                quit()

    def get_sequence_type(self, genes_and_alleles):
        '''
        This function takes a dictionary of genes->alleles and returns a
        sequence type for this MLST scheme.
        It also returns the status of the assignment:
           'existing' when the type was already in the MLST scheme
           'new' when the MLST scheme was expanded to make this type
           'fail' when it could not assign a type.
        It also returns the allele list in string form.
        '''
        allele_list = []
        failed = False
        for gene in self.genes:
            if gene not in genes_and_alleles:
                allele_list.append('?')
                failed = True
            else:
                allele_list.append(genes_and_alleles[gene])
        allele_list_string = ','.join(allele_list)

        assert len(allele_list) == len(self.genes)
        for allele in allele_list:
            if '*' in allele or '?' in allele or allele == '-':
                failed = True
        if failed:
            return ('?', 'fail', allele_list_string)

        if allele_list_string in self.sequence_types:
            return (self.sequence_types[allele_list_string], 'existing', allele_list_string)
        else:
            return (self.make_new_sequence_type(allele_list_string), 'new', allele_list_string)

    def make_new_sequence_type(self, allele_list_string):
        '''
        This function creates a new sequence type for the allele list given
        (in a comma-delimited string) and returns the number for the new type.
        It assumes that the allele list given does not match any existing
        sequence types.
        '''
        largest_current_sequence_type_num = 0
        if self.sequence_types:
            largest_current_sequence_type_num = max(self.sequence_types.values())
        new_num = largest_current_sequence_type_num + 1
        self.sequence_types[allele_list_string] = new_num
        return new_num

    def __str__(self):
        '''
        Produces a tab-delimited table for the MLST scheme.
        '''
        scheme_string = 'ST\t' + '\t'.join(self.genes) + '\n' # header line
        sequence_types_list = self.sequence_types.items()
        sequence_types_list = sorted(sequence_types_list, key=lambda x: x[1])
        for sequence_type in sequence_types_list:
            allele_list = sequence_type[0].split(',')
            num = sequence_type[1]
            scheme_string += str(num) + '\t' + '\t'.join(allele_list) + '\n'
        return scheme_string






if __name__ == '__main__':
    main()

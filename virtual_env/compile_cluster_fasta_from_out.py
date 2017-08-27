#!/usr/bin/env python
# This script will work on a .out directory, and read the .out files
# to produce an output fasta file of all the types of WP ids matching a given
# query name, or something (ie. pngc)

# The intent of this is to produce an output fasta file, a list of these related sequences,
# for further analysis by the perl script mkBlastClusters.pl

from __future__ import print_function

import argparse
import os
import sys
import pandas as pd

from collections import defaultdict

from util import parse_faa, write_fasta_sequences
# use the below line to add/change constants to be imported
from util import PPM_MATCH_LIST, PNGC_MATCH_LIST

# for reading outfiles
COLS = ['target_name', 't_accession', 'tlen', 'query_name' ]


def get_outfile_data(outfile_path, valid_target_names):
    # Extract just GCF_xyz from the given file path
    gcf_id = outfile_path.split('/')[-1].split('.')[0]

    with open(out_path, 'r') as out_file:
        # use first 4 columns, fourth column is query name, first column is target name

            input_df = pd.read_csv(os.path.join(INPUT_PATH, fname),
                comment='#',
                header=None,
                delimiter='\s+',
                usecols=range(4)
                )
            input_df.columns = COLS  # add the column names to the dataframe
            input_df = input_df[['query_name', 'target_name']]

            # Get the unique query names from input_df where target_name is in valid_target_names
            valid_rows = input_df[input_df['target_name'] in valid_target_names]
            valid_query_names = valid_rows['query_name'].unique()
            print('Hello! {} '.format(valid_query_names))

    # Now I want to build a dictionary of GCF_ to list of WP_ ids.
    wp_dict = defaultdict(set)


    return wp_dict
 

def write_fasta_output(outfile_path, fasta_path):
    wp_data = {}

    # So, first, we iterate over each .out file in the given path. (each!)
    out_files = [filename for filename in os.listdir(outfile_path) if filename.endswith('.out')]
    for filename in out_files:
        file_path = '/'.join([outfile_path, filename])

        # Now we need to READ the out file for the data we want.
        wp_data.update(get_outfile_data(file_path))

    for gcf_id, wp_list in distance_data.items():
        write_fasta_sequences(gcf_id, wp_list, fasta_path)


if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
            description=('Creates a new fasta file containing all sequences '
            'referenced in a given protein_distances.py output file.'))
    argparser.add_argument('outfile_path',
            help='Path to the input file, created by protein_distances.py')
    argparser.add_argument('fasta_path',
            help='Path to directory where relevant fasta files are located.')
    args = argparser.parse_args()

    write_fasta_output(args.outfile_path, args.fasta_path)

#!/usr/bin/env python
# This script will take a list of GCF/WPs, intended to be a list of similar WPs
# - similar, as in, they represent the SAME kind of protein (PPM, PNGC, LICC, etc)

# The intent of this is to produce an output fasta file, a list of these related sequences,
# for further analysis by the perl script mkBlastClusters.pl

from __future__ import print_function

import argparse
import os
import sys

from collections import defaultdict

from util import parse_faa


def get_distance_data(distance_path):
    with open(distance_path, 'r') as distance_file:
        distances_raw = distance_file.readlines()
    # Skip the first line, just column headers - 
    distances_raw = distances_raw[1:]

    # Now I want to build a dictionary of GCF_ to list of WP_ ids.
    wp_dict = defaultdict(set)
    gcf_id = None

    for line in distances_raw:
        if line.startswith('#'):
            # It's a comment line - and it either represents a GCF, or an error
            try:
                line_data = line.split()[1:]
            except Exception as e:
                print("Failed to parse comment line, empty comment. Error: {}".format(e),
                        file=sys.stderr)
                continue

            if line_data[0].startswith('GCF_'):
                # We're good to go - 
                gcf_id = line_data[0].split('.')[0] # remove the .out from the end

        else:
            if gcf_id is None:
                raise Exception("Error: Can't start adding data before finding a valid GCF id")

            # The line doesn't start with a comment symbol - it's a data line!
            line_data = line.split()
            for item in line_data[:2]: # up to 2 - ignore the distance value
                wp_dict[gcf_id].add(item)

    return wp_dict


def write_fasta_sequence(wp_id, sequence_data, line_length=80):
    print('>{}'.format(wp_id))
    for i in range(0, len(sequence_data), line_length):
        print(sequence_data[i:i+line_length])


def write_sequences(gcf_id, wp_ids, fasta_path):
    # step one - open the gcf fasta file
    try:
        faa_data = parse_faa('/'.join([fasta_path, gcf_id+'.faa']))
    except IOError as e:
        print('# {}'.format(e))
        return

    for wp_id in wp_ids:
        sequence = faa_data[wp_id]
        write_fasta_sequence(wp_id, sequence)


def write_fasta_output(distance_file_path, fasta_path):
    distance_data = get_distance_data(distance_file_path)

    for gcf_id, wp_list in distance_data.items():
        write_sequences(gcf_id, wp_list, fasta_path)


if __name__ == '__main__':
    argparser = argparse.ArgumentParser(
            description=('Creates a new fasta file containing all sequences '
            'referenced in a given protein_distances.py output file.'))
    argparser.add_argument('distances_path',
            help='Path to the input file, created by protein_distances.py')
    argparser.add_argument('fasta_path',
            help='Path to directory where relevant fasta files are located.')
    args = argparser.parse_args()

    write_fasta_output(args.distances_path, args.fasta_path)

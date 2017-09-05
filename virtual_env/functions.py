#!/usr/bin/env python

# 
# 

from __future__ import print_function

import argparse
import gzip
import sys

from protein_distances import parse_lastcol


PATH = '/home/kelly/Dropbox/gff/gff_files'


class GFFProteinData(object):
    def __init__(self, row_data):
        last_col_data = parse_lastcol(row_data[-1]) # call it on the last column

        # It's possible that we get a 'pseudo' protein, which has no Name/protein_id
        self.name = last_col_data.get('Name', last_col_data.get('protein_id', ''))
        # And Parent probably has the form, 'gene12345', so we convert it to 12345
        self.parent = convert_str_int(last_col_data['Parent'])
        self.id = last_col_data['ID']
        self.product = last_col_data['product']
        self.note = last_col_data.get('Note', '')

        # I have no idea what this is, but I don't think it's important.
        # self.data = row_data[0]


def convert_str_int(string):
    filtered = filter(lambda s: s.isdigit(), string)
    return int(filtered)

def collect_protein_data(gcf_id, protein_ids):
    """Returns 2 things: dict of all 'protein data' for given GCF_xyz.gff.gz,
    and also a dict of the indices for the given protein ids to access the data
    from the first dict.
    """
    protein_id_locations = {}
    protein_data = {}
    with gzip.open('/'.join([PATH, gcf_id]) + '.gff.gz', 'r') as f:
        for l in f:
            if l.startswith('#'):
                continue # it's a comment line, ignore

            data = l.split('\t') # we only care about data[0], data[2], data[-1]
            
            if data[2].lower() not in ['protein', 'cds']:
                # possibly: gene, sequence_feature, ...
                continue

            curr_protein = GFFProteinData(data)

            # Do this so we can instantly locate the desired proteins later
            if curr_protein.name in protein_ids:
                protein_id_locations[curr_protein.name] = curr_protein.parent

            protein_data[curr_protein.parent] = curr_protein

    return protein_data, protein_id_locations


def print_annotations(gcf_id, protein_ids, dist):

    protein_data, protein_id_locations = collect_protein_data(gcf_id, protein_ids)
    # Problem with how this works: produces duplicate outfit if desired proteins are within 2*dist of each other
    # cop-out solution: | sort | uniq
    # Proper solution probably involves using a set on i, or like a dictionary to produce output?
    for name, pos in protein_id_locations.items():
        print( '# {} {}'.format(name, pos))
        for i in range(pos-dist, pos):
            protein = protein_data[i]
            print (i, protein.name, protein.product)
        for i in range(pos+1, pos+1+dist):
            protein = protein_data[i]
            print (i, protein.name, protein.product)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=
            """Prints the 'product' value for proteins surrounding given protein ids in a given GCF.
            We understand this to be some form of functional annotation.\n
            output: position | id | product
            """)
    parser.add_argument('gcf', help='GCF identifier, eg. GCF_00024545')
    parser.add_argument('protein_list_file', help='Path to file consisting of a list of protein ids (WP_xyz)')
    parser.add_argument('--dist', default=10, type=int,
            help='The distance about each to produce annotations for, eg. 10')
    args = parser.parse_args()

    with open(args.protein_list_file, 'r') as plf:
        protein_ids = [l.strip() for l in plf.readlines()]

    print_annotations(args.gcf, protein_ids, args.dist)

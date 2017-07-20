#!/usr/bin/env python

# 
# 

from __future__ import print_function

import argparse
import gzip
import sys

from protein_distances import parse_lastcol


PATH = '/home/kelly/Dropbox/gff/gff_files'

def convert_str_int(string):
    filtered = filter(lambda s: s.isdigit(), string)
    return int(filtered)

def print_annotations(gcf_id, protein_ids, dist):
    with gzip.open('/'.join([PATH, gcf_id]) + '.gff.gz', 'r') as f:
        protein_id_locations = {}
        protein_data = {}
        for l in f:
            if l.startswith('#'):
                continue # it's a comment line

            data = l.split('\t') # we only care about data[0], data[2], data[-1]

            # we only care about protein or CDS rows?
            # other possible rows:
            # - gene, sequence_feature, ... ?? # would not be hard to collect
            if data[2].lower() != 'protein' and data[2].lower() != 'cds':
                continue

            lcd = parse_lastcol(data[-1]) # call it on the last column
            lcd['Parent'] = convert_str_int(lcd['Parent'])

            # Do this so we can instantly locate the desired proteins later
            if lcd.get('Name') in protein_ids:
                protein_id_locations[lcd['Name']] = lcd['Parent']
                
            protein_data[ lcd['Parent'] ] = (
                lcd.get( 'protein_id' ),
                lcd['ID'],
                data[0],
                lcd.get('Note'),
                lcd['product']
                )

    # Problem with how this works: produces duplicate outfit if desired proteins are within 2*dist of each other
    # cop-out solution: | sort | uniq
    # Proper solution probably involves using a set on i, or like a dictionary to produce output?
    for name, pos in protein_id_locations.items():
        print( '# {} {}'.format(name, pos))
        for i in range(pos-dist, pos):
            print (i, protein_data[i][1], protein_data[i][4])
        for i in range(pos+1, pos+1+dist):
            print (i, protein_data[i][1], protein_data[i][4])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=
            """Prints the 'product' value for proteins surrounding given protein ids in a given GCF.
            We understand this to be some form of functional annotation.

            The output format is: | position | id | product |

            Program Usage: cat input_list | python functions.py GCF_id 10
            """)
    parser.add_argument('gcf', help='GCF identifier, eg. GCF_00024545')
    parser.add_argument('--dist', default=10, type=int,
            help='The distance about each to produce annotations for, eg. 10')
    args = parser.parse_args()

    # call it with arguments gcf, dist
    # pass in the protein ids from stdin (easier to cat/awk/grep)
    protein_ids = [l.strip() for l in sys.stdin.readlines()]

    print_annotations(args.gcf, protein_ids, args.dist)

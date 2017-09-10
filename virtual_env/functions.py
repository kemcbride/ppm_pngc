#!/usr/bin/env python

# 
# 

from __future__ import print_function

import argparse
import gzip
import sys
import os

from util import ZIP_PATH, parse_faa, parse_lastcol, parse_distances_file


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


def get_neighbor_proteins(protein_data, protein_id_locations):
    """Given protein data and relevant id locations,
    produce a set of WP ids to look up in order to 'get respective sequences'
    """
    # Proper solution probably involves using a set on i, or a dictionary
    # It would be nice to have a way to denote the 'match' proteins in each set.
    neighbor_ids = {}
    for name, pos in protein_id_locations.items():
        for i in range(pos-dist, pos):
            neighbor_ids[i] = protein_data[i]
        for i in range(pos+1, pos+1+dist):
            neighbor_ids[i] = protein_data[i]
    return neighbor_ids


def get_respective_sequences(gcf_id, neighbor_protein_ids):
    # We know that all sequences MUST be present in this given fasta file.
    faa_data = parse_faa('/'.join([ZIP_PATH, gcf_id+'.faa.gz']))
    neighbor_seqs = {k: faa_data[v] for k, v in neighbor_protein_ids.items()}
    return neighbor_seqs


def produce_family_matches(gcf_id, neighbor_sequences):
    """In order to match a 'family' to each sequence identified, we have a few steps:
    - produce a faa file with 1 sequence contained
    - run hmmscan --domtblout on that one sequence
    - take the first result and get the family from that - our 'family' value
    - delete leftover files / be careful not to waste space
    """
    # Use made up name for temporary direcotry - we MUST delete everything within it
    # by the end of this function. (Just being considerate, really.)
    output_data = {}
    dirname = 'functions-temp/{}'.format(gcf_id)
    os.mkdirs(dirname)
    # Now, we want to write the fasta sequences here - one per file.
    for pos, faa_data in neighbor_sequences.items():
        faa_file_name = dirname + '/' + faa_data.wp + '.faa'
        out_file_name = dirname + '/' + faa_data.wp + '.out'
        with open(faa_file_name, 'w') as faa_file:
            write_fasta_sequence(faa_data, output_file=faa_file)

        # Now, we want to run hmmscan on each of these and - in memory, get the family id
        run_hmmscan(faa_file_name, out_file_name)
        family = get_family(out_file_name)

        # Now we need to delete the files we've created.
        os.remove(faa_file_name)
        os.remove(out_file_name)

        # Now, we need to create some kind of output data structure.
        output_data[pos] = (family, faa_data)
    os.removedirs(dirname)
    return output_data


def print_family_data(gcf_id, protein_ids, output_data):
    # the style of output I want is "pos / WP / family / source(bool)"
    print('# {} - Family Annotations'.format(gcf_id))
    print('\t'.join(['POS', 'WP_ID', 'FAMILY', 'SOURCE']))
    for pos, data in output_data.items():
        # NOTE/TODO: could be nice to make the 'wp in keys' lookup constant time.
        print('\t'.join([
            pos,
            data[1].wp,
            data[0],
            str(data[1].wp in protein_ids.keys())
            ]))


def print_gcf_family_data(gcf_id, match_data, dist):
    # Go through the steps needed to produce output
    # TODO/NOTE: need to use match data/match file instead of "protein_ids"
    protein_data, protein_id_locations = collect_protein_data(gcf_id, protein_ids)
    neighbor_protein_ids = get_neighbor_proteins(protein_data, protein_id_locations)
    neighbor_sequences = get_respective_sequences(gcf_id, neighbor_protein_ids)
    family_match_data = produce_family_matches(gcf_id, neighbor_sequences)

    print_family_data(gcf_id, protein_id_locations, family_match_data)


# We need a new main, that will parse out the matches PER gcf, and run our old main on each
def main(distances_path, dist):
    gcf_data_sets = parse_distances_file(distances_path)
    for gcf_id, match_data in gcf_data_sets.items():
        print_gcf_family_data(gcf_id, match_data, dist)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=
            """Prints the 'product' value for proteins surrounding given protein ids in a given GCF.
            We understand this to be some form of functional annotation.\n
            output: position | id | product
            """)
    parser.add_argument('distances_file', help='Path to distances file you want to use as input to find functions for the contained matches.')
    parser.add_argument('--dist', default=10, type=int,
            help='The distance about each to produce annotations for, eg. 10')
    args = parser.parse_args()

    main(args.distances_file, args.dist)

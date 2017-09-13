#!/usr/bin/env python

# 
# 

from __future__ import print_function

import argparse
import gzip
import sys
import os

from util import (
        ZIP_PATH,
        parse_faa,
        parse_lastcol,
        parse_distances_file,
        write_fasta_sequence,
        )


GFF_PATH = '/research/gmh/GENOME_DB/gff-Complete'



class MatchLocation(object):
    def __init__(self, match, left_pos, right_pos):
        self.match = match
        self.right_pos = right_pos
        self.left_pos = left_pos


class MatchNeighborhood(object):
    def __init__(self, match_loc, dist):
        self.match_loc = match_loc
        self.dist = dist
        if match_loc.left_pos < match_loc.right_pos:
            self.range = range(match_loc.left_pos - dist, match_loc.right_pos + dist)
        else:
            self.range = range(match_loc.right_pos - dist, match_loc.left_pos + dist)
        self.ids = {} # structure:  pos : id


class GFFProteinData(object):
    def __init__(self, row_data):
        last_col_data = parse_lastcol(row_data[-1]) # call it on the last column

        # It's possible that we get a 'pseudo' protein, which has no Name/protein_id
        # TODO/NOTE: if a protein_data has '' for a name, we should wonder what's up with that...
        self.name = last_col_data.get('Name', last_col_data.get('protein_id', ''))
        # And Parent has the form 'gene12345', so we convert it to 12345
        self.parent = convert_str_int(last_col_data['Parent'])
        self.id = last_col_data['ID']
        self.product = last_col_data['product']
        self.note = last_col_data.get('Note', '')


def convert_str_int(string):
    filtered = filter(lambda s: s.isdigit(), string)
    return int(filtered)


def collect_protein_data(gcf_id, match_data):
    """ rewrite
    """
    protein_data = {}
    with gzip.open('/'.join([GFF_PATH, gcf_id]) + '.gff.gz', 'r') as f:
        for l in f:
            if l.startswith('#'):
                continue # it's a comment line, ignore

            data = l.split('\t') # we only care about data[0], data[2], data[-1]
            if data[2].lower() not in ['protein', 'cds']:
                # possibly: gene, sequence_feature, ...
                continue
            curr_protein = GFFProteinData(data)
            protein_data[curr_protein.parent] = curr_protein

    wp_pos_map = {protein.name: protein.parent for _, protein in protein_data.items()}

    match_locations = []
    for match in match_data:
        match_locations.append(MatchLocation(
            match,
            wp_pos_map[match.left_wp],
            wp_pos_map[match.right_wp]
            ))

    return protein_data, match_locations, wp_pos_map


def get_neighbor_ids_and_neighborhoods(protein_data, match_locations, dist):
    """Given protein data and relevant id locations,
    produce a set of WP ids to look up in order to 'get respective sequences'
    """
    match_neighborhoods = []
    neighbor_ids = set()
    for match_loc in match_locations:
        neighborhood = MatchNeighborhood(match_loc, dist)

        for pos in neighborhood.range:
            if pos in protein_data:
                neighbor_ids.add(protein_data[pos].name)
                neighborhood.ids[pos] = protein_data[pos].name
        match_neighborhoods.append(neighborhood)
    return neighbor_ids, match_neighborhoods


def get_respective_sequences(gcf_id, neighbor_ids):
    # We know that all sequences MUST be present in this given fasta file.
    faa_data = parse_faa('/'.join([ZIP_PATH, gcf_id+'.faa.gz']))
    # TODO/NOTE: why WOULDN'T the WP id be present in the faa_data?
    # - one case is '' protein, where we failed to find a name/WP_id
    neighbor_seqs = {wp_id: faa_data[wp_id] for wp_id in neighbor_ids if wp_id in faa_data}
    return neighbor_seqs


def extract_family(dirname, faa_data):
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
    return family


def get_families(gcf_id, neighbor_sequences):
    """Use extract_family to match up WP_ids to their respective
    families via hmmscan
    """
    # Use made up name for temporary direcotry - we MUST delete everything within it
    # by the end of this function. (Just being considerate, really.)
    output_data = {}
    dirname = 'functions-temp/{}'.format(gcf_id)
    os.makedirs(dirname)
    # Now, we want to write the fasta sequences here - one per file.
    for wp_id, faa_data in neighbor_sequences.items():
        family = extract_family(dirname, faa_data)
        output_data[pos] = (family, faa_data)
    os.removedirs(dirname)
    return output_data


def print_family_data(gcf_id, neighborhoods, family_data):
    # the style of output I want is "GCF / pos / WP / family / source(bool) / match_id"
    print('\t'.join(['GCF', 'POS', 'WP_ID', 'FAMILY', 'SOURCE', 'MATCH_NEIGHBORHOOD']))
    for neighborhood in neighborhoods:
        for i in neighborhood.range:
            print('\t'.join([
                gcf_id,
                i,
                data[1].wp,
                data[0],
                str(data[1].wp in [neighborhood.match.left_wp, neighborhood.match.right_wp]),
                str(neighborhood.match.left_wp + '+' + neighborhood.match.right_wp),
                ]))


def print_gcf_family_data(gcf_id, match_data, dist):
    # TODO/NOTE: need to use match data/match file instead of "protein_ids"
    protein_data, match_locations, wp_pos_map = collect_protein_data(gcf_id, match_data)
    neighbor_ids, neighborhoods = get_neighbor_ids_and_neighborhoods(protein_data, match_locations, dist)
    neighbor_sequences = get_respective_sequences(gcf_id, neighbor_ids)
    families = get_families(gcf_id, neighbor_sequences)

    print_family_data(gcf_id, neighborhoods, family_match_data)


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

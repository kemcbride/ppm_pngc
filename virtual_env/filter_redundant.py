#WRITING A CODE TO FILTER OUT REDUNDANT GENOMES, AND TO NARROW RESULTS DOWN TO THOSE THAT HAVE PPM MOTIF FROM OUR LISTS

import os, re

LIST = '/Vagabundo/monica/notes/complete_list'
PPM = '/Vagabundo/monica/notes/virtual_env/ppm_motif_complete_matches'
OUTPUT = '/Vagabundo/monica/notes/virtual_env/shortest_complete_distance'

#OUTPUT FORMAT is GCF_000005225 for LIST_PATH and GCF_000005225.out for PPM_PATH

if __name__ == '__main__':
    with open(LIST, "r") as ins:


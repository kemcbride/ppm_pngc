#WRITING A CODE TO FILTER OUT REDUNDANT GENOMES, AND TO NARROW RESULTS DOWN TO THOSE THAT HAVE PPM MOTIF FROM OUR LISTS

import os, re

LIST = '/Vagabundo/monica/notes/complete_list'
PPM = '/Vagabundo/monica/notes/virtual_env/ppm_motif_complete_matches'
OUTPUT = '/Vagabundo/monica/notes/virtual_env/shortest_complete_distance'

#OUTPUT FORMAT is GCF_000005225 for LIST_PATH and GCF_000005225.out for PPM_PATH
if __name__ == '__main__':
    with open(LIST) as f:
       nonredundant_list = f.readlines()
    with open(PPM) as f:
       ppm_motif_list = f.readlines()
    with open(OUTPUT) as f:
        line = f.readline()
        while line != '':
             if line.startswith('#'):
                fname = line.split()[-1]
                ppm_name = line.split()[-1]
                nonredundant_name = ppm_name.split('.')[0]
                if motif_name in motif_list and nonredundant_name in nonredundant_list:
                   print()



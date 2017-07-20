#!/usr/bin/env python

# We need to count/get the distinct WP ids from outfiles with both PNGC/PPM, and just PNGC.
#
# We can get a ratio of Total :: JustPngC
# Run this program like so:
# python unique_query_names.py <path to matches.txt file> <path to directory of .out files>

from __future__ import print_function
from argparse import ArgumentParser
import pandas as pd
import os

COLS = ['target_name', 't_accession', 'tlen', 'query_name']


def unique_query_names(list_path, outfile_dir):
    with open(list_path, 'r') as listfile:
        outfile_list = listfile.readlines()
    outfile_list = [fname.strip() for fname in outfile_list]

    unique_query_names = set()
    for filename in outfile_list:
        fpath = os.path.join(outfile_dir, filename)
        if not os.path.isfile(fpath):
            print ("File {} does not exist".format(fpath))
            continue
        df = pd.read_csv(os.path.join(outfile_dir, filename),
            comment='#',
            header=None,
            delimiter='\s+',
            usecols=range(4)
            )
        df.columns = COLS  # add the column names to the dataframe
        unique_query_names.update(df.query_name.values())
    return unique_query_names

    
if __name__ == '__main__':
    ap = ArgumentParser()
    ap.add_argument('list_path')
    ap.add_argument('outfile_dir')
    args = ap.parse_args()
    unique_names = unique_query_names(args.list_path, args.outfile_dir)
    print( len(unique_names) )
    for item in unique_names:
        print( item )

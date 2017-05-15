# Some kind of program to measure distances between proteins that have hits
# based on a certain .out file (70-cutga-out output.)

import pandas as pd
import sys
import os
from StringIO import StringIO
from collections import defaultdict
import itertools
import gzip

# Lists of the protein names that are PNGC or PEP_MUTASE
from util import PPM_MATCH_LIST, PNGC_MATCH_LIST

# I'm assuming that we'll just know what the file paths we want are:
INPUT_PATH = '/Vagabundo/monica/temp/70-CUTGA-OUT-faa-Complete'
GFF_PATH = '/research/gmh/GENOME_DB/gff-Complete'
COLS = ['target_name', 't_accession', 'tlen', 'query_name', 'q_accession', 'qlen', 'e_full', 'score_full', 'bias_full', 'num_domain', 'of_domain', 'ie_domain', 'score_domain', 'bias_domain', 'from_hmm', 'to_hmm', 'from_ali', 'to_ali', 'from_env', 'to_env', 'acc', 'desc_target']

# This format string a - calls str() on each arg, and then
# formats it left-aligned with 20 spaces
DATA_FMT = '{!s:<20}'
        
def parse_lastcol(col_text):
    data = col_text.split(';')
    data = [segment.split('=') for segment in data]
    data = {segment[0]:segment[1] for segment in data}
    # Breaking the last column of the GFF file into segments that we care about i.e. Parent=gene14
    # So now data = {key:value for a=b;y=z} in col_text
    # We've turned col_text into a python dictionary
    return data


def main():
    # Print the column names
    print(''.join([DATA_FMT.format(colname) for colname in ['PPM', 'PNGC', 'Distance']]))

    # Reading all the .out files in gff-Complete
    out_files = os.listdir(INPUT_PATH)
    print >> sys.stderr, '\n'.join(out_files)
    for fname in out_files:
        try:
            input_df = pd.read_csv(os.path.join(INPUT_PATH, fname),
                comment='#',
                header=None,
                delimiter='\s+',
                usecols=range(22)
                )
            input_df.columns = COLS  # add the column names to the dataframe
            query_names = input_df[['query_name', 'target_name']]
            # We are matching the target name (i.e. PEP_mutase) with the query name (i.e. WP_05123543.01) 

            # Now we're only creating pairs between PPM and PNGC rows.
            pep_mutase_rows = query_names.loc[
                    query_names['target_name'].isin(PPM_MATCH_LIST)]
            pngc_rows = query_names.loc[
                    query_names['target_name'].isin(PNGC_MATCH_LIST)]
            # Dataframe.values turns it into a list of lists, where each inner list is the row of data
            # itertools.product does the cartesian product between the two sets.
            pairs = itertools.product(pep_mutase_rows.values, pngc_rows.values)

            # Load the GFF file
            gff_fname = fname.split('.')[0] + '.gff.gz'
            with gzip.GzipFile(os.path.join(GFF_PATH, gff_fname), 'r') as gff_file:
                gff_text = gff_file.read()
            gff_lines = [line for line in gff_text.split('\n') if 'CDS' in line]
            gff_lines = [line for line in gff_text.split('\n') if 'Protein' in line]
            gff_lines = [line for line in gff_lines if 'protein_id' in line]
            gff_text = '\n'.join(gff_lines)
            # Only interested in the lines containing gene# and the protein_id
            gff_df = pd.read_csv(StringIO(gff_text),
                comment='#',
                header=None,
                delimiter='\t',
                )

        except Exception as e:
            # idea: make all output that is non-relevant/useful a comment,
            # so any code that reads the output files can ignore it
            print('# {}'.format(e))
            continue

        last_col = gff_df[gff_df.columns[-1]]
        # this function turns the last column into a dictionary from its tag=value; format
        last_col_data = last_col.apply(parse_lastcol)

        # turn the keys/values from last_col into additional columns on gff_df
        new_dict = defaultdict(list)
        first_keys = []
        for el in last_col_data:
            if el == last_col_data[0]:
                first_keys = el.keys()
            for k, v in el.items():
                new_dict[k].append(v)

        parent_genes, protein_ids = new_dict['Parent'], new_dict['Name']
        # This is a "fancy/annoying" way of turning string -> int, eg. "gene6" -> 6
        gff_df['parent_gene'] = [int(filter(type(seq).isdigit, seq)) for seq in parent_genes]
        gff_df['protein_id'] = protein_ids
        # 'Parent' gives us the gene number we care about and 'Name' is the identifier
        # Then, we find the ones where the 'name' from last_col data == 'WP_xyz'
        # And then for thsoe columns, we subtract for the corresponding pairs
        # to get the gene distance using pairs

        results = []
        for pair in pairs:
            # We select 'left' and 'right' pair members based on whether they match
            # the target and query names of the pair that we have from input_df
            try:
                left = gff_df[gff_df.protein_id == pair[0][0]].iloc[0]
                right = gff_df[gff_df.protein_id == pair[1][0]].iloc[0]
            except:
                # print >> sys.stderr, 'HIT EXCEPTION IN PAIR PRINTING LOOP'
                continue

            # '0' here referring to the first column in the GFF file - the DNA/gene id
            # eg: NZ_xyzxyzxyz - if these are not from the same gene/DNA, skip it
            if left[0] != right[0]:
                continue

            distance = abs(left.parent_gene - right.parent_gene)
            results.append((pair, distance))

        print('# {}'.format(fname))
        for r in results:
            # r looks like ( (left_data, right_data), distance )
            left = r[0][0]
            right = r[0][1]
            distance = r[1]
            print(''.join(
                [DATA_FMT.format(x) for x in [left[0], right[0], distance]]
            ))

if __name__ == '__main__':
    main()

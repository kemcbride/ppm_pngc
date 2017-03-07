# Some kind of program to measure distances between proteins that have hits
# based on a certain .out file (70-cutga-out output.)
#

# Do your imports
import pandas as pd
import os
from StringIO import StringIO
from collections import defaultdict
import itertools
import zipfile

# I'm assuming that we'll just know what the file paths we want are:
INPUT_PATH = '/Vagabundo/monica/temp/70-CUTGA-OUT-faa-Complete'
GFF_PATH = '/research/gmh/GENOME_DB/gff-Complete'
COLS = ['target_name', 't_accession', 'tlen', 'query_name', 'q_accession', 'qlen', 'e_full', 'score_full', 'bias_full', 'num_domain', 'of_domain', 'ie_domain', 'score_domain', 'bias_domain', 'from_hmm', 'to_hmm', 'from_ali', 'to_ali', 'from_env', 'to_env', 'acc', 'desc_target']

GFF_PATH = '/home/kelly/Dropbox/Stuff/MSc/GCF_000439255.gff'
INPUT_PATH = '/home/kelly/Dropbox/Stuff/MSc/GCF_000439255.out'

# Reading all the zipped files in gff-Complete
for fname in os.listdir(GFF_PATH):
    with zipfile.Zipfile(os.path.join(GFF_PATH, fname, 'r')):
        
def parse_lastcol(col_text):
    data = col_text.split(';')
    data = [segment.split('=') for segment in data]
    data = {segment[0]:segment[1] for segment in data}
    return data


def main():
    # Load the input file as a dataframe.
    input_df = pd.read_csv(os.path.join(INPUT_PATH),
        comment='#',
        header=None,
        delimiter='\s+',
        usecols=range(22)
        )
    input_df.columns = COLS
    # We only care about the 'query_name' column - this might not be the right way of getting it, whatever. Idk.
    query_names = input_df[['query_name', 'target_name']]

    permutations = itertools.permutations(query_names.values, 2)

    # Load the GFF file
    with open(GFF_PATH, 'r') as gff_file:
        gff_text = gff_file.read()
    gff_lines = [line for line in gff_text.split('\n') if 'Protein' in line]
    gff_lines = [line for line in gff_text.split('\n') if 'protein_id' in line]
    gff_text = '\n'.join(gff_lines)

    gff_df = pd.read_csv(StringIO(gff_text),
        comment='#',
        header=None,
        delimiter='\t',
        )

    last_col = gff_df[gff_df.columns[-1]]
    last_col_data = last_col.apply(parse_lastcol)

    new_dict = defaultdict(list)
    first_keys = []
    for el in last_col_data:
        if el == last_col_data[0]:
            first_keys = el.keys()
        for k, v in el.items():
            new_dict[k].append(v)

    parent_genes, protein_ids = new_dict['Parent'], new_dict['protein_id']
    gff_df['parent_gene'] = [int(filter(type(seq).isdigit, seq)) for seq in parent_genes]
    gff_df['protein_id'] = protein_ids

    results = []
    for pair in permutations:
        left = gff_df[gff_df.protein_id == pair[0][0]]
        right = gff_df[gff_df.protein_id == pair[1][0]]

        lpos = left.get_value(left.index[0], 'parent_gene')
        rpos = right.get_value(right.index[0], 'parent_gene')
        distance = abs(lpos - rpos)

        results.append((pair, distance))

    for r in results:
        if r[0][1][1] != r[0][0][1]:
            print(r[0][1][1], r[0][0][1], r[1])

if __name__ == '__main__':
    main()

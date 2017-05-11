# Some kind of program to measure distances between proteins that have hits
# based on a certain .out file (70-cutga-out output.)
#

# Do your imports
import pandas as pd
import os
from StringIO import StringIO
from collections import defaultdict
import itertools
import gzip
import ipdb

# I'm assuming that we'll just know what the file paths we want are:
INPUT_PATH = '/Vagabundo/monica/temp/70-CUTGA-OUT-faa-Complete'
GFF_PATH = '/research/gmh/GENOME_DB/gff-Complete'
COLS = ['target_name', 't_accession', 'tlen', 'query_name', 'q_accession', 'qlen', 'e_full', 'score_full', 'bias_full', 'num_domain', 'of_domain', 'ie_domain', 'score_domain', 'bias_domain', 'from_hmm', 'to_hmm', 'from_ali', 'to_ali', 'from_env', 'to_env', 'acc', 'desc_target']

        
def parse_lastcol(col_text):
    data = col_text.split(';')
    data = [segment.split('=') for segment in data]
    data = {segment[0]:segment[1] for segment in data}
    # Breaking the last column of the GFF file into segments that we care about i.e. Parent=gene14
    # So now data = {key:value for a=b;y=z} in col_text
    # We've turned col_text into a python dictionary
    return data


def main():
    # Reading all the .out files in gff-Complete
    out_files = os.listdir(INPUT_PATH)
    for fname in out_files:
        try:
            input_df = pd.read_csv(os.path.join(INPUT_PATH, fname),
                comment='#',
                header=None,
                delimiter='\s+',
                usecols=range(22)
                )
            input_df.columns = COLS
            # We only care about the 'query_name' column - this might not be the right way of getting it, whatever. Idk.
            query_names = input_df[['query_name', 'target_name']]
    # We are matching the target name (i.e. PEP_mutase) with the query name (i.e. WP_05123543.01) 
            permutations = itertools.permutations(query_names.values, 2)

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
            print(e)
            continue

        last_col = gff_df[gff_df.columns[-1]]
        last_col_data = last_col.apply(parse_lastcol)
        # Only looking at the last column in the GFF file
        # Then... ideally, you'd turn the keys/values from last_col into additional columns on gff_df

        new_dict = defaultdict(list)
        first_keys = []
        for el in last_col_data:
            if el == last_col_data[0]:
                first_keys = el.keys()
            for k, v in el.items():
                new_dict[k].append(v)

        parent_genes, protein_ids = new_dict['Parent'], new_dict['Name']
        gff_df['parent_gene'] = [int(filter(type(seq).isdigit, seq)) for seq in parent_genes]
        gff_df['protein_id'] = protein_ids
    # 'Parent' gives us the gene number we care about and 'Name' is the identifier
    # Then, we find the ones where the 'name' from last_col data == 'WP_xyz'
    # And then for thsoe columns, we subtract for the corresponding pairs to get the gene distance using permutations

        results = []
        for pair in permutations:
            # We select 'left' and 'right' pair members based on whether they match
            # the target and query names of the pair that we have from input_df
            left = gff_df[gff_df.protein_id == pair[0][0]]
            right = gff_df[gff_df.protein_id == pair[1][0]]

            # Here we do .values[0][0] - values means we care only about the values,
            # not about the "dataframe", the first [0] indicates the first row of the
            # dataframe, and the second [0] indicates the first column (the DNA name?)
            if left.values[0][0] != right.values[0][0]:
                continue

            lpos = left.get_value(left.index[0], 'parent_gene')
            rpos = right.get_value(right.index[0], 'parent_gene')
            distance = abs(lpos - rpos)
            if distance == 0:
                continue
            results.append((pair, distance))

        print(fname)
        print('\n')
        for r in results:
            # r looks like ( (left_dataframe, right_dataframe), distance )
            left = r[0][0]
            right = r[0][1]
            distance = r[1]
            # Left[1], Right[1]
            if left[1] != right[1]:
                print(','.join([str(left[1]), str(left[0]), str(right[1]), str(right[0]), str(distance)]))
        print('\n')

if __name__ == '__main__':
    main()

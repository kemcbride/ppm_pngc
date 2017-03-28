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
    return data


def main():
    # Reading all the .out files in gff-Complete
    out_files = os.listdir(INPUT_PATH)
    for fname in out_files:
        try:
            # Load the input file as a dataframe
            input_fname = fname.split('.')[0] + '.out'
            #DIS DA PART TO IGNORE UNMATCHED GFF FILES
            input_df = pd.read_csv(os.path.join(os.path.join(INPUT_PATH, input_fname)),
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
            gff_fname = fname.split('.')[0] + '.gff.gz'
            with gzip.GzipFile(os.path.join(GFF_PATH, gff_fname), 'r') as gff_file:
                gff_text = gff_file.read()
            gff_lines = [line for line in gff_text.split('\n') if 'Protein' in line]
            gff_lines = [line for line in gff_lines if 'protein_id' in line]
            gff_text = '\n'.join(gff_lines)
            
            # Okay so this is super tricky I realize but this is what we need to figure out: 
            # So we open up the gff file, but we aren't only concerned with the "Protein" line
            # Basically, we need to look at the first column (NZ_CP01...) and AFTER we locate where the two genes are,
            # we need to check that they belong to the same "DNA" which is indicated by the first column identifier
            # So it's like Step #1 use Name= in that line to find both the genes
            # Step #2 find the identifier (from column 1) for gene #1 and then gene #2
            # If identifiers are = then continue process, find distances, record that shiz
            # If identifiers =/= then skip these genes because they are not from the same DNA, we don't care about it
            # The tricky part is that we can't do that line split thing with only looking at the lines with "Protein" in it because it's every single line we check
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
                print('{} {} {}'.format(left[1], right[1], distance))
        print('\n')

if __name__ == '__main__':
    main()

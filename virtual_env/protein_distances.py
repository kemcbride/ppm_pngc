# Some kind of program to measure distances between proteins that have hits
# based on a certain .out file (70-cutga-out output.)
#

# Do your imports
import pandas as pd
import ipdb # This is a cool import - if you put the line ipdb.set_trace() in your code,
# it will automatically open a debugger for you when you run it. then you can type in your variables,
# and see what their values are. IT's really useful. you can quit with 'q' and continue the program with 'c' or next line with 'n'

# I'm assuming that we'll just know what the file paths we want are:
INPUT_PATH = '/Vagabundo/monica/temp/70-CUTGA-OUT-faa-Complete'
GFF_PATH = '/research/gmh/GENOME_DB/gff-Complete'
COLS = ['target_name', 't_accession', 'tlen', 'query_name', 'q_accession', 'qlen', 'e_full', 'score_full', 'bias_full', 'num_domain', 'of_domain', 'ie_domain', 'score_domain', 'bias_domain', 'from_hmm', 'to_hmm', 'from_ali', 'to_ali', 'from_env', 'to_env', 'acc', 'desc_target']

def parse_lastcol(col_text):
    data = col_text.split(';')
    data = [segment.split('=') for segment in data]
    data = {segment[0]:segment[1] for segment in data}
    # So now data = {key:value for a=b;y=z in col_text
    # We've turned col_text into a python dictionary
    return data

def main():
    # Load the input file as a dataframe.
    input_df = pd.read_csv(os.path.join(INPUT_PATH),
        comment='#',
        header=None,
        delimiter='\s+',
        usecols=range(23)
        )
    # I'm not sure if this is right (or useful)
    input_df.columns = COLS
    # We only care about the 'query_name' column - this might not be the right way of getting it, whatever. Idk.
    query_names = input_df['query_name']

    # Somewhat tricky - create a list of tuples/lists of each permutation/pair in query_names
    # something like,
    permutations = []
    for item in query_names:
        for other_item in query_names:
            if item != other_item:
                permutations += (item, other_item)

    # Then we need to load the um, gff file.
    # If it loads and is weird, you might need to set usecols or change delimiter to maybe '\t'
    gff_df = pd.read_csv(os.path.join(GFF_PATH),
        comment='#',
        header=None,
        delimiter='\s+',
        )

    # The second column, we only care about the rows where the second column = 'Protein'
    # This is 100% wrong sample code, but it should look something like this... ish.
    gff_df = gff_df['2' == 'Protein']

    # Uhh so now you kinda wanna parse out the parent gene number
    # you need to do this from the last columns. so like, idk.
    # You'll need to actually figure out how to get last_col to be what you want
    last_col = gff_df[-1]
    last_col.apply(parse_lastcol)

    # Then... ideally, you'd turn the keys/values from last_col into additional columns on gff_df

    # Then, you'd find the ones where the 'name' from last_col data == 'WP_xyz' that you want
    # And then for thsoe columns, you'd do the subtraction for the corresponding pairs
    # in the permutations list


if __name__ == '__main__':
    main()

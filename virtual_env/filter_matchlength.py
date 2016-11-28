# coding: utf-8
import os, pandas as pd

# We're choosing 70% for starters, this may be too low.
MATCH_SIZE_THRESHOLD = 0.7
SETNAME = 'Complete'

# This variable isn't really "useful" but like, it's also not "useless"
COLS = ['target_name', 't_accession', 'tlen', 'query_name', 'q_accession', 'qlen', 'e_full', 'score_full', 'bias_full', 'num_domain', 'of_domain', 'ie_domain', 'score_domain', 'bias_domain', 'from_hmm', 'to_hmm', 'from_ali', 'to_ali', 'from_env', 'to_env', 'acc', 'desc_target']

TEMP = '/Vagabundo/monica/TEMP/'
IDIR = "CUTGA-OUT-faa-{}".format(SETNAME)
ODIR = "{}-CUTGA-OUT-faa-{}".format(str(MATCH_SIZE_THRESHOLD*100).split('.')[0], SETNAME)

# This is a dumb function to make writing the output easier...
def offset_idx(a, idx):
    c = 0
    for l in a:
        if l[0] != '#':
            break
        c+=1
    idx += c
    return idx


# Let's write a function that will do all the stuff I want to do: that is, read new file into df. create match_len column and filter df. write out new file based on df.
def write_filtered_file(fname):
    tfile = open(os.path.join(os.path.join(TEMP, IDIR, fname)), 'r')
    tfcontents = tfile.read()
    lines = tfcontents.split('\n')

    # Read the actual table contents into a pandas dataframe
    tdf = pd.read_csv(os.path.join(TEMP,IDIR, fname),
            comment='#', # Did you know pandas can automatically skip comment lines???
            header=None, # the header domtblout provides is virtually useless
            delimiter='\s+', # Unfortunately, domtblout is an idiot and chose NOT to use tabs to separate columns
            usecols=range(23) # Prevent it from freaking out about the description column (which is basically garbage anyways)
            )

    # Identify rows to be filtered out
    # 16 = from hmm. 15 = to hmm. 2 = tlen see: variable COLS
    # Worth noting that I'm trusting pandas to handle float division correctly... I belive it's doing its job.
    tdf['match_len'] = (tdf[16] - tdf[15]) / tdf[2]
    bad_row_list = tdf[tdf['match_len'] > MATCH_SIZE_THRESHOLD].index

    # Offset the "bad rows" by the size of the header comments - we're keeping that metadata bc why not.
    offset_bad_rows = offset_idx(lines, bad_row_list)
    filtered_lines = [ i for j, i in enumerate(lines) if j not in offset_bad_rows ]

    # Write the new file
    with open(os.path.join(TEMP, ODIR, fname), 'w') as new_file:
        new_file.write('\n'.join(filtered_lines))
    # And that's it!


if __name__ == '__main__':
    completed_files = os.listdir(os.path.join(TEMP, ODIR))
    input_files = os.listdir(os.path.join(TEMP, IDIR))
    input_files = [f for f in input_files if f not in completed_files]

    # Note! There's no like, responsive output yet. NOR is there parallelization!
    # But it does work. Tested it on those 3 first files.
    for f in input_files:
        write_filtered_file(f)

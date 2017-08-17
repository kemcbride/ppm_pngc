import os, pandas as pd

# Constants for now - From process_hmmscan
SETNAME = 'Complete'
DB_NAME = 'faa-{}'.format(SETNAME)
TEMP_PATH ='/Vagabundo/monica/temp/'
# TEMP_PATH ='/home/kelly/Dropbox/gff/temp/'
OUTPUT_DIR ='CUTGA-OUT-{}'.format(DB_NAME)
ZIP_PATH = os.path.join('/research/gmh/GENOME_DB/{}', DB_NAME)

# From/for sort_matches.py and protein_distances.py
PPM_MATCH_LIST =['PEP_mutase']
PNGC_MATCH_LIST =['NTP_transferase', 'NTP_transf_3', 'IspD']

HMM_FILE = '/Vagabundo/monica/Proteins/models.hmm'

MULTIPROCESSING_FACTOR = 100 # We'll run 100 per batch


class FastaData(object):

    def __init__(self, wp_id, sequence_data, extra_data=''):
        self.wp = wp_id
        self.sequence = sequence_data
        self.extra_data = extra_data


# You could pass in a lambda/function to cond and have that also filter out
# false results from idir (ie. is file in idir also in sort_matches both.txt)
def files_remaining(idir, odir, cond=lambda x: true):
    """Pass None if you want no filtering to happen...
    Although the default value is also a pass-thru
    """
    files = [f for f in os.listdir(idir) if f not in os.listdir(odir)]
    if cond is not None:
        files = [f for f in files if cond(f)]
    return files


BLASTDB = '/research/gmh/GENOME_DB/blastpDB-Complete/GCF_{}'
def blast_motif_match(qname):
    command = 'blastdbcmd -db ' + BLASTDB + ' -target_only -entry ' + qname + ' -outfmt %s'
    output = subprocess.check_output(command, shell = True)
    has_motif = re.search('EDK\w{3,7}NS',output)
    return bool(has_motif)


def parse_faa(faa_path):
    with open(faa_path, 'r') as f:
        file_contents = f.read()
    raw_sequence_list = file_contents.split('>')

    del file_contents # i hope this wouln't logically break the code

    # Now, I need to turn each item into the form WP_id : sequence
    faa_data = {}

    # We know that the WP_id is the first space delimited part of the first line.
    for raw_seq in raw_sequence_list:
        if not raw_seq: # we know for sure the first result will be empty
            continue
        seq_lines = raw_seq.split('\n')
        first_line_data = seq_lines[0].split()
        wp_id = first_line_data[0]
        extra_data = ' '.join(seq_lines[0].split()[1:]) if len(first_line_data) > 1 else ''
        sequence = ''.join(seq_lines[1:])
        fd = FastaData(wp_id, sequence, extra_data=extra_data)
        faa_data[wp_id] = fd

    return faa_data

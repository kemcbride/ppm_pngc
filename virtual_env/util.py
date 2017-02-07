import os, pandas as pd

# Constants for now - From process_hmmscan
SETNAME = 'Complete'
DB_NAME = 'faa-{}'.format(SETNAME)
TEMP_PATH ='/Vagabundo/monica/temp/'
OUTPUT_DIR ='CUTGA-OUT-{}'.format(DB_NAME)
ZIP_PATH = os.path.join('/research/gmh/GENOME_DB/{}', DB_NAME)

HMM_FILE = '/Vagabundo/monica/Proteins/models.hmm'

MULTIPROCESSING_FACTOR = 100 # We'll run 100 per batch

# You could pass in a lambda/function to cond and have that also filter out
# false results from idir (ie. is file in idir also in sort_matches both.txt)
def files_remaining(idir, odir, cond=lambda x: true):
    files = [f for f in os.listdir(idir) if f not in os.listdir(odir)]
    if cond is not None:
        files = [f for f in files if cond(f)]
    return files
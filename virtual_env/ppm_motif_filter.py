import os, re

database = "Complete"
#database = "Scaffold"

FAA_PATH = '/Vagabundo/monica/temp/faa-' + database
INPUT_PATH = '/Vagabundo/monica/temp/70-CUTGA-OUT-faa-' + database

#Check filenames in 70-CUTGA (.out) and scan for those filenames in faa-Complete(.faa)
#Read filenames in faa-Complete (amino acid sequences) to scan for conserved motif (EDKXXXXXNS)
#The XXXXX is a variable region where X can be any letter
#The variable region itself can vary in number of X - scan for 3 - 7X
#Output should be all of the filenames that have this motif from faa-Complete

#def multi_match(target, patterns):
#    """Return True if target begins with any of the strings in patterns"""
#    for pattern in patterns:
#        if target.startswith(pattern):
#            return True

#def multi_search(target,patterns):
#    """Return the first position in target where any of the strings in patterns is found"""
#    return min([target.find(pattern)
#                for pattern in patterns
#                if target.find(pattern >= 0)])

def file_has_motif(FAA_PATH, motif_regex):
    file_contents = open(FAA_PATH, 'r').read()
    has_motif = bool(re.search(motif_regex, file_contents))
    return has_motif

if __name__ == '__main__':
#read .out files to get the names that we want from .faa
    out_files = os.listdir(INPUT_PATH)
    for fname in out_files:
        faa_fname = fname.split('.')[0] + '.faa'
        #read .faa files?? probably?
        has_motif = file_has_motif(FAA_PATH + '/' + faa_fname, 'EDK\w{3,7}NS')
        if has_motif:
           print fname

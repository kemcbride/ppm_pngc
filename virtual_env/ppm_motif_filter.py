INPUT_PATH = /Vagabundo/monica/temp/faa-Complete/
DATA_PATH = /Vagabundo/monica/temp/70-CUTGA-faa-Complete/
#Check filenames in 70-CUTGA (.out) and scan for those filenames in faa-Complete(.faa)
#Read filenames in faa-Complete (amino acid sequences) to scan for conserved motif (EDKXXXXXNS)
#The XXXXX is a variable region where X can be any letter
#The variable region itself can vary in number of X - scan for 3 - 7X
#Output should be all of the filenames that have this motif from faa-Complete

def multi_match(target, patterns):
    """Return True if target begins with any of the strings in patterns"""
    for pattern in patterns:
        if target.startswith(pattern):
            return True

def multi_search(target,patterns):
    """Return the first position in target where any of the strings in patterns is found"""
    return min([target .find(pattern)
                for pattern in patterns
                if target.find(pattern >= 0])
     
multi_match(sequence, ('EDK,\D, \D, \D, \D, \D, NS'))
multi_search(sequence, ('EDK,\D, \D, \D, \D, \D, NS'))

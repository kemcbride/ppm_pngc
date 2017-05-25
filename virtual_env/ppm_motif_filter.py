DATA_PATA = /Vagabundo/monica/temp/faa-Complete/

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

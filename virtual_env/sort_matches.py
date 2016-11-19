from __future__ import print_function
import os

DATAPATH = '/Vagabundo/monica/temp/CUTGA-OUT-faa-{}'
SETNAME = 'Complete' # Or Complete
OUTPUTPATH = '/Vagabundo/monica/Proteins/'
PEP_MUTASE_MATCH_LIST =['PEP_mutase']
PNGC_MATCH_LIST =['NTP_transferase', 'NTP_transf_3', 'IspD']


def sort_matches():
    pngc_list = []
    pepmutase_list = []
    count = 0

    for datafile in os.listdir(DATAPATH.format(SETNAME)):
        print("On file {}...".format(count), end='\r')

        matched_pepmutase = False
        matched_pngc = False
        matched_both = False
        with open(os.path.join(DATAPATH.format(SETNAME), datafile), 'r') as df:
            # Header section - skip!
            line = df.readline()
            while(line.startswith('#')):
                line = df.readline()

            # Data section! Do not skip!
            line = df.readline()
            while(not line.startswith('#') and not matched_both):
                # All praise be to domtblout format for making this easy for us...
                match_name = line.split()[0] # the first "word" of each line
                if not matched_pepmutase and match_name in PEP_MUTASE_MATCH_LIST:
                    pepmutase_list.append(datafile)
                    matched_pepmutase = True
                elif not matched_pngc and match_name in PNGC_MATCH_LIST:
                    pngc_list.append(datafile)
                    matched_pngc = True

                matched_both = matched_pngc and matched_pepmutase

                # Reading the next line goes at the end of the loop body
                line = df.readline()
        count += 1

    # Now we're done reading the files
    # So now let's do some set arithmetic to determine both, and exclusively one
    pngc_set = set(pngc_list)
    pepmutase_set = set(pepmutase_list)
    both = pngc_set & pepmutase_set
    exclusive_pepmutase = pepmutase_set - pngc_set
    exclusive_pngc = pngc_set - pepmutase_set

    for name, fileset in [('both', both), ('pepmutase_only', exclusive_pepmutase), ('pngc_only', exclusive_pngc)]:
        with open(os.path.join(OUTPUTPATH, '{}_{}_matches.txt'.format(name, SETNAME)), 'w') as out:
            for filename in fileset:
                out.write(filename + '\n')
    # Done~!!!

def main():
    sort_matches()
    print("Find your sorted filename results in {}".format(OUTPUTPATH))

if __name__ == '__main__':
    main()

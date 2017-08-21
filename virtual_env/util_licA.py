Last login: Fri Aug 18 15:48:14 on ttys002
xook:~ monica$ cd notes/virtual_env/
xook:virtual_env monica$ nano util.py















































  GNU nano 2.0.6                  File: util.py                                           

import os, pandas as pd

# Constants for now - From process_hmmscan
SETNAME = 'Complete'
DB_NAME = 'faa-{}'.format(SETNAME)
TEMP_PATH ='/Vagabundo/monica/temp/'
# TEMP_PATH ='/home/kelly/Dropbox/gff/temp/'
#OUTPUT_DIR ='CUTGA-OUT-{}'.format(DB_NAME)
OUTPUT_DIR = 'licA-licC/CUTGA-OUT-{}'.format(DB_NAME)
ZIP_PATH = os.path.join('/research/gmh/GENOME_DB/', DB_NAME)

# From/for sort_matches.py and protein_distances.py
PPM_MATCH_LIST =['PEP_mutase']
PNGC_MATCH_LIST =['NTP_transferase', 'NTP_transf_3', 'IspD']
LICA_MATCH_LIST = ['Choline_kinase']

#HMM_FILE = '/Vagabundo/monica/Proteins/models.hmm'
HMM_FILE = '/Vagabundo/monica/Proteins/licA-licC_models.hmm'

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
    command = 'blastdbcmd -db ' + BLASTDB + ' -target_only -entry ' + qname + ' -outfmt %$
                                    [ Read 73 lines ]
^G Get Help    ^O WriteOut    ^R Read File   ^Y Prev Page   ^K Cut Text    ^C Cur Pos
^X Exit        ^J Justify     ^W Where Is    ^V Next Page   ^U UnCut Text  ^T To Spell
  [Restored Aug 18, 2017, 4:14:13 PM]
Last login: Fri Aug 18 16:14:13 on ttys001
/sw/share/dbus/launchd/org.finkproject.dbus-session.plist: service already loaded
kuum:virtual_env monica$ ls
#protein_distances.py#*               ppm_complete_list
Clusters/                             ppm_motif_filter_update.py*
bin/                                  ppm_motif_filter_update_no_filter.py
compile_cluster_fasta.py              process_hmmscan.py*
compile_fasta                         protein_distances_no_ppm.py
compile_fasta_no_ppm                  protein_distances_ppm.py
complete_distances                    requirements.txt
complete_distances_2                  runner.py*
complete_distances_no_motif           scaffold_distances
complete_list_no_motif                scaffold_ppm_matches
complete_ppm_matches                  share/
complete_shortest_distance            shortest_complete_distance
complete_species                      shortest_complete_distances_2
contig_distances                      shortest_contig_distance
contig_ppm_matches                    shortest_distance.py*
directory/                            shortest_scaffold_distance
filter_empties.py*                    sort_matches.py*
filter_matchlength.py                 sort_matches_complete.py*
filter_redundant.py*                  sort_matches_scaffold.py*
functions.py                          species_identifier.py
include/                              unique_query_names.py
lib/                                  util.py
mkBlastClusters.pl*                   util.pyc
pip-selfcheck.json
kuum:virtual_env monica$ more compile_fasta
kuum:virtual_env monica$ cp compile_fasta ~/Desktop
kuum:virtual_env monica$ 
  [Restored Aug 21, 2017, 4:44:48 PM]
Last login: Mon Aug 21 16:44:47 on ttys000
kuum:virtual_env monica$ nano compile_cluster_fasta.py
kuum:virtual_env monica$ rm protein_distances_no_ppm.py
kuum:virtual_env monica$ nano protein_distances_no_ppm.py
kuum:virtual_env monica$ python protein_distances_no_ppm.py > complete_distances_no_motif
  File "protein_distances_no_ppm.py", line 133
    
    ^
IndentationError: expected an indented block
kuum:virtual_env monica$ nano protein_distances_no_ppm.py 
kuum:virtual_env monica$ python protein_distances_no_ppm.py > complete_distances_no_motif
  File "protein_distances_no_ppm.py", line 132
    
                              ^
IndentationError: expected an indented block
kuum:virtual_env monica$ nano protein_distances_no_ppm.py 
kuum:virtual_env monica$ nano util.py
kuum:virtual_env monica$ nano protein_distances_no_ppm.py 
kuum:virtual_env monica$ python protein_distances_no_ppm.py 
PPM                 PNGC                Distance            
# GCF_000005825.out
# No columns to parse from file
# GCF_000006605.out
WP_011273297.1      WP_005292098.1      1231                
WP_011273297.1      WP_005293089.1      723                 
WP_011273297.1      WP_005293380.1      835                 
WP_011273297.1      WP_005293642.1      868                 
WP_011273297.1      WP_005293920.1      986                 
WP_011273297.1      WP_011273020.1      402                 
WP_011274226.1      WP_005292098.1      52                  
WP_011274226.1      WP_005293089.1      560                 
WP_011274226.1      WP_005293380.1      448                 
WP_011274226.1      WP_005293642.1      415                 
WP_011274226.1      WP_005293920.1      297                 
WP_011274226.1      WP_011273020.1      1685                
# GCF_000006645.out
# GCF_000006725.out
WP_010893741.1      WP_010892796.1      889                 
WP_010893741.1      WP_010893657.1      90                  
WP_010893741.1      WP_010893799.1      49                  
WP_010893741.1      WP_010894877.1      1129                
WP_010893741.1      WP_031337514.1      886                 
# No columns to parse from file
# No columns to parse from file
# Error - no matching ppm protein found in gff file: ['YP_208574.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['YP_208574.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['YP_208574.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['YP_208574.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['YP_208574.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# GCF_000006845.out
# No columns to parse from file
# No columns to parse from file
# No columns to parse from file
# Error - no matching ppm protein found in gff file: ['NP_384544.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_384544.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_384544.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_384544.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_384544.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_384544.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_384544.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_384544.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_384544.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_384544.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_384544.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_386455.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_386455.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_386455.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_386455.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_386455.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_386455.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_386455.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_386455.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_386455.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_386455.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_386455.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# GCF_000006965.out
# GCF_000007005.out
WP_009989295.1      WP_009988820.1      2132                
WP_009989295.1      WP_009988855.1      2456                
WP_009989295.1      WP_009989252.1      1523                
WP_009989295.1      WP_009989391.1      148                 
WP_009989295.1      WP_009989687.1      610                 
WP_009989295.1      WP_009990173.1      790                 
WP_009989295.1      WP_009990610.1      2193                
WP_009989295.1      WP_009991224.1      1832                
WP_009989295.1      WP_009991321.1      1766                
WP_009989295.1      WP_010923121.1      1665                
WP_009989295.1      WP_010923884.1      101                 
WP_009989295.1      WP_048054196.1      1681                
# GCF_000007125.out
# No columns to parse from file
# GCF_000007165.out
WP_003483895.1      WP_005922502.1      1225                
WP_003483895.1      WP_005929617.1      2595                
WP_003483895.1      WP_011051107.1      625                 
WP_003483895.1      WP_011051296.1      958                 
WP_003483895.1      WP_011052370.1      2590                
WP_003483895.1      WP_015463436.1      1857                
WP_003483895.1      WP_015463628.1      2659                
WP_003488965.1      WP_005922502.1      2155                
WP_003488965.1      WP_005929617.1      3525                
WP_003488965.1      WP_011051107.1      1555                
WP_003488965.1      WP_011051296.1      1888                
WP_003488965.1      WP_011052370.1      3520                
WP_003488965.1      WP_015463436.1      2787                
WP_003488965.1      WP_015463628.1      3589                
# GCF_000007305.out
# [Errno 2] No such file or directory: '/research/gmh/GENOME_DB/gff-Complete/GCF_000007385.gff.gz'
# GCF_000007405.out
# GCF_000007445.out
WP_000052186.1      WP_000011621.1      602                 
WP_000052186.1      WP_000079263.1      2134                
WP_000052186.1      WP_000183079.1      2127                
WP_000052186.1      WP_000246149.1      2868                
WP_000052186.1      WP_000253975.1      3745                
WP_000052186.1      WP_000676064.1      4214                
WP_000052186.1      WP_000718995.1      1235                
WP_000052186.1      WP_000933729.1      4158                
WP_000052186.1      WP_001052151.1      4306                
WP_000052186.1      WP_001272860.1      3004                
WP_000052186.1      WP_001297922.1      2118                
# GCF_000007505.out
# Error - no matching pngc protein found in gff file: ['WP_001533115.1' 'NTP_transferase'], single positional indexer is out-of-bounds
# GCF_000007545.out
WP_000052221.1      WP_000009019.1      1840                
WP_000052221.1      WP_000011583.1      586                 
WP_000052221.1      WP_000253993.1      1601                
WP_000052221.1      WP_000648783.1      1849                
WP_000052221.1      WP_000676077.1      944                 
WP_000052221.1      WP_000729450.1      899                 
WP_000052221.1      WP_000741649.1      348                 
WP_000052221.1      WP_000857536.1      1852                
WP_000052221.1      WP_000934853.1      1250                
WP_000052221.1      WP_000981469.1      1855                
WP_000052221.1      WP_001046874.1      1215                
# Error - no matching ppm protein found in gff file: ['NP_743548.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743548.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743548.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743548.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743548.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743548.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743548.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743548.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743548.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743548.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743548.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743667.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743667.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743667.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743667.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743667.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743667.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743667.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743667.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743667.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743667.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743667.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_744483.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_744483.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_744483.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_744483.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_744483.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_744483.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_744483.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_744483.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_744483.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_744483.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_744483.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_746235.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_746235.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_746235.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_746235.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_746235.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_746235.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_746235.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_746235.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_746235.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_746235.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_746235.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# GCF_000007565.out
# GCF_000007625.out
WP_011099893.1      WP_011098507.1      1609                
WP_011099893.1      WP_011098576.1      1515                
WP_011099893.1      WP_011100598.1      787                 
WP_011099893.1      WP_035111320.1      943                 
WP_011099893.1      WP_035124383.1      574                 
^CTraceback (most recent call last):
  File "protein_distances_no_ppm.py", line 132, in <module>
    main()
  File "protein_distances_no_ppm.py", line 76, in main
    delimiter='\t',
  File "/sw/lib/python2.7/site-packages/pandas/io/parsers.py", line 646, in parser_f
    return _read(filepath_or_buffer, kwds)
  File "/sw/lib/python2.7/site-packages/pandas/io/parsers.py", line 401, in _read
    data = parser.read()
  File "/sw/lib/python2.7/site-packages/pandas/io/parsers.py", line 939, in read
    ret = self._engine.read(nrows)
  File "/sw/lib/python2.7/site-packages/pandas/io/parsers.py", line 1508, in read
    data = self._reader.read(nrows)
  File "pandas/parser.pyx", line 848, in pandas.parser.TextReader.read (pandas/parser.c:9977)
  File "pandas/parser.pyx", line 898, in pandas.parser.TextReader._read_low_memory (pandas/parser.c:10660)
  File "pandas/parser.pyx", line 2046, in pandas.parser._concatenate_chunks (pandas/parser.c:26361)
  File "/sw/lib/python2.7/site-packages/pandas/types/common.py", line 102, in is_categorical_dtype
    def is_categorical_dtype(arr_or_dtype):
KeyboardInterrupt
kuum:virtual_env monica$ python protein_distances_no_ppm.py 
PPM                 PNGC                Distance            
# GCF_000005825.out
# No columns to parse from file
# GCF_000006605.out
WP_011273297.1      WP_005292098.1      1231                
WP_011273297.1      WP_005293089.1      723                 
WP_011273297.1      WP_005293380.1      835                 
WP_011273297.1      WP_005293642.1      868                 
WP_011273297.1      WP_005293920.1      986                 
WP_011273297.1      WP_011273020.1      402                 
WP_011274226.1      WP_005292098.1      52                  
WP_011274226.1      WP_005293089.1      560                 
WP_011274226.1      WP_005293380.1      448                 
WP_011274226.1      WP_005293642.1      415                 
WP_011274226.1      WP_005293920.1      297                 
WP_011274226.1      WP_011273020.1      1685                
# GCF_000006645.out
# GCF_000006725.out
WP_010893741.1      WP_010892796.1      889                 
WP_010893741.1      WP_010893657.1      90                  
WP_010893741.1      WP_010893799.1      49                  
WP_010893741.1      WP_010894877.1      1129                
WP_010893741.1      WP_031337514.1      886                 
# No columns to parse from file
# No columns to parse from file
# Error - no matching ppm protein found in gff file: ['YP_208574.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['YP_208574.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['YP_208574.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['YP_208574.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['YP_208574.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# GCF_000006845.out
# No columns to parse from file
# No columns to parse from file
# No columns to parse from file
# Error - no matching ppm protein found in gff file: ['NP_384544.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_384544.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_384544.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_384544.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_384544.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_384544.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_384544.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_384544.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_384544.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_384544.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_384544.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_386455.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_386455.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_386455.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_386455.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_386455.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_386455.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_386455.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_386455.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_386455.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_386455.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_386455.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# GCF_000006965.out
# GCF_000007005.out
WP_009989295.1      WP_009988820.1      2132                
WP_009989295.1      WP_009988855.1      2456                
WP_009989295.1      WP_009989252.1      1523                
WP_009989295.1      WP_009989391.1      148                 
WP_009989295.1      WP_009989687.1      610                 
WP_009989295.1      WP_009990173.1      790                 
WP_009989295.1      WP_009990610.1      2193                
WP_009989295.1      WP_009991224.1      1832                
WP_009989295.1      WP_009991321.1      1766                
WP_009989295.1      WP_010923121.1      1665                
WP_009989295.1      WP_010923884.1      101                 
WP_009989295.1      WP_048054196.1      1681                
# GCF_000007125.out
# No columns to parse from file
# GCF_000007165.out
WP_003483895.1      WP_005922502.1      1225                
WP_003483895.1      WP_005929617.1      2595                
WP_003483895.1      WP_011051107.1      625                 
WP_003483895.1      WP_011051296.1      958                 
WP_003483895.1      WP_011052370.1      2590                
WP_003483895.1      WP_015463436.1      1857                
WP_003483895.1      WP_015463628.1      2659                
WP_003488965.1      WP_005922502.1      2155                
WP_003488965.1      WP_005929617.1      3525                
WP_003488965.1      WP_011051107.1      1555                
WP_003488965.1      WP_011051296.1      1888                
WP_003488965.1      WP_011052370.1      3520                
WP_003488965.1      WP_015463436.1      2787                
WP_003488965.1      WP_015463628.1      3589                
# GCF_000007305.out
# [Errno 2] No such file or directory: '/research/gmh/GENOME_DB/gff-Complete/GCF_000007385.gff.gz'
# GCF_000007405.out
# GCF_000007445.out
WP_000052186.1      WP_000011621.1      602                 
WP_000052186.1      WP_000079263.1      2134                
WP_000052186.1      WP_000183079.1      2127                
WP_000052186.1      WP_000246149.1      2868                
WP_000052186.1      WP_000253975.1      3745                
WP_000052186.1      WP_000676064.1      4214                
WP_000052186.1      WP_000718995.1      1235                
WP_000052186.1      WP_000933729.1      4158                
WP_000052186.1      WP_001052151.1      4306                
WP_000052186.1      WP_001272860.1      3004                
WP_000052186.1      WP_001297922.1      2118                
# GCF_000007505.out
# Error - no matching pngc protein found in gff file: ['WP_001533115.1' 'NTP_transferase'], single positional indexer is out-of-bounds
# GCF_000007545.out
WP_000052221.1      WP_000009019.1      1840                
WP_000052221.1      WP_000011583.1      586                 
WP_000052221.1      WP_000253993.1      1601                
WP_000052221.1      WP_000648783.1      1849                
WP_000052221.1      WP_000676077.1      944                 
WP_000052221.1      WP_000729450.1      899                 
WP_000052221.1      WP_000741649.1      348                 
WP_000052221.1      WP_000857536.1      1852                
WP_000052221.1      WP_000934853.1      1250                
WP_000052221.1      WP_000981469.1      1855                
WP_000052221.1      WP_001046874.1      1215                
# Error - no matching ppm protein found in gff file: ['NP_743548.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743548.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743548.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743548.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743548.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743548.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743548.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743548.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743548.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743548.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743548.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743667.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743667.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743667.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743667.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743667.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743667.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743667.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743667.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743667.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743667.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_743667.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_744483.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_744483.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_744483.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_744483.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_744483.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_744483.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_744483.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_744483.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_744483.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_744483.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_744483.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_746235.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_746235.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_746235.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_746235.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_746235.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_746235.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_746235.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_746235.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_746235.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_746235.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_746235.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# GCF_000007565.out
# GCF_000007625.out
WP_011099893.1      WP_011098507.1      1609                
WP_011099893.1      WP_011098576.1      1515                
WP_011099893.1      WP_011100598.1      787                 
WP_011099893.1      WP_035111320.1      943                 
WP_011099893.1      WP_035124383.1      574                 
# Error - no matching ppm protein found in gff file: ['WP_001984053.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['WP_001984053.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['WP_001984053.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['WP_001984053.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['WP_001984053.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['WP_001984053.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['WP_001984053.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# GCF_000007685.out
# Error - no matching ppm protein found in gff file: ['WP_011133945.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['WP_011133945.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['WP_011133945.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['WP_011133945.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['WP_011133945.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['WP_011133945.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['WP_011133945.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['WP_011133945.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['WP_011133945.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['WP_011133945.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['WP_011133945.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching pngc protein found in gff file: ['WP_043596804.1' 'NTP_transferase'], single positional indexer is out-of-bounds
# Error - no matching pngc protein found in gff file: ['WP_043596804.1' 'NTP_transferase'], single positional indexer is out-of-bounds
# Error - no matching pngc protein found in gff file: ['WP_043596804.1' 'NTP_transferase'], single positional indexer is out-of-bounds
# GCF_000007705.out
WP_011135193.1      WP_011134229.1      972                 
WP_011135193.1      WP_011135177.1      16                  
WP_011135193.1      WP_011136891.1      1699                
WP_011135193.1      WP_011137440.1      2253                
WP_011135193.1      WP_011137448.1      2262                
WP_011135193.1      WP_011137577.1      2392                
WP_011135193.1      WP_043595631.1      395                 
WP_011135193.1      WP_043595850.1      34                  
WP_011135193.1      WP_043596477.1      1522                
WP_011135193.1      WP_043598400.1      2598                
WP_011135605.1      WP_011134229.1      1371                
WP_011135605.1      WP_011135177.1      415                 
WP_011135605.1      WP_011136891.1      1300                
WP_011135605.1      WP_011137440.1      1854                
WP_011135605.1      WP_011137448.1      1863                
WP_011135605.1      WP_011137577.1      1993                
WP_011135605.1      WP_043595631.1      794                 
WP_011135605.1      WP_043595850.1      365                 
WP_011135605.1      WP_043596477.1      1123                
WP_011135605.1      WP_043598400.1      2199                
WP_043595771.1      WP_011134229.1      840                 
WP_043595771.1      WP_011135177.1      116                 
WP_043595771.1      WP_011136891.1      1831                
WP_043595771.1      WP_011137440.1      2385                
WP_043595771.1      WP_011137448.1      2394                
WP_043595771.1      WP_011137577.1      2524                
WP_043595771.1      WP_043595631.1      263                 
WP_043595771.1      WP_043595850.1      166                 
WP_043595771.1      WP_043596477.1      1654                
WP_043595771.1      WP_043598400.1      2730                
# No columns to parse from file
# No columns to parse from file
# No columns to parse from file
# Error - no matching ppm protein found in gff file: ['NP_843617.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_843617.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_843617.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_843617.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_843617.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_843617.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_843617.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_844733.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_844733.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_844733.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_844733.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_844733.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_844733.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# Error - no matching ppm protein found in gff file: ['NP_844733.1' 'PEP_mutase'], single positional indexer is out-of-bounds
# GCF_000007845.out
# Error - no matching pngc protein found in gff file: ['WP_042632844.1' 'NTP_transf_3'], single positional indexer is out-of-bounds
# Error - no matching pngc protein found in gff file: ['WP_042632844.1' 'NTP_transf_3'], single positional indexer is out-of-bounds
# GCF_000007865.out
WP_003875874.1      WP_003872738.1      1313                
WP_003875874.1      WP_003874547.1      1199                
WP_003875874.1      WP_003877502.1      1375                
WP_003875874.1      WP_003878468.1      325                 
WP_003875874.1      WP_003879222.1      1669                
WP_003875874.1      WP_019683876.1      1839                
WP_003876996.1      WP_003872738.1      722                 
WP_003876996.1      WP_003874547.1      3234                
WP_003876996.1      WP_003877502.1      660                 
WP_003876996.1      WP_003878468.1      2360                
WP_003876996.1      WP_003879222.1      3704                
WP_003876996.1      WP_019683876.1      196                 
# GCF_000007885.out
# Error - no matching pngc protein found in gff file: ['WP_041309403.1' 'IspD'], single positional indexer is out-of-bounds
# GCF_000007905.out
WP_011115239.1      WP_011114920.1      318                 
WP_011115239.1      WP_011114943.1      296                 
WP_011115239.1      WP_011115517.1      269                 
WP_011115239.1      WP_011115533.1      285                 
WP_011115239.1      WP_011115859.1      606                 
WP_011115239.1      WP_011116279.1      998                 
WP_011115239.1      WP_034366326.1      173                 
# GCF_000008005.out
WP_000312648.1      WP_000049469.1      2468                
WP_000312648.1      WP_000057613.1      2579                
WP_000312648.1      WP_000071036.1      2333                
WP_000312648.1      WP_000288294.1      2292                
WP_000312648.1      WP_000648792.1      972                 
WP_000312648.1      WP_000676173.1      1005                
WP_000312648.1      WP_000679295.1      1915                
WP_000312648.1      WP_000757828.1      2972                
WP_000312648.1      WP_011040800.1      2636                
WP_000787356.1      WP_000049469.1      3579                
WP_000787356.1      WP_000057613.1      3690                
WP_000787356.1      WP_000071036.1      1222                
WP_000787356.1      WP_000288294.1      1181                
WP_000787356.1      WP_000648792.1      2083                
WP_000787356.1      WP_000676173.1      106                 
WP_000787356.1      WP_000679295.1      3026                
WP_000787356.1      WP_000757828.1      4083                
WP_000787356.1      WP_011040800.1      3747                
# Error - no matching pngc protein found in gff file: ['WP_001533115.1' 'NTP_transferase'], single positional indexer is out-of-bounds
# GCF_000008105.out
WP_000052229.1      WP_000011576.1      598                 
WP_000052229.1      WP_000253994.1      3342                
WP_000052229.1      WP_000632420.1      1852                
WP_000052229.1      WP_000729450.1      1480                
WP_000052229.1      WP_000741645.1      2693                
WP_000052229.1      WP_000981469.1      1857                
WP_000052229.1      WP_001048454.1      3809                
WP_000052229.1      WP_001541189.1      3737                
WP_000052229.1      WP_011264426.1      3673                
^CTraceback (most recent call last):
  File "protein_distances_no_ppm.py", line 132, in <module>
    main()
  File "protein_distances_no_ppm.py", line 103, in main
    ppm_gff = gff_df[gff_df.protein_id == pair[0][0]].iloc[0]
  File "/sw/lib/python2.7/site-packages/pandas/core/ops.py", line 863, in wrapper
    res = pd.Series(res, index=self.index, name=self.name, dtype='bool')
  File "/sw/lib/python2.7/site-packages/pandas/core/series.py", line 243, in __init__
    raise_cast_failure=True)
  File "/sw/lib/python2.7/site-packages/pandas/core/series.py", line 2867, in _sanitize_array
    subarr = _try_cast(data, True)
  File "/sw/lib/python2.7/site-packages/pandas/core/series.py", line 2842, in _try_cast
    subarr = _possibly_cast_to_datetime(arr, dtype)
  File "/sw/lib/python2.7/site-packages/pandas/types/cast.py", line 769, in _possibly_cast_to_datetime
    from pandas.tseries.timedeltas import to_timedelta
KeyboardInterrupt
kuum:virtual_env monica$ nano protein_distances_no_ppm.py 
kuum:virtual_env monica$ python protein_distances_no_lica.py
LICA                PNGC                Distance            
# GCF_000005825.out
Traceback (most recent call last):
  File "protein_distances_no_lica.py", line 132, in <module>
    main()
  File "protein_distances_no_lica.py", line 128, in main
    [DATA_FMT.format(x) for x in [lica[0], pngc[0], distance]]
NameError: global name 'lica' is not defined
kuum:virtual_env monica$ nano protein_distances_no_ppm.py 
kuum:virtual_env monica$ nano protein_distances_no_lica.py
kuum:virtual_env monica$ python protein_distances_no_lica.py
LICA                PNGC                Distance            
# GCF_000005825.out
Traceback (most recent call last):
  File "protein_distances_no_lica.py", line 132, in <module>
    main()
  File "protein_distances_no_lica.py", line 128, in main
    [DATA_FMT.format(x) for x in [lica[0], pngc[0], distance]]
NameError: global name 'lica' is not defined
kuum:virtual_env monica$ mv protein_distances_no_lica.py protein_distances_lica.py
kuum:virtual_env monica$ nano protein_distances_lica.py 
kuum:virtual_env monica$ nano util.py
kuum:virtual_env monica$ nano util.py

  GNU nano 2.0.6                                   File: util.py                                                                             

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

import os, re
import gzip

from util import parse_faa

database = "Complete"

FAA_PATH = '/research/gmh/GENOME_DB/faa-' + database
INPUT_PATH = '/Vagabundo/monica/temp/70-CUTGA-OUT-faa-' + database

MOTIF_REGEX = 'EDK\w{3,7}NS'

# Check filenames in 70-CUTGA (.out) and scan for those filenames in faa-Complete(.faa)
# Read filenames in faa-Complete (amino acid sequences) to scan for conserved motif (EDKXXXXXNS)
# The XXXXX is a variable region where X can be any letter
# The variable region itself can vary in number of X - scan for 3 - 7X

# Output formatted like:
# GCF_A:WP_A,WP_B,WP_C,WP_D
# GCF_B:WP_E,WP_F,WP_G,WP_H


def has_motif(content_string, motif_regex):
    return bool(re.search(motif_regex, content_string))


def file_has_motif(faa_path, motif_regex):
    try:
        if faa_path.endswith('.gz'):
            with gzip.open(faa_path, 'r') as f:
                file_contents = bytes(str(f.read()), 'utf-8')
        else:
            with open(faa_path, 'r') as f:
                file_contents = f.read()
        
        has_motif = bool(re.search(motif_regex, file_contents))
        return has_motif
    except Exception as e:
        print('# ERROR {}: {}'.format(e, faa_path))
        return False


if __name__ == '__main__':
    out_files = [f for f in os.listdir(INPUT_PATH) if f.endswith('.out')]
    for fname in out_files:
        gcf_id = fname.split('.')[0]

        # Use .faa.gz as extension, to read gzipped files with parse_faa
        faa_fname = fname.split('.')[0] + '.faa.gz'

        faa_path = FAA_PATH + '/' + faa_fname
        if file_has_motif(faa_path, MOTIF_REGEX):

            faa_sequences = parse_faa(faa_path)
            matching_wps = []
            for wp_id, fasta_data in faa_sequences.iteritems():
                if has_motif(fasta_data.sequence, MOTIF_REGEX):
                    matching_wps.append(wp_id)
            print(gcf_id + ':' + ','.join(matching_wps))

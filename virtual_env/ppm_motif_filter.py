import os, re

database = "Complete"

FAA_PATH = '/Vagabundo/monica/temp/faa-' + database
INPUT_PATH = '/Vagabundo/monica/temp/70-CUTGA-OUT-faa-' + database

MOTIF_REGEX = 'EDK\w{3,7}NS'

# Check filenames in 70-CUTGA (.out) and scan for those filenames in faa-Complete(.faa)
# Read filenames in faa-Complete (amino acid sequences) to scan for conserved motif (EDKXXXXXNS)
# The XXXXX is a variable region where X can be any letter
# The variable region itself can vary in number of X - scan for 3 - 7X

# Output formatted like:
# GCF_A
# WP_A,WP_B,WP_C,WP_D
# GCF_B
# WP_A,WP_B,WP_C,WP_D


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
        wp_id = seq_lines[0].split()[0]
        sequence = ''.join(seq_lines[1:])
        faa_data[wp_id] = sequence

    return faa_data


def has_motif(content_string, motif_regex):
    return bool(re.search(motif_regex, content_string))


def file_has_motif(faa_path, motif_regex):
    with open(faa_path, 'r') as f:
        file_contents = f.read()
    has_motif = bool(re.search(motif_regex, file_contents))
    return has_motif


if __name__ == '__main__':
    out_files = [f for f in os.listdir(INPUT_PATH) if f.endswith('.out')]
    for fname in out_files:
        gcf_id = fname.split('.')[0]
        faa_fname = fname.split('.')[0] + '.faa'

        faa_path = FAA_PATH + '/' + faa_fname
        if file_has_motif(faa_path, MOTIF_REGEX):

            faa_sequences = parse_faa(faa_path)
            matching_wps = []
            for wp_id, seq in faa_sequences.iteritems():
                if has_motif(seq, MOTIF_REGEX):
                    matching_wps.append(wp_id)
            print gcf_id + ':' + ','.join(matching_wps)

#!/sw/bin/python
import sys
import os
import subprocess

PSIBLAST = '/usr/local/ncbi/bin/psiblast'
NUM_ITERATIONS = 1
INPUT_DIR = '/Vagabundo/monica/psiblast_inputs/'
OUTPUT_DIR = '/Vagabundo/monica/psiblast_outputs/'

FAIL_LIST = 'failed_inputs.txt'

# NOTE: run this file like so: python -m psiblast_runner
# TODO: make the script check if there is already an output file named after it
# NOTE/TODO: shoot, this names the output files literally "filename.fasta" in the output dir but they're just txt files

def run_psiblast(INPUT_DIR, OUTPUT_DIR):
    failures = []
    inputs = [f for f in os.listdir(INPUT_DIR) if os.path.isfile(INPUT_DIR+f) and f.lower().endswith('.fasta')]
    print "Found {} fasta files in {}. Each one may take ... 15-20 minutes... Beware.".format(len(inputs), INPUT_DIR)

    for fname in inputs:
        try:
            output = subprocess.check_output('{} -query {} -num_iterations {}'.format(PSIBLAST, INPUT_DIR+fname, NUM_ITERATIONS))
            with open(OUTPUT_DIR+fname, 'w') as opath:
                opath.write(output)
        except:
            print "Psiblast failed for input file {}".format(fpath)
            failures.append(fpath)

    print "Finished doing everything I could!"
    if len(failures) > 0:
        with open(INPUT_DIR+FAIL_LIST, 'w') as failed_path:
            failed_path.write('\n'.join(failures))
        print "Check {}{} for a list of the inputs that failed being run in psiblast...".format(INPUT_DIR,FAIL_LIST)

def main():
    run_psiblast(INPUT_DIR, OUTPUT_DIR)

if __name__ == '__main__':
    main()

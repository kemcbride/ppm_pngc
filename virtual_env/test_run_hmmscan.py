from __future__ import print_function

from util import run_hmmscan, PFAM_DB

from argparse import ArgumentParser


if __name__ == '__main__':
    ap = ArgumentParser()
    ap.add_argument('input')
    ap.add_argument('output')
    args = ap.parse_args()


    print("I'm going to run hmmscan on {} and {} now!")

    run_hmmscan(args.input, args.output)

    print ( "I'm done running hmmscan.")

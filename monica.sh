#!/bin/bash
# Usage: ./monica.sh
set -e
mkdir -p ~/temp/faa-Scaffold
mkdir -p ~/temp/OUT-faa-Scaffold

total_files=26909 # Scaffold is too big, ls-1 | wc-l fails
curr_file=0

pushd /research/gmh/GENOME_DB/faa-Scaffold/ > /dev/null
for f in *.gz
do
	pushd ~/temp/faa-Scaffold > /dev/null
	gunzip -c /research/gmh/GENOME_DB/faa-Scaffold/$f > ~/temp/faa-Scaffold/${f%.gz}
	echo -e -n "Processing $f file. $curr_file/$total_files\\r"
	hmmscan --domtblout ~/temp/OUT-faa-Scaffold/${f%.*}.out ~/Proteins/models.hmm ${f%.gz} > /dev/null
	curr_file=$(($curr_file+1))
	popd > /dev/null
done
popd > /dev/null

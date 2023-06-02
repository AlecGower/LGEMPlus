#!/bin/bash
printf "Creating theory files for ... yeastGEM" &&
	python3 experiments/create_theory_files.py "yeastGEM" &&
	print "Done!"
printf "Creating theory files for ... iMM904" &&
	python3 experiments/create_theory_files.py "iMM904" &&
	print "Done!"
printf "Creating theory files for ... iND750" &&
	python3 experiments/create_theory_files.py "iND750" &&
	print "Done!"

for model in $(ls experiments/theories/); do
	python3 experiments/single_gene_deletions.py \
		experiments/theories/$model/$(ls experiments/theories/$model | tail -n 1)
done

#!/bin/bash
printf "Creating theory files for yeastGEM... " &&
	python3 helpers/create_theory_files.py "yeastGEM" &&
	printf "Done!\n"
printf "Creating theory files for iMM904... " &&
	python3 helpers/create_theory_files.py "iMM904" &&
	printf "Done!\n"
printf "Creating theory files for iND750... " &&
	python3 helpers/create_theory_files.py "iND750" &&
	printf "Done!\n"

for model in $(ls experiments/theories/); do
	latest=$model/$(ls experiments/theories/$model | tail -n 1)
	mkdir -p experiments/results/$latest
	python3 helpers/single_gene_deletions.py \
		experiments/theories/$latest \
		>experiments/results/$latest/single_gene_deletions.txt
	python3 helpers/lethality_classification.py \
		experiments/results/$latest/single_gene_deletions.txt \
		>experiments/results/$latest/single_gene_deletions_evaluation.txt
done

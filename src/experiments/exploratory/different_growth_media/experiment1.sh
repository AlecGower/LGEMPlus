#!/bin/bash
printf "Creating theory files for yeastGEM... " &&
	python3 helpers/create_theory_files.py "yeastGEM" $1 &&
	printf "Done!\n"
# printf "Creating theory files for iMM904... " &&
# 	python3 helpers/create_theory_files.py "iMM904" $1 &&
# 	printf "Done!\n"
# printf "Creating theory files for iFF708... " &&
# 	python3 helpers/create_theory_files.py "iFF708" $1 &&
# 	printf "Done!\n"

for model in $(ls -t experiments/theories/ | head -n 1); do
	latest=$model/$(ls experiments/theories/$model | tail -n 1)
	mkdir -p experiments/results/$latest
	python3 helpers/single_gene_deletions.py \
		experiments/theories/$latest \
		>experiments/results/$latest/single_gene_deletions.txt
	python3 helpers/lethality_classification.py \
		experiments/results/$latest/single_gene_deletions.txt \
		>experiments/results/$latest/single_gene_deletions_evaluation.txt
	python3 helpers/single_gene_deletions_pathways.py \
		experiments/theories/$latest \
		>experiments/results/$latest/single_gene_deletions_pathways.txt
	(cd experiments/theories/$latest &&
		$IPROVER_HOME"/iproveropt" reactions.p compound_synonyms.p genes.p media_compounds.p \
			ubiquitous_compounds.p abduced_extra_compounds.p query.p |
		grep "__in_genome" |
			sed -E "s/^.*g_([A-Z0-9_]+)_([A-Z0-9_]+)__in_genome.*$/\1 \2/") \
		>experiments/results/$latest/wt_pathway.txt
	python3 helpers/compare_pathways.py \
		experiments/results/$latest \
		>experiments/results/$latest/wt_deviation_pathways.txt
done

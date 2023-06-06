#!/bin/bash

###
# Experiment 2: Abduction of hypotheses
# For each NGG error in the list:
#   Run iProver and extract hypotheses to fix the error
###

for model in $(ls experiments/theories/); do
	if [[ $model =~ "yeastGEM" ]]; then
		latest=$model/$(ls experiments/theories/$model | tail -n 1)
		cut -d " " -f 1 experiments/results/$latest/log_ngg.txt |
			xargs -n 1 -P 8 -I{} /bin/bash \
				helpers/abduction/ngg_abduction_preparation.sh \
				{} \
				experiments/theories/$latest/genes.p \
				experiments/theories/$latest/reactions.p \
				experiments/theories/$latest/media_compounds.p \
				experiments/theories/$latest/ubiquitous_compounds.p \
				experiments/theories/$latest/abduced_extra_compounds.p \
				experiments/theories/$latest/query.p \
		1>experiments/results/$latest/"abduction.out" \
			2>experiments/results/$latest/"abduction.log"
	fi
done

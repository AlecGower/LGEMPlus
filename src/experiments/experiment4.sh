#!/bin/bash

###
# Experiment 4: Investigation of GNG errors
# For each GNG error in the list:
#   Run iProver and extract reactions in proof
#   Run FBA constraint using reactions
###

for model in $(ls experiments/theories/); do
	if [[ $model =~ "yeastGEM" ]]; then
		latest=$model/$(ls experiments/theories/$model | tail -n 1)
		printf "BASE_STRAIN\n$(cut -d " " -f 1 experiments/results/$latest/log_gng.txt)"|
			xargs -n 1 -P 8 -I{} timeout 300 /bin/bash \
				helpers/gng_investigation.sh \
				{} \
				experiments/theories/$latest/genes.p \
				experiments/theories/$latest/reactions.p \
				experiments/theories/$latest/media_compounds.p \
				experiments/theories/$latest/ubiquitous_compounds.p \
				experiments/theories/$latest/abduced_extra_compounds.p \
				experiments/theories/$latest/query.p \
		1>experiments/results/$latest/"gng_investigation.out" \
			2>experiments/results/$latest/"gng_investigation.log"
	fi
done

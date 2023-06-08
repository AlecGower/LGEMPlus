for model in $(ls experiments/theories/); do
	if [[ $model =~ "yeastGEM" ]]; then
		latest=$model/$(ls experiments/theories/$model | tail -n 1)

		for hypothesis_problem_file in $(ls experiments/results/$latest/ngg_hypotheses/*.p); do
			while read orf hypothesis; do
				echo $orf $hypothesis $(python3 helpers/single_gene_deletions.py -c \
					-g experiments/results/$latest/log_ngng.txt \
					-p $hypothesis_problem_file \
					experiments/theories/$latest) 
			done <<<$(echo $hypothesis_problem_file | xargs -I{} basename {} | sed -E "s/^([^_]+)_hypothesis_([0-9]+)\.p$/\1 \2/")
		done

	fi
done

#!/bin/bash
standardise_hypothesis() {
	echo $1 | sed "s/~(/~/g ; s/);/;/g ; s/)}/}/g ; s/^{// ; s/}$// ; " | tr ";" "\n" | sort | tr "\n" ";" | sed "s/;$/}\n/; s/^/{/"
}
export -f standardise_hypothesis

for model in $(ls experiments/theories/); do
	if [[ $model =~ "yeastGEM" ]]; then
		latest=$model/$(ls experiments/theories/$model | tail -n 1)
		for proof_fp in $(ls experiments/results/$latest/ngg_hypotheses/*_proof.flattened); do
			orf=$(echo $proof_fp | xargs -I{} basename {} | sed -E "s/^([^_]+)_proof.flattened$/\1/")
            hyp_count=0
            for hyp in $(grep "{~(" $proof_fp|grep -o "{.*}$"); do
                ((hyp_count++))
                hyp_sorted=$(echo $hyp | xargs -I{} /bin/bash -c 'standardise_hypothesis "$@"' _ {})
                echo $orf $hyp_count $hyp_sorted
            done
		done
	fi
done

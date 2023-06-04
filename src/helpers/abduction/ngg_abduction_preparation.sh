#!/bin/bash
alias iprover=$IPROVER_HOME"/iproveropt"
alias iprover_simple="iprover --proof_out false --sat_out_model none --sat_out_clauses false --sat_out_model none"
alias iprover_dot="iprover --schedule none --preprocessing_flag false --instantiation_flag false --superposition_flag false --resolution_flag true  --res_prop_simpl_given false --res_to_prop_solver none --prop_solver_per_cl 0  --proof_reduce_dot \"[all_neg;unit]\" --proof_dot_file"

###
# The purpose of this script is to take an ORF as input as well as the problem files relevant and return: hypotheses for each gene in the ORF that ``fix'' the deletant; for each of these hypotheses a list of reactions that are activated in the proof.
###

echo >&2 "Starting abduction for "$1"."

output_directory=$(dirname $2|sed "s/theories/results/")/ngg_hypotheses
mkdir -p $output_directory

##
# Run knockout
out_fp=""$output_directory"/"$1"_proof.flattened"
helpers/gene_knockout.sh $1 ${@:2} | helpers/flatten_iprover_output.sh >$out_fp

## Extract hypotheses and run tests
hypCount=1 # Set the count to 1
# Create empty files for the list of hypotheses and for the first problem file
cat /dev/null >""$output_directory"/"$1"_hypothesis_"$hypCount".p"
cat /dev/null >""$output_directory"/"$1".hypotheses"
while read -r hyp add_clause_component; do
	echo $hyp $add_clause_component >>""$output_directory"/"$1".hypotheses"
	# If we have started a new
	if [[ $hypCount > $hyp ]]; then
		continue
	elif [[ $hypCount == $hyp ]]; then
		if [[ ! $add_clause_component =~ "met\(" ]]; then
			echo >&2 "Non-metabolite component found in clause ("$1"): \""$add_clause_component"\""
			echo >&2 "Skipping hypothesis "$hypCount" ("$1")..."
            echo >&1 $1 $hypCount $add_clause_component "NonMet"
			((hypCount++))
			continue
		fi
		echo $add_clause_component >>""$output_directory"/"$1"_hypothesis_"$hypCount".p"
	elif [[ $hypCount < $hyp ]]; then # we have reached the first component of the next hypothesis
		echo >&2 ""$1" - Proving hypothesis "$hypCount" ("$1")..."
		# iprover_simple ${@:2} ""$output_directory"/"$1"_hypothesis_"$hypCount".p" |grep "% SZS status"
		iprover_dot ""$output_directory"/"$1"_hypothesis_"$hypCount".dot" ${@:2} ""$output_directory"/"$1"_hypothesis_"$hypCount".p" >/dev/null
		grep -v " -> " ""$output_directory"/"$1"_hypothesis_"$hypCount".dot" | grep rxn | perl -pe 's|^.*(r_\d{4})\w+?(reverse)?.*$|\1 \2|' | sort >""$output_directory"/"$1"_hypothesis_"$hypCount".rxns"

		# Solve FBA with the gene knockout
		python helpers/abduction/LMtoFBA_proof.py --reactions ""$output_directory"/"$1"_hypothesis_"$hypCount".rxns" --ko_gene $1 --model "ModelFiles/yeastGEM.xml" --mode "promote" --threshold=1e-9 --verbose=True
	    echo >&1 $1 $hypCount $add_clause_component ""$output_directory"/"$1"_hypothesis_"$hypCount".rxns"

		# Increment
		((hypCount++))
		cat /dev/null >""$output_directory"/"$1"_hypothesis_"$hypCount".p"
		echo $add_clause_component >>""$output_directory"/"$1"_hypothesis_"$hypCount".p"
	fi
done < <(helpers/abduction/get_hypotheses_only_from_sat.pl $out_fp | sort -in)
if [[ $hypCount == $hyp ]]; then
	# iProver for last g hypo hypothesis
	echo >&2 ""$1" - Proving hypothesis "$hypCount" ("$1")..."
	# iprover_simple ${@:2} ""$output_directory"/"$1"_hypothesis_"$hypCount".p" |grep "% SZS status"
	iprover_dot ""$output_directory"/"$1"_hypothesis_"$hypCount".dot" ${@:2} ""$output_directory"/"$1"_hypothesis_"$hypCount".p" >/dev/null
	grep -v " -> " ""$output_directory"/"$1"_hypothesis_"$hypCount".dot" | grep rxn | perl -pe 's|^.*(r_\d{4})\w+?(reverse)?.*$|\1 \2|' | sort >""$output_directory"/"$1"_hypothesis_"$hypCount".rxns"

    # Solve FBA with the gene knockout
    python helpers/abduction/LMtoFBA_proof.py --reactions ""$output_directory"/"$1"_hypothesis_"$hypCount".rxns" --ko_gene $1 --model "ModelFiles/yeastGEM.xml" --mode "promote" --threshold=1e-9 --verbose=True
    echo >&1 $1 $hypCount $add_clause_component ""$output_directory"/"$1"_hypothesis_"$hypCount".rxns"
fi

echo >&2 "Finishing abduction for "$1"."

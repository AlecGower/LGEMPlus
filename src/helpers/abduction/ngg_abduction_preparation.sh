#!/bin/bash

###
# The purpose of this script is to take an ORF as input as well as the problem files relevant and return: hypotheses for each gene in the ORF that ``fix'' the deletant; for each of these hypotheses a list of reactions that are activated in the proof.
###

echo >&2 "Starting abduction for "$1"."

output_directory=$(dirname $2 | sed "s/theories/results/")/ngg_hypotheses
info_file=$(dirname $2)/info.txt
base_gem_filepath=$(grep "Base GEM filepath" $info_file | cut -d ":" -f 2 | awk '{$1=$1};1')
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
while read -r hyp clause_count add_clause_component; do
	echo $hyp $clause_count $add_clause_component >>""$output_directory"/"$1".hypotheses"
	# If we have started a new
	if [[ $hypCount > $hyp ]]; then
		continue
	elif [[ $hypCount == $hyp ]]; then
		if [[ ! $add_clause_component =~ met\( ]]; then
			echo >&2 "Non-metabolite component found in clause ("$1"): \""$add_clause_component"\""
			echo >&2 "Skipping hypothesis "$hypCount" ("$1")..."
			echo >&1 $1 $hypCount $add_clause_component "NonMet"
			((hypCount++))
			continue
		fi
		echo $add_clause_component >>""$output_directory"/"$1"_hypothesis_"$hypCount".p"
	elif [[ $hypCount < $hyp ]]; then # we have reached the first component of the next hypothesis
		echo >&2 ""$1" - Proving hypothesis "$hypCount" ("$1")..."
		echo >&2 "Checking that the hypothesis does fix it: " $(helpers/gene_knockout_simple.sh $1 $2 ${@:3} \
			""$output_directory"/"$1"_hypothesis_"$hypCount".p" | grep "% SZS status")
		dot_path=""$output_directory"/"$1"_hypothesis_"$hypCount".dot"
		helpers/gene_knockout_dot.sh $dot_path $1 $2 ${@:3} \
			""$output_directory"/"$1"_hypothesis_"$hypCount".p" >/dev/null
		grep -v " -> " $dot_path |
			grep rxn | perl -pe 's|^.*(r_\d{4})\w+?(reverse)?.*$|\1 \2|' | sort \
			>""$output_directory"/"$1"_hypothesis_"$hypCount".rxns"

		# Solve FBA with the gene knockout
		python3 helpers/abduction/LMtoFBA_proof.py \
			--reactions ""$output_directory"/"$1"_hypothesis_"$hypCount".rxns" \
			--ko_gene $1 --model $base_gem_filepath \
			--mode "promote" --threshold=1e-9 --verbose=True
			# $(grep "Base deletants" $info_file | cut -d ":" -f 2 | awk '{$1=$1};1' | sed "s/ / --ko_gene /g ; s/^/--ko_gene /") \
			# --ignore_rxn "r_3547" --ignore_rxn "r_4395" --ignore_rxn "r_0473" \
			# --ignore_rxn "r_0982" --ignore_rxn "r_0984" --ignore_rxn "r_2096" \
			# --ignore_rxn "r_3532" --ignore_rxn "r_3546"
		echo >&1 $1 $hypCount $add_clause_component ""$output_directory"/"$1"_hypothesis_"$hypCount".rxns"

		# Increment
		((hypCount++))
		cat /dev/null >""$output_directory"/"$1"_hypothesis_"$hypCount".p"
		echo $add_clause_component >>""$output_directory"/"$1"_hypothesis_"$hypCount".p"
	fi
done < <(helpers/abduction/get_hypotheses_only_from_sat.pl $out_fp | sort -k1,1 -k2,2 -n)
if [[ $hypCount == $hyp ]]; then
	# iProver for last g hypo hypothesis
	echo >&2 ""$1" - Proving hypothesis "$hypCount" ("$1")..."
	echo >&2 "Checking that the hypothesis does fix it: " $(helpers/gene_knockout_simple $1 $2 ${@:3} \
		""$output_directory"/"$1"_hypothesis_"$hypCount".p" | grep "% SZS status")
	dot_path=""$output_directory"/"$1"_hypothesis_"$hypCount".dot"
	helpers/gene_knockout_dot.sh $dot_path $1 $2 ${@:3} \
		""$output_directory"/"$1"_hypothesis_"$hypCount".p" >/dev/null
	grep -v " -> " $dot_path |
		grep rxn | perl -pe 's|^.*(r_\d{4})\w+?(reverse)?.*$|\1 \2|' | sort \
		>""$output_directory"/"$1"_hypothesis_"$hypCount".rxns"

	# Solve FBA with the gene knockout
	python3 helpers/abduction/LMtoFBA_proof.py \
		--reactions ""$output_directory"/"$1"_hypothesis_"$hypCount".rxns" \
		$(grep "Base deletants" $info_file | cut -d ":" -f 2 | awk '{$1=$1};1' | sed "s/ / --ko_gene /g ; s/^/--ko_gene /") \
		--ko_gene $1 --model $base_gem_filepath \
		--mode "promote" --threshold=1e-9 --verbose=True
		# --ignore_rxn "r_3547" --ignore_rxn "r_4395" --ignore_rxn "r_0473" \
		# --ignore_rxn "r_0982" --ignore_rxn "r_0984" --ignore_rxn "r_2096" \
		# --ignore_rxn "r_3532" --ignore_rxn "r_3546"
	echo >&1 $1 $hypCount $add_clause_component ""$output_directory"/"$1"_hypothesis_"$hypCount".rxns"
fi

echo >&2 "Finishing abduction for "$1"."

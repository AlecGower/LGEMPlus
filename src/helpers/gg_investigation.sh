#!/bin/bash

###
# The purpose of this script is to take an ORF as input as well as the problem files relevant and return: hypotheses for each gene in the ORF that ``fix'' the deletant; for each of these hypotheses a list of reactions that are activated in the proof.
###

echo >&2 "Starting investigation for "$1"."

output_directory=$(dirname $2 | sed "s/theories/results/")/gg_investigation
info_file=$(dirname $2)/info.txt
base_gem_filepath=$(grep "Base GEM filepath" $info_file | cut -d ":" -f 2 | awk '{$1=$1};1')
mkdir -p $output_directory

## Extract reactions from proof
dot_path=""$output_directory"/"$1"_proof.dot"
helpers/gene_knockout_dot.sh $dot_path $1 $2 ${@:3} >/dev/null
grep -v " -> " $dot_path |
	grep rxn | perl -pe 's|^.*(r_\d{4})\w+?(reverse)?.*$|\1 \2|' | sort \
	>""$output_directory"/"$1"_proof.rxns"

# Solve FBA with the gene knockout
python3 helpers/abduction/LMtoFBA_proof.py \
	--reactions ""$output_directory"/"$1"_proof.rxns" \
	--ko_gene $1 --model $base_gem_filepath \
	--mode "promote" --threshold=1e-9 --verbose=True
	# $(grep "Base deletants" $info_file | cut -d ":" -f 2 | awk '{$1=$1};1' | sed "s/ / --ko_gene /g ; s/^/--ko_gene /") \
	# --ignore_rxn "r_3547" --ignore_rxn "r_4395" --ignore_rxn "r_0473" \
	# --ignore_rxn "r_0982" --ignore_rxn "r_0984" --ignore_rxn "r_2096" \
	# --ignore_rxn "r_3532" --ignore_rxn "r_3546"
echo >&1 $1 $hypCount $add_clause_component ""$output_directory"/"$1"_proof.rxns"
echo >&2 "Finishing investigation for "$1"."
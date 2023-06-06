#!/bin/bash
sed "/"$2"/ s/gn(g_/~gn(g_/;/"$2"/ s/_in_genome/_deletion/" $3 |
	cat - "${@:4}" |
	$IPROVER_HOME"/iproveropt" \
		--stdin=true --schedule none --preprocessing_flag false \
		--instantiation_flag false --superposition_flag false \
		--resolution_flag true --res_prop_simpl_given false \
		--res_to_prop_solver none --prop_solver_per_cl 0 \
		--proof_reduce_dot "[all_neg;unit]" \
		--proof_dot_file $1

#!/bin/zsh
source ~/.zshrc
sed "/"$1"/ s/gn(g_/~gn(g_/;/"$1"/ s/_in_genome/_deletion/" $2|cat - "${@:3}"|$IPROVER_HOME/iproveropt --stdin=true --proof_out false --sat_out_model none --sat_out_clauses false --sat_out_model none

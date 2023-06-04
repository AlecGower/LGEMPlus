#!/bin/bash
sed "/"$1"/ s/gn(g_/~gn(g_/;/"$1"/ s/_in_genome/_deletion/" $2|cat - "${@:3}"|$IPROVER_HOME/iproveropt --stdin=true
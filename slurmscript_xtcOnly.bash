#!/bin/bash

accnt=lcls
part=roma

dark=$1
bright=$2
exp=$3
for ((i=4;i<=$#;i++)); do 
	runfile="${!i}"
	echo "exp=${exp} dark=${dark} bright=${bright} runfile=${runfile}"
	if [[ "$runfile" =~ r([[:digit:]]+)-s([[:digit:]]+) ]];then 
		sbatch --account=$accnt --partition=$part --time=0-04:00:00 --job-name="bt_${BASH_REMATCH[1]}" ./runscript_xtcOnly.bash $dark $bright $exp $runfile
	fi
done


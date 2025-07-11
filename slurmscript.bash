#!/bin/bash

accnt=lcls
part=roma

dark=$1
bright=$2
exp=$3
for ((i=4;i<=$#;i++)); do 
	run="${!i}"
	echo "exp=${exp} dark=${dark} bright=${bright} run=${run}"
	sbatch --account=$accnt --partition=$part --time=0-02:00:00 --job-name="brt_$run" ./runscript.bash $dark $bright $exp $run
done


#!/bin/bash

accnt=lcls
part=roma

path=$1
for ((i=2;i<=$#;i++)); do 
	run="${!i}"
	echo "path=${path} run=${run}"
	sbatch --account=$accnt --partition=$part --time=0-02:00:00 --job-name="eig_$run" ./runeigen.bash $path $run
done


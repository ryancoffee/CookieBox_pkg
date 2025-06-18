#!/bin/bash

accnt=lcls
part=roma

dark=$1
exp=$2
for ((i=2;i<=$#;i++)); do 
	run="${!i}"
	echo "exp=$exp" "\ndark=$dark" "\nrun=$run"; 
	sbatch --account=$accnt --partition=$part --time=0-01:00:00 ./runscript.bash $dark $exp $run
done


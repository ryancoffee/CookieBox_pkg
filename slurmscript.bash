#!/bin/bash

accnt=lcls
part=roma

dark=$1
brigth=$2
exp=$3
for ((i=3;i<=$#;i++)); do 
	run="${!i}"
	echo "exp=$exp" "\ndark=$dark" "\nbright=$bright" "\nrun=$run"; 
	sbatch --account=$accnt --partition=$part --time=0-01:00:00 ./runscript.bash $dark $bright $exp $run
done


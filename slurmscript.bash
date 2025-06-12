#!/bin/bash

dark=$1
exp=$2
for ((i=2;i<=$#;i++)); do 
	run="${!i}"
	echo "exp=$exp" "\ndark=$dark" "\nrun=$run"; 
	sbatch ./runscript.bash $dark $exp $run
done


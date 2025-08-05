#!/bin/bash

exp=$1
echo "exp=$exp , runs = ${@:2}"; 
for r in ${@:2}; do
	echo $r
	sbatch --time=0-01:00:00 --account=lcls --partition=roma --job-name="dark_$r" ./rundark.bash $exp $r
done


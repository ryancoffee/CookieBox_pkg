#!/bin/bash

exp=$1
echo "exp=$exp , runs = ${@:2}"; 
for xtc in ${@:2}; do
	echo $xtc
	sbatch --time=0-01:00:00 --account=lcls --partition=roma --job-name="dk_${xtc: -17:13}" ./rundark_xtcOnly.bash $exp $xtc
done


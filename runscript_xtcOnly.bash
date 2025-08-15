#!/bin/bash

export darkpath=$1
export brightpath=$2
export exp=$3
export runfile=$4
echo "exp = $exp"
echo "runfile = $runfile"
if [[ "$runfile" =~ r([[:digit:]]+)-s([[:digit:]]+) ]]; then 
	echo "runnum ${BASH_REMATCH[1]}"
	export run=${BASH_REMATCH[1]}

	export darklist=`ls ${darkpath}/xtcav_dark_images_*.h5`
	export brightfile=${brightpath}/xtcav_bright_images_${run}.h5
	echo "using following darkfiles:"
	for v in ${darklist}; do echo $v;done
	echo "using the following brightfile"
	echo $brightfile

	source init.bash
	python3 ./python/XtcavOnly_preproc.py 

else
	echo "failed to match runnum in for runfile"

fi

#!/bin/bash

export darkpath=$1
export brightpath=$2
export exp=$3
export run=$4
echo "exp = $exp"
echo "run = $run"
export darklist=`ls ${darkpath}/xtcav_dark_images_*.h5`
export brightfile=${brightpath}/xtcav_bright_images_run${run}.h5
echo "using following darkfiles:"
for v in ${darklist}; do echo $v;done
echo "using the following brightfile"
echo $brightfile

source init.bash
python3 ./python/Xtcav_preproc.py 


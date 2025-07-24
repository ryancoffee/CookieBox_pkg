#!/bin/bash

export path=$1
export run=$4
echo "path = $path"
echo "run = $run"
export file=${path}/xtcav_bright_images_${run}.h5
echo "using the following file"
echo $file

source init_psconda2.bash
python/Xtcav_eigenimages.py $file


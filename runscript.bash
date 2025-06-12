#!/bin/bash

export dark=$1
export exp=$2
export run=$3
source init.bash
python3 ./python/Xtcav_preproc.py $dark $exp $run


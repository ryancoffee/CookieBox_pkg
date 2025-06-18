#!/bin/bash

export exp=$1
export run=$2
source init.bash
python3 ./python/Xtcav_store_dark.py $exp $run


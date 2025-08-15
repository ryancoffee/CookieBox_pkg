#!/bin/bash

export exp=$1
export xtc=$2
source init.bash
python3 ./python/XtcavOnly_store_dark.py $exp $xtc


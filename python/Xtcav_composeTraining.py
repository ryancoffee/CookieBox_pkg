#!/usr/bin/python3
import h5py
import numpy as np
import os
import random
import matplotlib.pyplot as plt
import sys
import math
import DataSelector

def main(flist):
    plotting = False
    for fname in flist:
        trainselect = DataSelector.DataSelector(32)
        trainselect.Load(fname,1<<6) #1<<9)
        if plotting:
            trainselect.PlotDistributions()
        outname,ext = os.path.splitext(fname)
        outname += '_train' + ext
        trainselect.StoreH5(outname)
    return

if __name__ == '__main__':
    if len(sys.argv)<2:
        print('give me at least one file to build into a training set')
    else:
        main(sys.argv[1:])


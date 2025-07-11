#!/sdf/group/lcls/ds/ana/sw/conda2/inst/envs/ps-4.6.3/bin/python3
'''#!/sdf/group/lcls/ds/ana/sw/conda1/inst/envs/ana-4.0.66-py3/bin/python'''

import numpy as np
import matplotlib.pyplot as plt
import random
import h5py
import sys
import os
import math

def main():
    print('REMEMBER, keep the number of bright files to a power of two for the sake of bit shift divisions\nThis will get enforced as truncating files beyond nearest power of 2\t\t... maybe')
    print('for now, running single bright file.  Think of federating the eigenimage calculation across bright files in the future')
    print('in future, also consider using the SSIM metric or Wasserstein to determine the maximally different images to include in the "training" set')
    nimages:int = 0
    batchsize:int = 1<<7
    datalist = None
    for fname in sys.argv[1:]:
        with h5py.File(fname,'r') as f:
            xt = f['xtcav']
            evtkeys = list(xt.keys())
            random.shuffle(evtkeys)
            for k in evtkeys[:1<<7]:
                if type(datalist) == type(None):
                    datalist = [xt[k]]
                else:
                    datalist += [xt[k]]
                nimages += 1
            data = np.array(datalist,dtype=int)
            print(data.shape)


    return

if __name__ == '__main__':
    if len(sys.argv)<2:
        print('give me at least one bright file')
    else:
        print(len(sys.argv))
        main()

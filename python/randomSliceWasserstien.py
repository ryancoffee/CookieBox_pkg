#!/usr/bin/python3
import h5py
import numpy as np
import os
import random
import matplotlib.pyplot as plt
import sys
import math

def main(aname:str,bname:str):
    with h5py.File(aname,'r') as a:
        with h5py.File(bname,'r') as b:
            akeys = [k for k in a['xtcav'].keys()]
            bkeys = [k for k in b['xtcav'].keys()]
            ima = a['xtcav'][akeys[0]][()]
            imb = b['xtcav'][bkeys[0]][()]
            random.shuffle(akeys)
            random.shuffle(bkeys)
            for indx in range(10):
                ima = a['xtcav'][akeys[indx]][()]
                imb = b['xtcav'][bkeys[indx]][()]
                rows = [i for i in range(ima.shape[0])]
                random.shuffle(rows)
                mse:float = float(0)
                for r in rows[:ima.shape[1]>>4]:
                    acumsum = np.cumsum(ima[r,:])
                    bcumsum = np.cumsum(imb[r,:])
                    mse += np.mean(np.power(acumsum-bcumsum,int(2)))

                    '''
                    plt.plot(acumsum,label = 'a')
                    plt.plot(bcumsum,label = 'b')
                    plt.show()
                    '''
                print(math.sqrt(mse)/(ima.shape[1]>>4))

        
        #grab two images
        
        #select random rows
        #evaluate cumsum and then MSE for matched random rows.
    return

if __name__ == '__main__':
    if len(sys.argv)<2:
        print('give me two files to work with for now')
    else:
        main(sys.argv[1],sys.argv[2])

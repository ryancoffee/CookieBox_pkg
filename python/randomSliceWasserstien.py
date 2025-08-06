#!/usr/bin/python3
import h5py
import numpy as np
import os
import random
import matplotlib.pyplot as plt
import sys
import math
import DataSelector

def main(aname:str,bname:str):
    plotting = False
    samplerows:int = 64
    maxhist = 1<<11
    rmse_hist = [0]*((maxhist>>2)+1)
    rmse_leveled_hist = [0]*((maxhist>>2)+1)
    with h5py.File(aname,'r') as a:
        with h5py.File(bname,'r') as b:
            akeys = [k for k in a['xtcav'].keys()]
            bkeys = [k for k in b['xtcav'].keys()]
            ima = a['xtcav'][akeys[0]][()]
            imb = b['xtcav'][bkeys[0]][()]
            random.shuffle(akeys)
            random.shuffle(bkeys)
            for indx in range(min(len(akeys),len(bkeys))):
                ima = a['xtcav'][akeys[indx]][()]
                imb = b['xtcav'][bkeys[indx]][()]
                rows = [i for i in range(ima.shape[0])]
                random.shuffle(rows)
                mse:float = float(0)
                for r in rows[:samplerows]:
                    acumsum = np.cumsum(ima[r,:])
                    bcumsum = np.cumsum(imb[r,:])
                    mse += np.mean(np.power(acumsum-bcumsum,int(2)))

                    '''
                    plt.plot(acumsum,label = 'a')
                    plt.plot(bcumsum,label = 'b')
                    plt.show()
                    '''
                updateind = min(np.uint16(math.sqrt(mse)/samplerows)>>2,maxhist>>2)
                rmse_hist[updateind] += 1
                if rmse_leveled_hist[updateind]==0:
                    rmse_leveled_hist[updateind] += 1
                else:
                    rmse_leveled_hist[updateind] += int(random.uniform(0,1)<1.0/rmse_leveled_hist[updateind])


        

        if plotting:
            histbins = [(1<<2)*i for i in range((maxhist>>2)+2)]
            plt.stairs(rmse_hist,histbins,label='histogram')
            plt.stairs(rmse_leveled_hist,histbins,label='leveled histogram')
            plt.xlabel('rmse value')
            plt.ylabel('occurence')
            plt.legend()
            plt.show()
        #grab two images
        #select random rows
        #evaluate cumsum and then MSE for matched random rows.
    return

if __name__ == '__main__':
    if len(sys.argv)<2:
        print('give me two files to work with for now')
    else:
        main(sys.argv[1],sys.argv[2])


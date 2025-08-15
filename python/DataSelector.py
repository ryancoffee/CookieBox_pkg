#!/usr/bin/python3
import numpy as np
import random
import h5py
import time
import hashlib
import math
import matplotlib.pyplot as plt
import utils

class DataSelector:
    def __init__(self,subsamplerows:int=1<<6,style:str='randomWasserstein'):
        random.seed(time.time_ns()%(1<<10))
        self.style = style
        self.subsamplerows:int = subsamplerows
        self.measurerows = []
        self.trainSet = []
        self.testSet = []
        self.trainKeys = []
        self.testKeys = []
        self.trainMeasures = []
        self.testMeasures = []
        self.maxhist = 1<<11
        self.rmse_hist = [0]*(self.maxhist)
        self.rmse_leveled_hist = [0]*(self.maxhist)
        return

    def Load(self,fname:str,nsamples:int=1<<10,shuffle=True): # fill a training set of images up to nsamples deep
        print('loading %s'%fname)
        with h5py.File(fname,'r') as f:
            allkeys = [k for k in f['xtcav'].keys()]
            if shuffle:
                random.shuffle(allkeys)
            
            self.trainSet += [f['xtcav'][allkeys[0]][()]]
            self.trainKeys += [allkeys[0]]
            self.trainMeasures = [0]
            rows = [i for i in range(self.trainSet[0].shape[0])]
            random.shuffle(rows)
            self.subsamplerows = min(self.subsamplerows,self.trainSet[0].shape[0]>>2) # use at most 1/4 of the image for measuring Wasserstein
            self.measurerows = rows[:self.subsamplerows]
            for k in allkeys[1:]: # this runs from 1: since we used the 0th already to start the list.
                if len(self.trainSet)==nsamples: # first correct neighbor distance measure for the 0th image and then break
                    self.trainMeasures[0] = self.getNeighborDistance_zeroth(self.trainSet[0])
                    self.rmse_leveled_hist[ min(self.trainMeasures[0],len(self.rmse_leveled_hist)-1) ] += 1
                    break
                im = f['xtcav'][k][()]
                d = self.getNeighborDistance(im)
                self.rmse_hist[ min(d,len(self.rmse_hist)-1) ] += 1
                if random.uniform(0.0,1.0)<(1.0/(1.0+self.rmse_hist[min(d,len(self.rmse_hist)-1)])):
                    self.trainSet += [f['xtcav'][k][()]]
                    self.trainKeys += [k]
                    self.trainMeasures += [d]
                    self.rmse_leveled_hist[ min(d,len(self.rmse_leveled_hist)-1) ] += 1
                    l = len(self.trainSet)
                    if (utils.bit_count(l)==1 or utils.bit_count(l%64)==0):
                        print('\ttraining set legth = %i'%l)
        return self


    def StoreH5(self,outname:str):
        grpstr = 'train_%s'%(hashlib.blake2b(time.asctime().encode(),digest_size=2).hexdigest()) 
        print('Writing to %s training set %s'%(outname,grpstr) )
        with h5py.File(outname,'a') as o:
            if 'xtcav' not in o.keys():
                o.create_group('xtcav')
            if grpstr in o['xtcav'].keys():
                del o['xtcav'][grpstr]
            grp = o['xtcav'].create_group(grpstr)
            grp.create_dataset('images',data=self.trainSet)
            grp.create_dataset('keys',data=self.trainKeys)
            grp.create_dataset('measures',data=self.trainMeasures)
            grp.create_dataset('measurerows',data=self.measurerows)
            grp.create_dataset('rmseHist',data=self.rmse_hist)
            grp.create_dataset('rmseLeveledHist',data=self.rmse_leveled_hist)
        return self

    def getNeighborDistance_zeroth(self,im): # return a the nearest neighbor distance, but only used for correcting the 0th
        d:int = 1<<16
        for sample in self.trainSet[1:]:
            mse = 0.0
            for r in self.measurerows:
                mse += np.mean(np.power(np.cumsum(im[r,:])-np.cumsum(sample[r,:]),int(2)))
            tmp = int(math.sqrt(mse)/float(len(self.measurerows)))
            if tmp < d:
                d = tmp
        return d

    def getNeighborDistance(self,im): # return a the nearest neighbor distance
        d:int = 1<<16
        for sample in self.trainSet:
            mse = 0.0
            for r in self.measurerows:
                mse += np.mean(np.power(np.cumsum(im[r,:])-np.cumsum(sample[r,:]),int(2)))
            tmp = int(math.sqrt(mse)/float(len(self.measurerows)))
            if tmp < d:
                d = tmp
        return d

    def PlotDistributions(self):
        histbins = [i for i in range(self.maxhist+1)]
        plt.stairs(self.rmse_hist,histbins,label='histogram')
        y = np.array(self.rmse_leveled_hist) + 10
        plt.stairs(y,histbins,label='leveled histogram')
        plt.xlabel('rmse value')
        plt.ylabel('occurence')
        plt.legend()
        plt.show()
        return self

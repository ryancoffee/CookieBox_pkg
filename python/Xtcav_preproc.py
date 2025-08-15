#!/sdf/group/lcls/ds/ana/sw/conda1/inst/envs/ana-4.0.66-py3/bin/python

import psana
import sys
import matplotlib.pyplot as plt
import numpy as np
import cv2
import utils
import os
import re
import h5py
import math
import random

plotting = False # use this to supress plotting

def main():
    darklist = os.environ['darklist'].split('\n')
    darkruns = [int(re.search(r'xtcav_dark_images_(\d+).h5$',n).group(1)) for n in darklist]
    darkpath = os.environ['darkpath']

    if not os.path.exists(os.path.join(darkpath,darklist[0])):
        print('not a single dark file at all')
        print(os.path.join(darkpath,darklist[0]))
        return

    
    #with h5py.File(darkfilename,'r') as df: 
        #darkruns = [int( re.search(r'(\d+)$',k).group(1) ) for k in df.keys()]

    print('using darkruns ' + ' '.join([str(r) for r in darkruns]))

    runstr = os.environ['run']
    dsourcestr = 'exp=' + os.environ['exp'] + ':run=' + runstr
    runnum = int(runstr)
    darkind = np.argmin(np.abs(runnum-np.array(darkruns)))
    darkname = darklist[darkind]
    with h5py.File(os.path.join(darkpath,darklist[darkind])) as df:

        ## this is always starting the seed as the runnum, so sequence is the same every time ##
        random.seed(runnum)
        darkkey = [k for k in df.keys() if re.search(r'\d+',k)][0]
        avimkeys = [k for k in df[darkkey].keys()]
        random.shuffle(avimkeys)
        _=[print(v) for v in avimkeys] 
        darkimg = np.zeros(df[darkkey][avimkeys[0]][()].shape,dtype=np.uint16)
        pow2 = int(math.log2(len(avimkeys)))
        for k in avimkeys[:(1<<pow2)]:
            darkimg += df[darkkey][k][()]
        darkimg >>= (pow2-1) 

        if plotting:
            plt.imshow(darkimg,origin='lower')
            plt.colorbar()
            plt.show()

        print("running:\t",dsourcestr)
        ds = psana.DataSource(dsourcestr)
        print(psana.DetNames('detectors'))
        xt = psana.Detector('xtcav')
        #for nevent,evt in enumerate(ds.events()):
        nevent = 0
        nlimit = 1<<16

        print("opening " + os.environ['brightfile'])
        with h5py.File(os.environ['brightfile'],'w') as of:
            grp = of.create_group('xtcav')
            of.create_dataset('dark',data=darkimg,dtype=np.uint16)
                

            while nevent < nlimit:
                if utils.bit_count(nevent)==1: # ^ (utils.bit_count(nevent)==2 and utils.bit_count(nevent&(1<<8))==1):
                    print('working event: %i \tfor %s'%(nevent,dsourcestr))
                #_ = [next(ds.events) for i in range(nskip)]
                evt = next(ds.events())
                raw = xt.raw(evt)
                if raw is not None:
                    ny,nx = utils.getshape(raw)
                    im = np.zeros((ny,nx),dtype=int)
                    for i in range(ny):
                        row = raw[i]
                        for j in range(nx):
                            im[i,j] = row[j]
                    img = np.array(im,dtype=np.int16)
                    img = cv2.erode(img,kernel=(3,3))
                    img = cv2.blur(img,ksize=(3,3))
                    img -= darkimg.astype(np.int16)
                    boolimg = np.zeros((ny,nx),dtype=bool)
                    thresh = max(np.max(img)>>3,1<<7)
                    threshinds = np.where(img>thresh)
                    # print('nx = %i'%nx)
                    # print('ny = %i'%ny)
                    # print(threshinds[0])
                    # horizontal is ind 1, vertical is ind 0
                    if len(threshinds[0]) == 0:
                        vert = ny>>1
                    else:
                        vert = sum(threshinds[0])//len(threshinds[0])
                    if len(threshinds[1])==0:
                        horiz = nx>>1
                    else:
                        horiz = sum(threshinds[1])//len(threshinds[1])

                    rollv = (ny>>1) - vert
                    rollh = (nx>>1) - horiz
                    img = np.roll(img,rollv,axis=0)
                    img = np.roll(img,rollh,axis=1)

                    evtimg = grp.create_dataset('evt_%08i'%(nevent),data=img,dtype=np.int16)
                    #evtimg.attrs.create('vertthresh',data=threshinds[0])
                    #evtimg.attrs.create('horizthresh',data=threshinds[1])
                    evtimg.attrs.create('centroid',data=(horiz,vert))
                    evtimg.attrs.create('roll_img',data=(rollh,rollv))
                    evtimg.attrs.create('erode_blur_kernel',data=(3,3))


    
                    # boolimages[-1][threshinds] = True
                        
                        
                    if plotting:
                        plt.imshow(img,origin='lower')
                        plt.colorbar()
                        plt.clim(0,1<<9)
                        plt.show()
                        
                nevent += 1
    return

if __name__ == '__main__':

    m = re.search(r'^(\w{3})(.*)$',os.environ['exp'])
    if m:
        print('hutch = %s'%m.group(1))
        print('expnum = %s'%m.group(2))
        hutchname = m.group(1)

    print(os.environ['darkpath'])
    print(os.environ['brightpath'])

    m = re.search(r'^/sdf/data/lcls/.*',os.environ['darkpath'])
    if m:
        print('dark image path = \t%s'%(m.group(0)))
    else:
        print('no dark path given')
        os.environ['darkpath'] = 'none'

    m = re.search(r'^/sdf/data/lcls/.*',os.environ['brightpath'])
    if m:
        print('bright image path = \t%s'%(m.group(0)))
    else:
        print('no dark path given')
        os.environ['brightpath'] = 'none'


    # /sdf/data/lcls/ds/amo/amo86815/scratch/
    #os.system('echo $brightpath')
    main()

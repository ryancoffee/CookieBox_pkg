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

def main(expstr,runs):
    darkfilename = os.environ['darkfile']
    m = re.search(r'dark.*h5$',darkfilename)
    if not m:
        print('no dark file given, aborting')
        return

    
    with h5py.File(darkfilename,'r') as df: 
        darkruns = [int( re.search(r'(\d+)$',k).group(1) ) for k in df.keys()]
        print('using darkruns ' + ' '.join([str(r) for r in darkruns]))

        for r in (runs):
            dsourcestr = 'exp=' + expstr + ':run=' + r
            runnum = int(r)
            darkind = np.argmin(np.abs(runnum-np.array(darkruns)))
            darkkey = list(df.keys())[darkind]

            ## this is always starting the seed as the runnum, so sequence is the same every time ##
            random.seed(runnum)
            avimkeys = [k for k in df[darkkey].keys()]
            random.shuffle(avimkeys)
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
            nlimit = 1<<14

            with h5py.File(os.environ['outfile'],'w') as of:
                grp = of.create_group('xtcav')

                while nevent < nlimit:
                    if utils.bit_count(nevent)==1: # ^ (utils.bit_count(nevent)==2 and utils.bit_count(nevent&(1<<8))==1):
                        print('working event: %i'%(nevent))
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
    hutchname = 'tmo'
    if len(sys.argv)<4:
        runstrs = ['69','70']
        expstr = str('amo86815') #str('amoi0314');
        print('running %s: runs '%(expstr) + ' '.join(runstrs))
    else:
        expstr = str(sys.argv[2])
        runstrs = [str(v) for v in sys.argv[3:]]
        print('running %s: runs '%(expstr) + ' '.join(runstrs))
    m = re.search(r'^(\w{3})(.*)$',expstr)
    if m:
        print('hutch = %s'%m.group(1))
        print('expnum = %s'%m.group(2))
        hutchname = m.group(1)

    m = re.search(r'^/sdf/.*dark.*\.h5$',sys.argv[1])
    if m:
        print('dark image filename = \t%s'%(m.group(0)))
        os.environ['darkfile'] = m.group(0)
    else:
        print('no dark image given')
        os.environ['darkfile'] = 'none'


    # /sdf/data/lcls/ds/amo/amo86815/scratch/
    outpath = '/sdf/data/lcls/ds/' + hutchname + '/' + expstr + '/scratch/xtcav_bright_images_' + '_'.join(runstrs) + '.h5'
    os.environ['outfile'] = outpath
    os.system('echo $outfile')
    main(expstr,runstrs)

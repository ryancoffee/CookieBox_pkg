#!/sdf/group/lcls/ds/ana/sw/conda1/inst/envs/ana-4.0.66-py3/bin/python

import psana
import sys
import matplotlib.pyplot as plt
import numpy as np
import h5py
import utils 
import re
import os

def main(expstr,runs):
    fname = os.getenv('outfile')
    with h5py.File(fname,'a') as f:
        rungrp = None
        for r in (runs):
            rkey = 'run_%03i'%int(r)
            if rkey in list(f.keys()):
                rungrp = f[rkey]
            else:
                rungrp = f.create_group(rkey)
            dsourcestr = 'exp=' + expstr + ':run=' + r
            print("running:\t",dsourcestr)
            print('Need to subtract a running estimate of dark imates')
            ds = psana.DataSource(dsourcestr)
            print(psana.DetNames('detectors'))
            xt = psana.Detector('xtcav')
            avgimg = None #np.zeros((1,1),dtype=np.uint16)
            nevent = 0
            stepevent = 0

            # seems to be the limit is 1<<19 shots in dark runs
            nlimit = 1<<12
            navglim_shift = 7
            nskip = 0
            avgidx = 0
            while nevent < nlimit:
                try:
                    evt = next(ds.events())
                    _ = [next(ds.events()) for i in range(nskip)]
                    raw = xt.raw(evt)
                    if raw is not None:
                        nx,ny = utils.getshape(raw)
                        im = np.zeros((nx,ny),dtype=np.uint16)
                        if (avgimg is None) or (stepevent == 0):
                            avgimg = np.zeros((nx,ny),dtype=np.uint16)
                        for i in range(nx):
                            row = raw[i]
                            for j in range(ny):
                                im[i,j] = row[j]
                        im = utils.massage(im)
                        if nevent%(1<<navglim_shift) == 0:
                            avgimg = np.copy(im)
                        else:
                            avgimg += im

                    nevent += 1
                    stepevent += 1
    
                    if nevent%(1<<navglim_shift) == 0:
                        avgimg >>= navglim_shift 
                        akey = 'avg_%03i'%avgidx
                        if akey in list(rungrp.keys()):
                            del rungrp[akey]
                        print('storing %s'%akey)
                        avimds = rungrp.create_dataset(akey,data=avgimg)
                        avimds.attrs['avgindex'] = avgidx
                        avimds.attrs['navg'] = 1<<navglim_shift
                        avgidx += 1
                        stepevent = 0

                
                except:
                    next
    return

if __name__ == '__main__':
    hutchname = 'tmo'
    if len(sys.argv)<3:
        runstrs = ['69','70']
        expstr = str('amo86815') #str('amoi0314');
        print('running %s: runs '%(expstr) + ' '.join(runstrs))
    else:
        expstr = str(sys.argv[1])
        runstrs = [str(v) for v in sys.argv[2:]]
        print('running %s: runs '%(expstr) + ' '.join(runstrs))
    m = re.search(r'^(\w{3})(.*)$',expstr)
    if m:
        print('hutch = %s'%m.group(1))
        print('expnum = %s'%m.group(2))
        hutchname = m.group(1)

    # /sdf/data/lcls/ds/amo/amo86815/scratch/
    outpath = '/sdf/data/lcls/ds/' + hutchname + '/' + expstr + '/scratch/xtcav_dark_images_run' + '_'.join(runstrs) + '.h5'
    os.environ['outfile'] = outpath
    os.system('echo $outfile')
    main(expstr,runstrs)

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
        for r in (runs):
            dsourcestr = 'exp=' + expstr + ':run=' + r
            print("running:\t",dsourcestr)
            print('Need to subtract a running estimate of dark imates')
            ds = psana.DataSource(dsourcestr)
            print(psana.DetNames('detectors'))
            xt = psana.Detector('xtcav')
            images = []
            nevent = 0
            nlimit = 10
            nskip = 100
            while nevent < nlimit:
                evt = next(ds.events())
                _ = [next(ds.events()) for i in range(nskip)]
                raw = xt.raw(evt)
                if raw is not None:
                    nx,ny = utils.getshape(raw)
                    im = np.zeros((nx,ny),dtype=np.uint16)
                    for i in range(nx):
                        row = raw[i]
                        for j in range(ny):
                            im[i,j] = row[j]
                    img = np.copy(im)
                    img = utils.massage(img)
                    images += [img]
                nevent += 1
                 
                plt.imshow(images[-1],origin='lower')
                plt.show()
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
    outpath = '/sdf/data/lcls/ds/' + hutchname + '/' + expstr + '/scratch/xtcav_dark_images_' + '_'.join(runstrs) + '.h5'
    os.environ['outfile'] = outpath
    os.system('echo $outfile')
    main(expstr,runstrs)

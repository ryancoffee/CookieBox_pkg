#!/sdf/data/lcls/ds/amo/amox28216/scratch/all_output/

import psana
import sys
import matplotlib.pyplot as plt
import numpy as np


def getshape(raw):
    nx = len(raw)
    ny = len(raw[0])
    print(nx,ny)
    return nx,ny

def main(expstr,runs):
    for r in (runstrs):
        dsourcestr = 'exp=' + expstr + ':run=' + r
        print("running:\t",dsourcestr)
        ds = psana.DataSource(dsourcestr)
        print(psana.DetNames('detectors'))
        xt = psana.Detector('xtcav')
        images = []
        #for nevent,evt in enumerate(ds.events()):
        nevent = 0
        while nevent < 10:
            evt = next(ds.events())
            raw = xt.raw(evt)
            if raw is not None:
                nx,ny = getshape(raw)
                im = np.zeros((nx,ny),dtype=int)
                for i in range(nx):
                    row = raw[i]
                    for j in range(ny):
                        im[i,j] = row[j]
                images += [im]
            nevent += 1
                    
            plt.imshow(images[0],origin='lower')
            plt.show()



            
    return

if __name__ == '__main__':
    if len(sys.argv)<3:
        runstrs = ['69','70']
        expstr = str('amo86815') #str('amoi0314');
        print('running %s: runs '%(expstr) + ' '.join(runstrs))
    else:
        expstr = str(sys.argv[1])
        runstrs = [str(v) for v in sys.argv[2:]]
        print('running %s: runs '%(expstr) + ' '.join(runstrs))
    main(expstr,runstrs)

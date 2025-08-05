#!/sdf/group/lcls/ds/ana/sw/conda1/inst/envs/ana-4.0.66-py3/bin/python

import psana
import sys
import matplotlib.pyplot as plt
import numpy as np
import cv2


def getshape(raw):
    nx = len(raw)
    ny = len(raw[0])
    print(nx,ny)
    return nx,ny

def main(expstr,runs):
    for r in (runs):
        dsourcestr = 'exp=' + expstr + ':run=' + r
        print("running:\t",dsourcestr)
        print('Need to subtract a running estimate of dark imates')
        ds = psana.DataSource(dsourcestr)
        print(psana.DetNames('detectors'))
        xt = psana.Detector('xtcav')
        images = []
        residuals = []
        #for nevent,evt in enumerate(ds.events()):
        nevent = 0
        nlimit = 10
        nskip = 100
        while nevent < nlimit:
            #_ = [next(ds.events) for i in range(nskip)]
            evt = next(ds.events())
            raw = xt.raw(evt)
            if raw is not None:
                nx,ny = getshape(raw)
                im = np.zeros((nx,ny),dtype=int)
                for i in range(nx):
                    row = raw[i]
                    for j in range(ny):
                        im[i,j] = row[j]
                img = np.array(im,dtype=np.uint16)
                img = cv2.erode(img,kernel=(3,3))
                img = cv2.blur(img,ksize=(3,3))
                images += [img]
                residuals += [np.array(im,dtype=np.uint16) - img]
            nevent += 1
                    
            plt.imshow(images[-1],origin='lower')
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

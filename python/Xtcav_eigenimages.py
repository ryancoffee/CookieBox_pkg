#!/sdf/group/lcls/ds/ana/sw/conda2/inst/envs/ps_20241122/bin/python3
'''#!/sdf/group/lcls/ds/ana/sw/conda1/inst/envs/ana-4.0.66-py3/bin/python'''

import numpy as np
import matplotlib.pyplot as plt
import random
import h5py
import sys
import os
import math
import cv2

def oversample(im,shape,factor=4,plotting=False):
    lv,lh = shape
    newshape = (factor*lv,factor*lh)
    res = np.zeros(newshape,dtype=type(im[0]))
    RES = np.zeros(newshape,dtype=complex)
    IM = np.fft.fft2(im.reshape(shape))
    if plotting:
        plt.imshow(np.fft.ifft2(IM).real)
        plt.show()
        print(IM.shape)
        print(RES.shape)
    RES[:lv>>1,:lh>>1] = IM[:lv>>1,:lh>>1]
    RES[-(lv>>1):,:lh>>1] = IM[-(lv>>1):,:lh>>1]
    RES[:lv>>1,-(lh>>1):] = IM[:lv>>1,-(lh>>1):]
    RES[-(lv>>1):,-(lh>>1):] = IM[-(lv>>1):,-(lh>>1):]
    res = np.fft.ifft2(RES).real
    if plotting:
        plt.imshow(res)
        plt.show()
    return res

def main():
    print('REMEMBER, keep the number of bright files to a power of two for the sake of bit shift divisions\nThis will get enforced as truncating files beyond nearest power of 2\t\t... maybe')
    print('for now, running single bright file.  Think of federating the eigenimage calculation across bright files in the future')
    print('in future, also consider using the SSIM metric or Wasserstein to determine the maximally different images to include in the "training" set')
    print('Good god, I have to crop into the image significantly and change stride to 4 in both dimensions in order to fit the memory on iana!  Otherwise in the tens of gigs.')
    nimages:int = 0
    batchsize:int = 1<<10
    datalist = None
    lenv = 1
    lenh = 1
    for fname in sys.argv[1:]:
        with h5py.File(fname,'r') as f:
            xt = f['xtcav']
            evtkeys = list(xt.keys())
            random.shuffle(evtkeys)
            for k in evtkeys[:batchsize]:
                im = cv2.blur(np.array(xt[k][64:-64,200:-200]),ksize=(3,3))[::4,::4]
                (lenv,lenh) = im.shape
                if type(datalist) == type(None):
                    datalist = [im.flatten()]
                else:
                    datalist += [im.flatten()]
                nimages += 1
            data = np.array(datalist,dtype=int)
            avgimg = np.mean(data,axis=0)
            resdata = np.array([d - avgimg.astype(int) for d in datalist],dtype=int)
            eigvecs,eigvals,_ = np.linalg.svd(resdata.T)
            '''
            print(data.shape)
            plt.imshow(data[-1].reshape((lenv,lenh)),origin='lower')
            plt.colorbar()
            plt.show()
            print(np.max(data),np.min(data))
            plt.imshow(avgimg.reshape((lenv,lenh)))
            plt.colorbar()
            plt.title('lenv = %i, lenh = %i'%(lenv,lenh))
            plt.show()

            plt.loglog(eigvals,'.')
            plt.xlabel('eigenimage')
            plt.ylabel('eigenvalue')
            plt.savefig('../eigenimages/eigenvalues.png')
            plt.show()
            for i in range(1<<4):
                plt.imshow((eigvecs[i,:]).reshape((lenv,lenh)))
                plt.title('eigenimage %i'%(i))
                plt.savefig('../eigenimages/eigenimage_%03i.png'%(i))
            '''

            with h5py.File('../eigenimages/eigenfile.h5','a') as out:
                if fname in out.keys():
                    del out[fname]
                grp = out.create_group(fname)
                grp.create_dataset('eigenvalues',data=eigvals)
                grp.create_dataset('avgimage',data=avgimg.reshape((lenv,lenh)))
                grp.create_dataset('avgimage_oversampled',data=oversample(avgimg,shape=(lenv,lenh),factor=4,plotting=False))
                for i in range(1<<8):
                    eigim = grp.create_dataset('eigenimage_%03i'%i,data=eigvecs[i,:].reshape((lenv,lenh))) 
                    eigim.attrs.create('note',data='# this is undersampled by factor of 4, so really could inflate back wiht fft oversample method.')
                    eigimover = grp.create_dataset('eigenimage_oversampled_%03i'%i,data=oversample(eigvecs[i,:],shape=(lenv,lenh),factor=4,plotting=False)) 
                    eigimover.attrs.create('note',data='# this is artificially oversampled by padding with zeros to get up to the original shape before the undersampling for sake of SVD')

            for k in evtkeys[batchsize:batchsize<<1]:
                im = np.array(xt[k][64:-64:4,200:-200:4]).flatten()
                coeffs = [np.inner(im,eigvecs[i,:]) for i in range(1<<8)]
                plt.scatter(np.arange(len(coeffs)),coeffs,label=k)
            plt.show()
            



    return

if __name__ == '__main__':
    if len(sys.argv)<2:
        print('give me at least one bright file')
    else:
        print(len(sys.argv))
        main()

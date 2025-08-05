#!/usr/bin/python3

import numpy as np
import sys
import re

def loadfile(fname):
    return np.loadtxt(fname)*1e-6

def savefile(fname,data):
    np.savetxt(fname,data,fmt='%f')
    return

def convfilter(nbins):
    f = np.fft.fftfreq(nbins)
    t0 = 20.
    gamma = 1./17
    bwd = 0.3
    Q = 1./(gamma + 1j*f)*np.exp(-1j*f*2.*np.pi*t0)*np.exp(-(f/bwd)**2)
    invQ = (gamma + 1j*f)*np.exp(-1j*f*2.*np.pi*t0)*np.exp((f/bwd)**2)
    q = np.fft.ifft(Q)
    return (Q,invQ,q)

def lowhighpass(f,bwd):
    lowfilt = np.zeros(f.shape[0],dtype=float)
    highfilt = np.zeros(f.shape[0],dtype=float)
    linds = np.where(np.abs(f)<bwd)
    highbwd = 2*bwd
    hinds = np.where(np.abs(f)<highbwd)
    c2=np.power(np.cos(f/bwd*np.pi/2.),int(2))
    lowfilt[linds] = c2[linds]
    s2 = np.power( np.sin((np.abs((f)/highbwd)*np.pi)) , int(2) ) 
    #highfilt[hinds] = s2[hinds]
    highfilt = 1-lowfilt
    np.savetxt('filters.dat',np.column_stack((f,c2,s2,lowfilt,highfilt)))
    return lowfilt,highfilt

def main():
    if len(sys.argv)< 2:
        print('need a filelist to analyze')
        return
    filelist = sys.argv[1:]
    data = loadfile(filelist[0])
    sumdata = np.zeros(data.shape)
    sumresult = np.zeros(data.shape)
    nsums = 0
    f = np.fft.fftfreq(data.shape[0],1./data.shape[0])
    bwd = data.shape[0]/3.
    (Q,iQ,q) = convfilter(data.shape[0])
    Q = np.fft.fft(q)
    LOGABSQ = np.fft.fft(np.log(np.abs(q)))
    DLOGABSQ = 1j*f*np.fft.fft(np.log(np.abs(q)))
    QMAT = np.tile(Q,(data.shape[1],1)).T
    LOGABSQMAT = np.tile(LOGABSQ,(data.shape[1],1)).T
    DLOGABSQMAT = np.tile(DLOGABSQ,(data.shape[1],1)).T
    iQMAT = np.tile(iQ,(data.shape[1],1)).T
    lpf,hpf = lowhighpass(f,bwd)
    LPFMAT = np.tile( lpf,(data.shape[1],1)).T
    HPFMAT = np.tile( hpf,(data.shape[1],1)).T
    FMAT = np.tile(f,(data.shape[1],1)).T
    DQMAT = 1j*FMAT * np.tile(Q,(data.shape[1],1)).T
    for fname in filelist:
        m = re.search('\d$',fname)
        if m:
            nsums +=1
            data = loadfile(fname)
            sumdata += data
            DATA = np.fft.fft(data,axis=0)
            signdata = np.sign(data)
            #logabsdata = np.log(np.abs(data)+1)
            #LOGABSDATA = np.fft.fft(logabsdata,axis=0) 
            DATALPF = np.copy(DATA) * LPFMAT
            DDATALPF = 1j * FMAT* np.copy(DATA) * LPFMAT 
            DDDATALPF = -1 * FMAT*FMAT* np.copy(DATA) * LPFMAT 
            #LOGABSDATAHPF = np.copy(LOGABSDATA) * HPFMAT
            #recover = np.exp( np.fft.ifft(DLOGABSDATALPF,axis=0).real + np.fft.ifft(DLOGABSQMAT,axis=0).real ) #-LOGABSQMAT,axis=0).real ) -1. # * signdata
            #remove = np.exp( np.fft.ifft(LOGABSDATAHPF,axis=0).real) -1 #-LOGABSQMAT,axis=0).real ) -1. # * signdata
            #recover = ( np.fft.ifft(LOGABSDATA-LOGABSQMAT,axis=0).real ) # * signdata
            #recover = np.fft.ifft(DATA*QMAT,axis=0).real 
            #oname = fname + '.recover'
            #savefile(oname, recover)
            #oname = fname + '.remove'
            #savefile(oname, remove)
            oname = fname + '.result'
            datalpf = np.fft.ifft(DATALPF,axis=0).real
            inds = np.where(datalpf < 0)
            datalpf[inds] = 0.
            savefile(oname, datalpf )#* np.fft.ifft(DQMAT,axis=0).real)
            oname = fname + '.result2'
            ddatalpf =  - np.fft.ifft(DDDATALPF,axis=0).real
            inds = np.where(ddatalpf<0)
            ddatalpf[inds] = 0.
            result = datalpf*ddatalpf
            finalresult = np.zeros(result.shape)
            thresh = 5e-5
            inds = np.where(result>thresh)
            finalresult[inds] = np.log(1/thresh*result[inds])
            sumresult += finalresult
            savefile(oname, finalresult)
    sumdata /= nsums
    sumresult /= nsums
    oname = 'average_waveforms.dat'
    np.savetxt(oname,sumdata)
    oname = 'average_results.dat'
    np.savetxt(oname,sumresult)
    return



if __name__ == "__main__":
    main()

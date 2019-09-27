#!/usr/bin/python3

import numpy as np
import sys
import re

#newprojections_2019_sept_05x05/projections-amoi0314-r0164_e%i_g%i_l%i_tspan1.083.dat
def main():
    if len(sys.argv) < 2:
        print('give me some files to process, e.g. ./getprojections.py newprojections_2019_sept_05x05/projections-amoi0314-r0164_e*_g1_l2_tspan1.083.dat')
        return
    rotorlist = [0,1,2]
    enlist = []
    fnames = sys.argv[1:]
    datastore = np.zeros((0,0,0),dtype = float)
    l = 0
    datadir = './'
    procdir = './'
    for fname in fnames: 
        m = re.match('(.+/)projections(-amoi0314-r\d+)_e(\d+)_g(\d+)_l(\d+)_tspan(.+).dat',fname)
        if m:
            data = np.loadtxt(fname,dtype=float)
            if datastore.shape[1]<1:
                datastore.resize((len(rotorlist),len(fnames),data.shape[0]))

            datadir = m.group(1)
            procdir = datadir + 'processed/'

            exprunname = m.group(2)
            e = int(m.group(3))
            enlist += [e] 
            g = int(m.group(4))
            l = int(m.group(5))

            for r,rotor in enumerate(rotorlist):
                datastore[r,e,:] = data[:,r]

    for r,rotor in enumerate(rotorlist):
        procout = procdir + 'projection{}{}_g{}_l{}.dat'.format(rotor,exprunname,g,l) 
        np.savetxt(procout,datastore[r,:,:],fmt='%.3e')
        print(procout)
    return

if __name__ == '__main__':
    main()

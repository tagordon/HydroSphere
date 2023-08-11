import numpy as np
import os

inds, rho, f = np.loadtxt('key.dat').T
inds = np.int64(inds)

f0 = [0.005, 0.01, 0.015, 0.02, 0.025, 0.03]
outdir = 'output_g'

for mf in f0:
    ids = inds[np.where(abs(f - mf) == np.min(np.abs(f - mf)))]
    files = ['out' + str(i) + '.dat' for i in ids]
    zprofile = []
    for i, file in enumerate(files):
        z, dens, alpha, Cp, dT_dz, phase, T, P, g = np.loadtxt(outdir + '/' + file).T
        zprofile.append(phase)
   
    zprofile = np.array(zprofile)
    #zprofile = zprofile[np.argsort(rho[inds])] 
    np.savetxt('zprofile{0}'.format(np.int64(1e3 * mf)) + '_g.out', zprofile)

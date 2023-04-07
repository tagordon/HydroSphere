import numpy as np
import os

i, rho, f = np.loadtxt('key.dat').T

rho0 = 0.9
outdir = 'output'

inds = i[np.where(rho - rho0 == np.min(rho - rho0))]
files = ['out' + str(ind) + '.dat' for ind in inds]

zprofile = np.empty(len(files), dtype=np.ndarray)

for i, f in enumerate(files):
    z, rho, alpha, Cp, dT_dz, phase, T, P, g = np.loadtxt(outdir + '/' + file).T
    zprofile[i] = z
    
np.savetxt('zprofile.out', zprofile)
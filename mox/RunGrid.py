import numpy as np
import sys
sys.path.append('/gscratch/astro/tagordon/hydrosphere')
from HydroSphere import HydroSphere
from multiprocessing import Pool, cpu_count

# helper function to compute structure as a function of bulk density and equilibrium surface temperature
# assuming core density, water fraction, radius, and surface pressure
def compute_structure(
    rho, 
    T, 
    fw,
    Mc,
    Ps, 
    resolution=100, 
    M_thresh=0.01, 
    savefile=None
):
    
    rho *= 5514
    Mc *= 5.972e24
    Mw = Mc / (1 / fw - 1)
    
    # estimate Rc given estimated density of water layer 
    Rc = ((3./4.) * Mc / (np.pi * rho)) ** (1./3.)

    res = np.array(HydroSphere(Ps, T, Mw, Mc, rho))
    
    if savefile is not None:
        
        # columns are z, rho, alpha, Cp, dT_dz, phase, T, P, g
        np.savetxt(savefile, res.T)
    else:
        return res

rho = float(sys.argv[1])
Mc = float(sys.argv[2])
Ps = 0.1 # MPa 
outpath = sys.argv[4]
    
def f(c):
    i, fw, T = c
    savefile = outpath + '/out{0}.dat'.format(i)
    compute_structure(rho, T, fw, Mc, Ps, savefile=savefile)
    
N = int(sys.argv[3])
fw = np.linspace(0.0001, 0.025, N)
T = np.linspace(180, 300, N)
fw, T = np.meshgrid(fw, T)
coords = [(i, fw, T) for i, (fw, T) in enumerate(zip(fw.flatten(), T.flatten()))]
np.savetxt(outpath + '/key.dat', coords)

if __name__ == '__main__':
    Pool(cpu_count()).map(f, coords)

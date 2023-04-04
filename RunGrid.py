import numpy as np
from HydroSphere import HydroSphere
from multiprocessing import Pool

# helper function to compute structure as a function of bulk density and equilibrium surface temperature
# assuming core density, water fraction, radius, and surface pressure
def compute_structure(
    rho_bulk, 
    Teq, 
    rho_core, 
    fwater, 
    R, 
    Psurf, 
    resolution=100, 
    M_thresh=0.01, 
    savefile=None
):
    
    rho_bulk *= 5514
    rho_core *= 5514
    R *= 6.678e6
    M = (4/3) * np.pi * R**3 * rho_bulk
    Mw = M * mf
    Rc = (M * (1 - mf) * 3 / (4 * np.pi)) ** (1/3)
    
    res = np.array(HydroSphere(Psurf, Teq, Mw, Rc, rho_core))
    
    if savefile is not None:
        
        # columns are z, rho, alpha, Cp, dT_dz, phase, T, P, g
        np.savetxt(savefile, res.T)
    else:
        return res
    
# trappist-1 G params

rhoc = 1.5 # core density in earth units
mf = 0.1 # water / core
R = 1.129 # in earth units
Ps = 0.1 # MPa 
    
def f(c):
    i, D, T = c
    print(i, D, T)
    savefile = 'out{0}.dat'.format(i)
    compute_structure(D, T, rhoc, mf, R, Ps, savefile=savefile)
    
D = np.linspace(0.5, 1.5, 5)
T = np.linspace(150, 250, 5)
D, T = np.meshgrid(D, T)
coords = [(i, D, T) for i, (D, T) in enumerate(zip(D.flatten(), T.flatten()))]

if __name__ == '__main__':
    Pool(2).map(f, coords)

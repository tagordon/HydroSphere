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
    Mw = M * fwater
    Rc = (M * (1 - fwater) * 3 / (4 * np.pi * rho_core)) ** (1/3)
    
    res = np.array(HydroSphere(Psurf, Teq, Mw, Rc, rho_core))
    
    if savefile is not None:
        
        # columns are z, rho, alpha, Cp, dT_dz, phase, T, P, g
        np.savetxt(savefile, res.T)
    else:
        return res
    
# trappist-1 G params

rhoc = 1.0 # core density in earth units
R = 1.129 # in earth units
Ps = 0.1 # MPa 
# temps b, c, d, e, f, g: 400, 341, 288, 251, 219, 198, 168 from Agol 2020
T = 198
    
def f(c):
    i, D, fw = c
    #print(i, D, T)
    savefile = 'out{0}.dat'.format(i)
    compute_structure(D, T, rhoc, fw, R, Ps, savefile=savefile)
    
D = np.linspace(0.5, 1.0, 20)
fw = np.linspace(0.001, 0.01, 20)
D, fw = np.meshgrid(D, fw)
coords = [(i, D, fw) for i, (D, fw) in enumerate(zip(D.flatten(), fw.flatten()))]
np.savetxt('key.dat', coords)

if __name__ == '__main__':
    Pool(28).map(f, coords)

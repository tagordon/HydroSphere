import numpy as np
import sys
sys.path.append('/gscratch/astro/tagordon/hydrosphere')
from ReverseHydro import T0, melting_curve
from multiprocessing import Pool, cpu_count

outpath = sys.argv[1]
N = int(sys.argv[2])

M = float(sys.argv[3])
Rhoc = float(sys.argv[4])

M *= 5.972e24
Rhoc *= 5514

P = np.linspace(0, 3000, N)
T = melting_curve(P)

PT = [(P, T) for P, T in zip(P, T)]

def f(P, T):
    print(P)
    Tx, Mx, _, _ = T0(0.1, M, Rhoc, P=P, T=T, dz=5e-5)
    return Tx, Mx / (Mx + M)
    
if __name__ == '__main__':
    with Pool(cpu_count()) as pool:
        result = pool.starmap(f, PT)    

np.savetxt(outpath, result)

import numpy as np
import os

files = os.listdir('output')
i, rho, f = np.loadtxt('key.dat').T
dx, dy = len(np.unique(rho)), len(np.unique(f))

depth = np.zeros_like(rho)
contact_phase = np.zeros_like(rho)
water_temp_min = np.zeros_like(rho)
water_temp_max = np.zeros_like(rho)
base_pressure = np.zeros_like(rho)

for file in files:
    j = np.int64(f[3:-4])
    z, rho, alpha, Cp, dT_dz, phase, T, P, g = np.loadtxt('output/' + file).T
    
    iswater = np.where(phase == 0)
    water_layer = z[iswater]
    if len(water_layer) > 0:
        water_depth = water_layer[-1] - water_layer[0]
        water_temp_min[j] = np.min(T[iswater])
        water_temp_max[j] = np.max(T[iswater])
        depth[j] = water_depth
        
    contact_phase[j] = phase[-1]
    base_pressure[j] = P[-1]
    
depth = depth.reshape(dx, dy)
contact_phase = contact_phase.reshape(dx, dy)
water_temp_min = water_temp_min.reshape(dx, dy)
water_temp_max = water_temp_max.reshape(dx, dy)
base_pressure = base_pressure.reshape(dx, dy)

np.savetxt(
    'summary.dat', 
    np.concatenate([
        [depth], 
        [contact_phase], 
        [water_temp_min], 
        [water_temp_max],
        [base_pressure]
    ]).T
)
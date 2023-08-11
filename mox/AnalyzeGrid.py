import numpy as np
import os
import sys

outdir = sys.argv[1]
files = os.listdir(outdir)
i, x, y = np.loadtxt(outdir + '/key.dat').T
dx, dy = len(np.unique(x)), len(np.unique(y))

structure = np.zeros_like(x)
base_water_pressure = np.zeros_like(x)
depth = np.zeros_like(x)
contact_phase = np.zeros_like(x)
water_temp_min = np.zeros_like(x)
water_temp_max = np.zeros_like(x)
base_pressure = np.zeros_like(x)
total_depth = np.zeros_like(x)
has_water = np.zeros_like(x)

for file in files:
    if ('.dat' in file) & (file != 'key.dat'):
        j = np.int64(file[3:-4])
        z, rho, alpha, Cp, dT_dz, phase, T, P, g = np.loadtxt(outdir + '/' + file).T
    
        iswater = np.where(phase == 0)
        water_layer = z[iswater]
        if len(water_layer) > 0:
            water_depth = water_layer[-1] - water_layer[0]
            water_temp_min[j] = np.min(T[iswater])
            water_temp_max[j] = np.max(T[iswater])
            depth[j] = water_depth
            has_water[j] = 1.0
            base_water_pressure[j] = P[np.where(z == water_layer[-1])]
        
        contact_phase[j] = phase[-1]
        base_pressure[j] = P[-1]
        total_depth[j] = z[-1]
        
        _, idx = np.unique(phase, return_index=True)
        layers = phase[np.sort(idx)]

        if np.all(layers == [0, 7]):
        # water then ice-vii
            struc = 1
        elif np.all(layers == [1, 0, 7]):
        # ice, water, ice-vii
            struc = 2
        elif (layers[0] == 0) & (len(layers) == 1):
        # only water
            struc = 3
        elif not 0 in layers:
        # no water
            struc = 4
        elif np.all(layers == [1, 0]):
        #ice then water
            struc = 5
        elif (len(layers) > 2) & (np.all(layers[:2] == [1, 0])):
        #ice, water, hp ice
            struc = 6
        elif (layers[0] == 0) & (len(layers) > 1):
        #water then hp ice"
            struc = 7
        structure[j] = struc    

structure = structure.reshape(dx, dy)
base_water_pressure = base_water_pressure.reshape(dx, dy)
depth = depth.reshape(dx, dy)
contact_phase = contact_phase.reshape(dx, dy)
water_temp_min = water_temp_min.reshape(dx, dy)
water_temp_max = water_temp_max.reshape(dx, dy)
base_pressure = base_pressure.reshape(dx, dy)
total_depth = total_depth.reshape(dx, dy)
has_water = has_water.reshape(dx, dy)

np.savetxt(outdir + '/structure.out', structure)
np.savetxt(outdir + '/base_water_pressure.out', base_water_pressure)
np.savetxt(outdir + '/depth.out', depth)
np.savetxt(outdir + '/contact_phase.out', contact_phase)
np.savetxt(outdir + '/water_temp_min.out', water_temp_min)
np.savetxt(outdir + '/water_temp_max.out', water_temp_max)
np.savetxt(outdir + '/base_pressures.out', base_pressure)
np.savetxt(outdir + '/total_depth.out', total_depth)
np.savetxt(outdir + '/has_water.out', has_water)

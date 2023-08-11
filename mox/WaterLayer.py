import numpy as np
import os

files = os.listdir('output')
i, rho, f = np.loadtxt('key.dat').T
dx, dy = len(np.unique(rho)), len(np.unique(f))
depth = np.zeros_like(rho)
contact = np.zeros_like(rho)
for file in files:
    j = np.int64(file[3:-4])
    z, _, _, _, _, phase, _, _, _ = np.loadtxt('output/' + file).T
    water_layer = z[np.where(phase == 0.0)]
    if len(water_layer) > 0:
        water_depth = water_layer[-1] - water_layer[0]
    else:
        water_depth = 0
    in_contact = phase[-1] == 0.0
    depth[j] = water_depth
    contact[j] = in_contact
    
depth = depth.reshape(dx, dy)
contact = contact.reshape(dx, dy)

np.savetxt('depth.dat', depth)
np.savetxt('contact.dat', contact)

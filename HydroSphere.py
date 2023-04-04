import numpy as np
import seafreeze as sf
import warnings

warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

Temp_dataset = [[ 6.39162953e+02,  7.51104318e-02, -1.96063855e-04],
                [ 4.86357760e+02, -5.05830073e-01,  7.59695073e-04],
                [ 1.81045287e+02,  3.69540967e-01, -3.90757776e-04],
                [ 1.60798154e+02,  8.00941009e-01, -5.96848739e-04],
                [ 1.89382439e+02,  1.30683834e+00, -1.32469061e-03],
                [ 6.58847367e+02,  1.00020949e+00, -7.84685730e-04]]

E = [1.60, 1.25, -0.02, 0.16, 0.37, 0.65]
F = [-0.44, 0.2, 0.2, 0.2, 0.16, 0.2]

direc1 = {1:0, 2:1, 3:2, 5:3, 6:4, 7:5}

Rg = 8.314
c1 = 1.43
c2 = -0.03
# The following arrays follow 'water','Ih','II','III','V','VI'
nu0 = [0, 5e13, 1e18, 5e12, 5e14, 5e14]
E_Jmol = [0, 60e3, np.mean([98, 55])*1e3, np.mean([103, 151])*1e3, 136e3, 110e3]
direc2 = {0:0, 1:1, 2:2, 3:3, 5:4, 6:5, 7:4} # {phase:position in array} using ice V values for ice VII

# Kappa fix
# I, II, III, V, VI, VII
T_intcpt = [130, 120, 240, 246, 246, 286]   # In K
P_intcpt = [0.1, 0.24, 0.24, 0.53, 1, 2.4]  # In GPa

def _compute_conductivity(P, T, rho, alpha, grav, z, phase):
    
    boundary = np.where(np.diff(phase) >= 1)[0]
    boundary = np.concatenate([[0], boundary, [len(z)-1]])
    
    Ra = np.zeros(len(boundary))
    Conductivity = np.zeros(len(boundary))
    
    k_Tfix = np.zeros(6)
    k_Pfix = np.zeros(6)
    for i in range(6):
        k_Tfix[i] = Temp_dataset[i][0]/T_intcpt[i] + Temp_dataset[i][1] + Temp_dataset[i][2]*T_intcpt[i]
        k_Pfix[i] = np.exp(E[i] + F[i]*P_intcpt[i])
    
    for i in range(1, len(boundary)):
        upper = boundary[i]
        lower = boundary[i-1]
        if (phase[upper] == 0):
            Ra[i] = -1
            Conductivity[i] = -1
        else:
            dir1 = direc1[phase[upper]]
            dataset = Temp_dataset[dir1]
            k = dataset[0]/T[lower]+dataset[1]+dataset[2]*T[lower]+np.exp(E[dir1]+F[dir1]*P[lower]*1e-3)
            Kappa = k/rho[lower]/Cp[lower]
        
            dir2 = direc2[phase[upper]]
            A = E_Jmol[dir2]/Rg/T[lower]
            B = E_Jmol[dir2]/2/Rg/c1
            C = c2*(T[lower]-T[i])
            Tc = B*(np.sqrt(1+2/B*(T[lower]-C))-1)
            nu = nu0[dir2]*np.exp(A*(T[lower]/Tc-1))

            Ra[i] = alpha[lower]*rho[lower]*grav[lower]*(T[lower]-T[upper])*z[lower]**3/Kappa/nu
            Conductivity[i] = -k*np.abs(T[lower]-T[upper])
            
    return Ra, Conductivity

def HydroSphere(Ps, Ts, Mw, Rc, Rhoc, resolution=100, G_iter=1, M_thresh=0.01, compute_conductivity=False):
        
    Mass_core = 4./3 * np.pi * Rc**3 * Rhoc
    g_s = 6.67430e-11 * Mass_core / Rc**2 # Gravity at the Hydrosphere Mantle Boundary
    depth = (Mw / 800 / (4./3. * np.pi) + Rc**3)**(1./3.) - Rc # Depth approximation in m

    # initializing the grids
    z = np.linspace(0, depth, num=resolution)  # depth grid
    dz = np.diff(z)
    
    rho = np.empty(resolution)  # Density grid
    alpha = np.empty(resolution)  # thermal epansivity grid
    Cp = np.empty(resolution)  # Heat capacity grid
    dT_dz = np.empty(resolution)  # thermal gradient grid
    phase = np.empty(resolution)  # phase grid
    T = np.empty(resolution)  # Temperature grid
    P = np.empty(resolution)  # Presolutionsure grid
    grav = np.empty(resolution) # gravity grid
    PT = np.empty((1,), np.object)
    
    Mass_diff = 1
    
    while (Mass_diff > M_thresh):
        
        # For mass loop the factor being iterated is /depth/
        
        grav = g_s * np.ones(resolution)
    
        # Gravity conversion loop
        for k in range(G_iter): 

            # For gravity loop the factor being iterated is /grav/
            g = grav[::-1]
            PT[0] = (Ps, Ts)
            
            phase_s = 7 if Ps > 2200 else sf.whichphase(PT)[0]
            out = sf.seafreeze(PT,sf.phasenum2phase[phase_s]) 

            T[0] = Ts
            P[0] = Ps
            rho[0] = out.rho
            alpha[0] = out.alpha
            Cp[0] = out.Cp
            dT_dz[0] = out.alpha * g[0] * Ts / out.Cp
            phase[0] = phase_s
            
            for i in range(1, z.size):  # Integration with depth
                T[i] = T[i - 1] + dT_dz[i - 1] * dz[i - 1]
                P[i] = P[i - 1] + rho[i - 1] * g[i - 1] * dz[i - 1] * 1e-6
                PT[0] = (P[i], T[i])
                
                phase[i] = 7 if P[i] > 2200 else sf.whichphase(PT)
                out = sf.seafreeze(PT,sf.phasenum2phase[phase[i]])

                rho[i] = out.rho
                alpha[i] = out.alpha
                Cp[i] = out.Cp
                dT_dz[i] = out.alpha * g[i] * T[i] / out.Cp
                
            M_L = rho[1:] * 4/3 * np.pi * ((Rc + z[:-1] + (depth / resolution)) ** 3 - (Rc + z[:-1])**3)
            M_L = np.insert(M_L, 0, 0)
            Mass_Shells = np.cumsum(M_L[::-1])
            grav = 6.67430e-11 * (Mass_core + Mass_Shells) / (Rc + z)**2

        Mass_WL = np.sum(M_L)
        Mass_diff = Mw - Mass_WL
        depth_diff = (Mass_diff / (np.mean(rho) * 1.8) / (4/3 * np.pi) + Rc**3) ** (1/3) - Rc
        depth += depth_diff
        
        Mass_diff = np.abs(Mass_diff / Mw)
    
    if compute_conductivity:
        
        Ra, Conductivity = _compute_conductivity(P, T, rho, alpha, grav, z, phase)    
        return z, rho, alpha, Cp, dT_dz, phase, T, P, grav, Ra, Conductivity
    
    else:
        
        return z, rho, alpha, Cp, dT_dz, phase, T, P, grav
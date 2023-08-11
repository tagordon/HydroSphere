import numpy as np
import seafreeze as sf
import warnings

warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

def melting_curve(P):
    
    a0, b0, c0 = 2.61817424e-01, 4.98751174e+01, 2.92900510e+02
    a1, b1, c1 = 1.10248081,  32.12066695, 236.20399626
    a2, b2, c2 = 1.12198087,  90.57968215, 201.11686541
    a3, b3, c3 = 1.17655548, 163.58744868, 147.48886742
    
    model = lambda P, a, b, c: ((P * 1e-3 - a) + 1)**(1/3) * b + c
    model0 = lambda P, a, b, c: -((P * 1e-3 - a) + 1)**3 * b + c
    
    trip1 = (208.566, 251.165) # ice-1h, III, liquid
    trip2 = (350.1, 256.164) # ice-III, V, liquid
    trip3 = (632.4, 273.31) # ice-V, VI, liquid
    
    if (P <= trip1[0] - 4.8):
        return model0(P, a0, b0, c0)
    if (P <= trip2[0]):
        return model(P, a1, b1, c1)
    if (P <= trip3[0]):
        return model(P, a2, b2, c2)
    if (P <= 2146.5):
        return model(P, a3, b3, c3)
    if P > 2146.5:
        return ((P * 1e-3 - 2.17) / 1.253 + 1)**(1/3) * 354.8
    
melting_curve = np.vectorize(melting_curve)

# X_fe between 0.15 and 0.8, and M in Earth-masses, returns radius in Earth-radii
def R(M, X_fe):
    return ((7030 - 1840 * X_fe) * M ** 0.282) / 6371

def T0(Ps, Mc, Rhoc, P=208.566, T=251.165, dz=0.001):
    
    PT = np.empty((1,), np.object)
    
    Rc = ((3 * Mc) / (4 * np.pi * Rhoc)) ** (1./3.)
    dz *= Rc

    z = Rc
    M = Mc
    g = 6.67430e-11 * M / (Rc * Rc)
    
    Pz = []
    Tz = []
    
    while P > Ps:
        
        PT[0] = (P, T)
        
        if P > 2200:
            Tm = ((Ps * 1e-3 - 2.17) / 1.253 + 1)**(1/3) * 354.8
            if T < Tm:
                phase_s = 7
            else:
                phase_s = 0
        else:
            phase_s = sf.whichphase(PT)[0]
                
        out = sf.seafreeze(PT,sf.phasenum2phase[phase_s]) 
                
        M += (4./3.) * np.pi * ((z + dz)**3 - z**3) * out.rho[0]
        T -= out.alpha[0] * g * T * dz / out.Cp[0]
        P -= out.rho[0] * g * dz * 1e-6
        
        z += dz
        g = 6.67430e-11 * M / z**2
        
        Pz.append(P)
        Tz.append(T)

    return T, M - Mc, Pz, Tz

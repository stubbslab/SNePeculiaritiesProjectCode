import numpy as np
import math
from scipy.integrate import odeint
from scipy.optimize import fsolve
from scipy.integrate import quad
from scipy.interpolate import interp1d 
from CosmologicalParameterArchive import CosmologicalParameterArchive 
import scipy
from scipy.optimize import minimize_scalar
from scipy.optimize import minimize
from cantrips import insertListElements
from computeMuForDifferentCosmologies import computeMuForCosmology 
#import sympy as sp 

def BoatFishSystem(state, t):
    fish, boat = state
    d_fish = fish * (2 - boat - fish)
    d_boat = -boat * (1 - 1.5 * fish)
    return [d_fish, d_boat]

def getBoatFishState(init_state, t):
    return odeint(BoatFishSystem, init_state, t)

def getV(phi, C, A, alpha, nu, phi_0):
    V = C * ((phi / phi_0) ** (- alpha)) * (1.0 - A * np.sin(nu * phi) )
    #print 'V = ' + str(V) 
    return V

def getDV(phi, C, A, alpha, nu, phi_0):
    
    dV =  ( C * (- alpha) / phi_0 * (phi / phi_0) ** (- alpha - 1.0) * (1.0 - A * np.sin(nu * phi) )
             + C * (phi / phi_0) ** (- alpha) * (- A * nu * np.cos(nu * phi) ) )
    #print 'dV = ' + str(dV)
    return dV 

def getH(zs, phi, theta, C, A, alpha, nu, phi_0, Omega_m0, flip_z = 0):
    #if type(zs) in [int, float] and (zs * (flip_z - 0.5) * (-2.0)) < 0.0:
    #    print 'Warning: tried to use invalid z value of ' + str(zs) + ' when flip_z = ' + str(flip_z) 
    if flip_z:
        H = np.sqrt(Omega_m0 * (1.0 - zs) ** 3.0 + (8.0 * math.pi / 3.0) * (theta ** 2.0 / 2.0 + getV(phi, C, A, alpha, nu, phi_0) ) )
    else: 
        H = np.sqrt(Omega_m0 * (1.0 + zs) ** 3.0 + (8.0 * math.pi / 3.0) * (theta ** 2.0 / 2.0 + getV(phi, C, A, alpha, nu, phi_0) ) )
    #print 'H = ' + str(H)
    return H 

def getDPhi(zs, phi, theta, C, A, alpha, nu, phi_0, Omega_m0, flip_z = 0):
    #print 'flip_z = ' + str(flip_z) 
    if flip_z: 
        dPhi = -1.0 * theta / ( getH(zs, phi, theta, C, A, alpha, nu, phi_0, Omega_m0, flip_z = flip_z) * (1.0 - zs) )
    else:
        dPhi = -1.0 * theta / ( getH(zs, phi, theta, C, A, alpha, nu, phi_0, Omega_m0, flip_z = flip_z) * (1.0 + zs) )
    #print 'dPhi = ' + str(dPhi)
    if flip_z:
        return -dPhi
    else:
        return dPhi 

def getDTheta(zs, phi, theta, C, A, alpha, nu, phi_0, Omega_m0, flip_z = 0):
    #print 'flip_z = ' + str(flip_z) 
    if flip_z:
        dTheta = 3.0 * theta / (1.0 - zs) + 1.0 / (getH(zs, phi, theta, C, A, alpha, nu, phi_0, Omega_m0, flip_z = flip_z) * (1.0 - zs) ) * getDV(phi, C, A, alpha, nu, phi_0)
    else:
        dTheta = 3.0 * theta / (1.0 + zs) + 1.0 / (getH(zs, phi, theta, C, A, alpha, nu, phi_0, Omega_m0, flip_z = flip_z) * (1.0 + zs) ) * getDV(phi, C, A, alpha, nu, phi_0)
    #print 'dTheta = ' + str(dTheta) 
    if flip_z:
        return -dTheta
    else:
        return dTheta

def MonodromeDEDifferentialEqns(state, zs, C, A, alpha, nu, phi_0, Omega_m0, flip_z):
    phi = state[0]
    theta = state[1]
    #print 'phi = ' + str(phi) 
    #V = getV (phi, C, A, alpha, nu, phi_0) 
    #H = getH (zs, V, theta, A, Omega_m0, gamma)
    #d_phi = - theta / ( H * (1.0 + zs) )
    #d_theta = 3.0 * theta / (1.0 + zs) + 1 / ( H * (1.0 + zs) ) * getDV(phi, C, A, alpha, nu, phi_0)
    d_phi = getDPhi(zs, phi, theta, C, A, alpha, nu, phi_0, Omega_m0, flip_z = flip_z)
    d_theta = getDTheta(zs, phi, theta, C, A, alpha, nu, phi_0, Omega_m0, flip_z = flip_z)
    return [d_phi, d_theta]

def getMonodromeDEState(init_state, zs, C, A, alpha, nu, phi_0, Omega_m0, flip_z = 0): 
    state = odeint(MonodromeDEDifferentialEqns, init_state, zs, args = (C, A, alpha, nu, phi_0, Omega_m0, flip_z), hmax = 1.0,mxstep=5000000 )
    return state

def equationForRedshift0(p, C, A, alpha, nu, phi_0, Omega_m0):
    init_phi = p[0]
    #calculated_values = (getH(0.0, init_phi, init_theta, C, A, alpha, nu, phi_0, Omega_m0) - 1.0,
    #        1.0 / 2.0 * (Omega_m0 * 3.0 + 8.0 * math.pi / 3.0 * (init_theta * getDTheta(0.0, init_phi, init_theta, C, A, alpha, nu, phi_0, Omega_m0) * 2.0
    #                                                             + getDV(init_phi, C, A, alpha, nu, phi_0) * getDPhi(0.0, init_phi, init_theta, C, A, alpha, nu, phi_0, Omega_m0) )
    #                    ) - 3.0 / 2.0 * Omega_m0
    #        )
    calculated_value = (1 - Omega_m0 - 8.0 * math.pi / 3.0 * C * (init_phi / phi_0) ** (-alpha) * (1.0 - A * np.sin(nu * init_phi)))
    #print 'calculated_value = ' + str(calculated_value)
    return calculated_value
                                                                 
def getHFromMonodromeDE(zs, A, alpha, nu, phi_0, Omega_m0, forced_init_phi = None, init_z = 1000.0, flip_z = 1): #C, init_H = 1.0 taken out 
    #We use the ansatz that phi(z) = phi_0 * ((1.0 + z)/(1.0 + z_0)) ^ power for some very early z
    power = -3.0 / (2.0 + alpha)
    init_phi = phi_0 * (1.0 + init_z) ** power
    init_theta = - 1.0 * np.sqrt(Omega_m0) * phi_0 * power * (1.0 + init_z) ** (power + 3.0 / 2.0)  
    C = phi_0 ** 2.0 * Omega_m0 / alpha * (power ** 2.0 - 3.0 / 2.0 * power)
    init_H = np.sqrt(Omega_m0 * (1.0 + init_z) ** 3.0)

    #init_phi_funct = lambda init_phi: 1 - Omega_m0 - 8.0 * math.pi / 3.0 * C * (init_phi / phi_0) ** (-alpha) * (1.0 - A * np.sin(nu * init_phi) )
    #if alpha == 0.0 and (nu == 0.0 or A == 0.0): 
    #    if C >= (1.0 - Omega_m0): C = (1.0 - Omega_m0)
    #    init_theta = math.sqrt(2.0) * math.sqrt(1.0 - C - Omega_m0)
    #    print 'init_theta = ' + str(init_theta) 
    #elif alpha == 0.0:
    #    init_guess = [np.arcsin((((1.0 - Omega_m0) * (-3.0) / (8.0 * math.pi) * 1.0 / C + 1.0)/ A) % (math.pi / 2.0)) / nu ]
    #else:
    #    init_V = [(3.0 / (8.0 * math.pi * C) * (1.0 - Omega_m0)) ** (- 1.0 / alpha) * phi_0]
    #print 'init_guess = ' + str(init_guess) 
    #init_phi = scipy.optimize.broyden1(init_phi_funct, init_guess)
    #if not(forced_init_phi is None): 
    #    init_phi = forced_init_phi
    #else:
    #    init_V = init_H ** 2.0 - Omega_m0 * (1.0 + init_z) ** 3.0 - init_theta ** 2.0 / 2.0
    #    init_phi = (3.0 / (8.0 * math.pi * C) * (init_V)) ** (- 1.0 / alpha) * phi_0
    #fsolve(equationForRedshift0, (0.5), args = (C, A, alpha, nu, phi_0, Omega_m0) )
    
    init_state = [init_phi, init_theta] #calculate initial state by constraints on H and d H /d z today
    #print 'init_state = ' + str(init_state)
    
    if flip_z:
        if not(np.shape(zs) is ()): 
            zs = -1.0 * np.flip(zs, 0)
        else:
            zs = -1.0 * zs 
        init_z = -1.0 * init_z

    if type(zs) is list:
        zs = np.array([init_z] + zs)
    elif type(zs) is float or type(zs) is int or type(zs) is np.float64 :
        zs = np.array([init_z] + [zs])
    else:
        zs = np.array([init_z] + zs.tolist())
        
    state = getMonodromeDEState(init_state, zs, C, A, alpha, nu, phi_0, Omega_m0, flip_z = flip_z)
    phi = state[:,0]
    theta = state[:,1] 
    
    V = getV (phi, C, A, alpha, nu, phi_0)
    dV = getDV(phi, C, A, alpha, nu, phi_0)
    dTheta = getDTheta(zs, phi, theta, C, A, alpha, nu, phi_0, Omega_m0, flip_z = flip_z)
    H = getH(zs, phi, theta, C, A, alpha, nu, phi_0, Omega_m0, flip_z = flip_z)
    w = ( ((H * (1 + zs) * theta) ** 2.0) / 2.0 - V ) / ( ((H * (1 + zs) * theta) ** 2.0) / 2.0 + V)# Pressure / Energy Density
    rho = ( (H ** 2.0) * (1 + zs) ** 2.0 * theta ** 2.0 + V)
    #print 'The followng array should be 0: '
    #print getDTheta(zs, phi, theta, C, A, alpha, nu, phi_0, Omega_m0) * (- H * (1.0 + zs)) + 3.0 * H * theta + dV 
    #print 'getH(0.0, init_phi, init_theta, C, A, alpha, nu, phi_0, Omega_m0) = '
    #print (getH(0.0, init_phi, init_theta, C, A, alpha, nu, phi_0, Omega_m0))

    #Truncate the first value returned, since that was the initial z, which was not part of the input z array
    if flip_z:
        #print 'len(phi[1:]) = ' + str(len(phi[1:]))
        return [np.flip(phi[1:],0), np.flip(theta[1:],0), np.flip(H[1:],0), np.flip(V[1:],0), np.flip(w[1:],0), np.flip(rho[1:],0)]
    else:
        return [phi[1:], theta[1:], H[1:], V[1:], w[1:], rho[1:]]

#For some given parameters, returns parameters according to the differential eqns for Monodrome DE,
# with the other parameters scaled such that the average difference between the calculated H and the expected H
# is minimized.
# These parameters can be fed into the appropriate function below to get Hs/dLs/mus.  
def getParamsForFlattenedHForMonodromeDE(zs, fixed_params, minimization_params_dict, Omega_m0 = 0.3, bounds_dict = {}, forced_init_phi = None, init_z = 1000.0, flip_z = 1, tol = 1e1-5, use_canonical_cosmology = 1, n_zs_to_compute_H = 1000):
    if use_canonical_cosmology:
        cosmo_archive = CosmologicalParameterArchive()
        Omega_m0 = cosmo_archive.getOmegaM()[0]
    #Deal with fact that I may pass dictionary or fixed params with extraneous parameters not considered here
    trimmed_minimization_params_dict = {}
    trimmed_bounds_dict = {}
    for key in minimization_params_dict.keys():
        if key <= 3:
            trimmed_minimization_params_dict[key] = minimization_params_dict[key]
            if len(bounds_dict) != 0: trimmed_bounds_dict[key] = bounds_dict[key]
    bounds_dict = trimmed_bounds_dict
    minimization_params_dict = trimmed_minimization_params_dict

    n_allowed_params = 4
    fixed_params = fixed_params[0:n_allowed_params - len(minimization_params_dict)]

    if len(bounds_dict) != 0:
        bounds = tuple([bounds_dict[key] for key in sorted(bounds_dict.keys()) ])
    else:
        bounds = ()
        
    max_z = max(zs)
    min_z = min(zs) 
    zs_to_compute_H = np.linspace(0.0, max_z, n_zs_to_compute_H)
    
    minimization_params_keys = list(sorted(minimization_params_dict.keys()))
    getMinimizationParamsDict = lambda minimization_params: dict(zip(minimization_params_keys, minimization_params))
    canonical_H = np.sqrt(Omega_m0 * (1.0 + np.array(zs_to_compute_H)) ** 3.0 + (1.0 - Omega_m0))
    functToMinimize = lambda minimization_params: abs(np.mean(np.array(getHFromMonodromeDE(zs_to_compute_H, *(insertListElements(fixed_params, getMinimizationParamsDict(minimization_params) ) + [Omega_m0]), init_z = init_z)[2]) / canonical_H - 1.0))

    init_minimization_params = [minimization_params_dict[key] for key in sorted(minimization_params_dict.keys()) ]
    fixed_params_results = minimize(functToMinimize, init_minimization_params, bounds = bounds)

    print 'minimized average difference between measured H(z) and canonical H(z) is ' + str(functToMinimize(fixed_params_results['x'])) + ' achieved at ' + str(fixed_params_results['x'])
    return dict(zip(minimization_params_keys, fixed_params_results['x']) )
    #return fixed_params_results

#Returns H according to the DEs for Monodrome differential eqns, with phi_0 scaled such that H(z = 0) = H0
# (this is the technique employed in the papar)
def getNormalizedHForMonodromeDE(zs, A, alpha, nu, Omega_m0, forced_init_phi = None, init_z = 1000.0, flip_z = 1, min_z = 0.00001):
    minim_phi_0 = minimize_scalar(lambda phi_0: (getHFromMonodromeDE([min_z], A, alpha, nu, phi_0, Omega_m0, forced_init_phi = forced_init_phi, init_z = init_z, flip_z = flip_z)[2][0] - 1.0 ) ** 2.0, bounds = [0.001, 10.0], method = 'bounded', options  = {'xatol':1e-07})['x']
    print 'minim_phi_0 = ' + str(minim_phi_0) 
    return getHFromMonodromeDE(zs, A, alpha, nu, minim_phi_0, Omega_m0, forced_init_phi = forced_init_phi, init_z = init_z, flip_z = flip_z)


def getDLForMonodromeDE(zs, A, alpha, nu, Omega_m0 = 0.3, phi_0 = 0.1, forced_init_phi = None, init_z = 1000.0, flip_z = 0, n_zs_to_compute_H = 1000, H0 = 70.0, use_canon_params = 1, normalize_H_to_H0 = 1):
    #Would be extremelly expensive to redo every integral.  So define an interpolator to compute over
    cosmo_archive = CosmologicalParameterArchive()
    if use_canon_params:
        H0 = cosmo_archive.getH0()[0]
        Omega_m0 = cosmo_archive.getOmegaM()[0]
    c = cosmo_archive.getc()
    
    if type(zs) in [float, int]:
        zs = [zs]
    max_z = max(zs)
    min_z = min(zs) 
    zs_to_compute_H = np.arange(0.0, max_z + max_z / n_zs_to_compute_H * 1.5, max_z / n_zs_to_compute_H)
    if normalize_H_to_H0:
        computed_Hs = getNormalizedHForMonodromeDE(zs_to_compute_H, A, alpha, nu, Omega_m0, forced_init_phi = forced_init_phi, init_z = init_z, flip_z = flip_z)[2]
    else:
        computed_Hs = getHFromMonodromeDE(zs_to_compute_H, A, alpha, nu, phi_0, Omega_m0, forced_init_phi = forced_init_phi, init_z = init_z, flip_z = flip_z)[2]
    H_interp = interp1d(zs_to_compute_H, computed_Hs, kind = 'cubic') 

    scaling = c / H0 # units of Mpc
    dLs = [scaling * (1.0 + z) * quad(lambda z: (1.0 / H_interp(z)) , 0.0, z)[0] for z in zs ]
    return dLs 

def getMuForMonodromeDE(zs, A, alpha, nu, Omega_m0 = 0.3, phi_0 = 0.1, forced_init_phi = None, init_z = 1000.0, flip_z = 1, n_zs_to_compute_H = 1000, normalize_H_to_H0 = 1):

    mus = (25.0 + 5.0
           * np.log10(getDLForMonodromeDE(zs, A, alpha, nu, Omega_m0,
                                          phi_0 = phi_0, 
                                          forced_init_phi = forced_init_phi, init_z = init_z, flip_z = flip_z, n_zs_to_compute_H = n_zs_to_compute_H, normalize_H_to_H0 = normalize_H_to_H0))
           )
    return mus

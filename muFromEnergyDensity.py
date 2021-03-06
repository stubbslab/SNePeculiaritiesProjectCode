import numpy as np
import math
from scipy.integrate import quad
from CosmologicalParameterArchive import CosmologicalParameterArchive

def getIntegrableFunction(function, params_to_fix):
    return lambda x: function(x, *params_to_fix)

def muFromEnergyDensity(energy_density_function):
    cosmo_archive =  CosmologicalParameterArchive()
    OmM = cosmo_archive.getOmegaM()[0]
    OmL = cosmo_archive.getOmegaLambda()[0]
    Om0 = cosmo_archive.getOmega0()[0]
    H0 = cosmo_archive.getH0()[0]
    c = cosmo_archive.getc()
    scaling = c / H0
    canonicalH = lambda z: np.sqrt(OmM * (1.0 + z) ** 3.0 + OmL)
    numer = lambda z, *fit_params: ( quad ((lambda x: scaling * (getIntegrableFunction(energy_density_function, fit_params)(x)/canonicalH(x))) , 0, z) 
                                             )
    denom = lambda z: (quad (lambda x: 1/canonicalH(x)**3.0 , 0, z))
    
    return lambda z, *fit_params: ( (numer(z, *fit_params)[0] / denom(z)[0]) if (type(z) is float or type(z) is int) else np.array( [numer(z_elem, *fit_params)[0] / denom(z_elem)[0] for z_elem in z]) )

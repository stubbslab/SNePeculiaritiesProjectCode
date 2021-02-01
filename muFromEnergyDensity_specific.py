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
    numer = lambda z, a, b, c, d: ( quad ((lambda x: scaling * (energy_density_function(x, a, b, c, d)/canonicalH(x))) , 0, z) 
                                             )
    denom = lambda z: (quad (lambda x: 1/canonicalH(x)**3.0 , 0, z))
    return lambda z, a, b, c, d: numer(z, a, b, c, d)[0] / denom(z)[0]

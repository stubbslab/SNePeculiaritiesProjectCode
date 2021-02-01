import math
import numpy as np
import CosmologicalParameterArchive as cpa
import scipy.integrate as integrate

def lum_dist(z, OmM = 0.3, OmL = 0.7, H0 = 70.0, Om0 = 1.0, use_canon_params = 1):
    #dp = c * integr.quad(1 / a(t) )
    cosmo_archive = cpa.CosmologicalParameterArchive(params_source = 'pantheon')
    if use_canon_params:
        H0 = cosmo_archive.getH0()[0]
        OmM = cosmo_archive.getOmegaM()[0]
        OmL = cosmo_archive.getOmegaLambda()[0]
        Om0 = cosmo_archive.getOmega0()[0]
    c = cosmo_archive.getc()
    scaling = c / H0 #units of MpC
    dl = scaling * (1.0+z) * integrate.quad(lambda x: 1.0 / math.sqrt( OmM * ((1+x) ** (3.0)) + OmL + (1.0 - Om0) * ((1+x) ** (2.0)) ), 0.0, z)[0]
    #print 'for z = ' + str(z) + ', dl = ' + str(dl)
    return dl

def computeMuForCosmology(z, OmM = 0.3, OmL = 0.7, H0 = 70.0, Om0 = 1.0, OmR = 0.0, use_canon_params = 1):
    mu = 5.0 * np.log10(lum_dist(z, OmM = OmM, OmL = OmL, H0 = H0, Om0 = Om0, use_canon_params = use_canon_params)) + 25.0 # distance modulus
    return mu

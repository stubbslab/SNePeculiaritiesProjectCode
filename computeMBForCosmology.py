import math
import numpy as np
from CosmologicalParameterArchive import CosmologicalParameterArchive
import scipy.integrate as integrate

def lum_dist(z):
    #dp = c * integr.quad(1 / a(t) ) 
    cosmo_archive = CosmologicalParameterArchive()
    H0 = cosmo_archive.getH0()[0]
    OmM = cosmo_archive.getOmegaM()[0]
    OmL = cosmo_archive.getOmegaLambda()[0]
    Om0 = cosmo_archive.getOmega0()[0]
    c = cosmo_archive.getc()
    #print 'Om0 = ' + str(Om0)
    #print type(Om0) 
    #print 'OmL = ' + str(OmL)
    #print type(OmL) 
    #print 'OmM = ' + str(OmM)
    #print type(OmM) 
    #print 'z = ' + str(z)
    #print type(z) 
    scaling = c / H0 #unite of MpC 
    dl = scaling * (1.0+z) * integrate.quad(lambda a: 1.0 / math.sqrt(OmM * a **3.0 + OmL * a ** 6.0 + (1 - Om0) * a ** 4.0 ), 1.0 / (1.0 + z), 1.0)[0]
    #print 'dl = ' + str(dl) 
    return dl 

def computeMBForCosmology(z):
    mB = 5 * np.log(lum_dist(z)) + 25.0# distance modulus 
    return mB 

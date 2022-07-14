#Computes cosmological parameters (z, a, t) from a a conanical cosmology
# Does allow for some nonstandard computations as well.

import math
import numpy as np
from CosmologicalParameterArchive import CosmologicalParameterArchive
import scipy.integrate as integrate

class CosmologicalCalculator:

    def getTsFromZs(self, zs):
        #to = (1.0 / self.H0_invSec) * integrate.quad(lambda x: 1.0 / (1.0 + x) * 1.0 / (self.H_funct(x)), 0, np.inf)[0]
        ts = [(1.0 / self.H0_invSec) * integrate.quad(lambda x: 1.0 / (1.0 + x) * 1.0 / (self.H_funct(x)), 0, z)[0] for z in zs]
        if self.verbose: print ('Cosmological calculator computed ts from zs. ')
        return ts

    def getTausFromZs(self, zs):
        #ts = [self.t0 - (1.0 / self.H0_invSec) * integrate.quad(lambda x: 1.0 / (self.H_funct(x)),0, z)[0] for z in zs]
        #tau0 = (1.0 / self.H0_invSec) * integrate.quad(lambda x: 1.0 / (self.H_funct(x)), 0, np.inf)[0]
        taus = [(1.0 / self.H0_invSec) * integrate.quad(lambda x: 1.0 / (self.H_funct(x)), 0, z)[0] for z in zs]
        if self.verbose: print ('Cosmological calculator computed taus from zs. ')
        return taus

    def getTsFromAs(self, a_params):
        zs = self.getZsFromAs(a_params)
        ts = getTsFromZs(zs)
        return ts

    def getTausFromAs(self, a_params):
        zs = self.getZsFromAs(a_params)
        ts = getTausFromZs(zs)
        return taus

    def getZsFromAs(self, a_params):
        zs = [1.0 / a - 1.0 for a in a_params]
        return zs

    def getAsFromZs(self, zs):
        a_params = [1.0/ (1.0 + z) for z in zs]
        return a_params

    #Gives rho_w(z) for some eos parameter, w, that can be a function of redshift
    def getRhoOfZFunctFromWOfZFunct(self, wOfZFunct, rho_today = 1):
        exponent = lambda zs: ( 3.0 * np.array(integrate.quad(lambda x: (1.0 + wOfZFunct(x)) / (1.0 + x) ,0.0, zs)[0])
                               if np.shape(np.array(zs)) is ()
                               else 3.0 * np.array([integrate.quad(lambda x: (1.0 + wOfZFunct(x)) / (1.0 + x) ,0.0, z)[0] for z in np.array(zs)])
                              )
        rhoOfZFunct = lambda zs: rho_today * np.exp(exponent(zs))
        return rhoOfZFunct

    def getRhoOfZValsFromWOfZFunct(self, wOfZFunct, rho_today = 1):
        rhoOfZFunct = self.getRhoOfZFunctFromWOfZFunct(wOfZFunct, rho_today = rho_today)
        return rhoOfZFunct(self.zs)

    def __init__(self, zs = None, ts = None, a_params = None, use_canon_params = 1, H_funct = 'standard', H0 = 70.0, OmM = 0.3, OmL = 0.7, Om0 = 1.0, OmR = 0.00001 * 8.413, verbose = 1):
        cosmo_archive = CosmologicalParameterArchive()
        self.verbose = verbose
        if use_canon_params:
            self.H0 = cosmo_archive.getH0()[0]
            self.H0_invSec = cosmo_archive.getH0_invSec()[0]
            self.OmM = cosmo_archive.getOmegaM()[0]
            self.OmL = cosmo_archive.getOmegaLambda()[0]
            self.Om0 = cosmo_archive.getOmega0()[0]
            self.t0 = cosmo_archive.gett0()[0]
        self.c = cosmo_archive.getc()
        self.OmR = OmR
        if H_funct in ['standard', 'usual', None]:
            self.H_funct = lambda z: np.sqrt(self.OmL + self.OmM * (1.0 + z) ** 3.0 + self.OmR * (1.0 + z) ** 4.0 + (1.0 - self.Om0) * (1.0 + z) ** 2.0 )
        else:
            self.H_funct = H_funct
        if not(a_params is None):
            self.a_params = a_params
            self.zs = self.getZsFromAs(a_params)
            self.ts = self.getTsFromAs(a_params)
            self.taus = self.getTausFromAs(a_params)
        elif not(zs is None):
            self.zs = zs
            self.a_params = self.getAsFromZs(zs)
            self.ts = self.getTsFromZs(zs)
            self.taus = self.getTausFromZs(zs)
        elif not(ts is None):
            print ('Not yet sure how to compute t from z. ')
            return 0
        else:
            print ('Must specify one of zs, ts, a_params.  got all of them as None')

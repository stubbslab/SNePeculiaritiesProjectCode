#Solves the set of coupled ODEs acquired when w_lambda is allowed to be an arbitrary function of time.
# The coupled ODEs are of a(t) and rho_Lambda(t)
# H(t) can then be calculated numerically from this result. 


import numpy as np
import math
from scipy.integrate import odeint
from scipy.optimize import fsolve
from scipy.integrate import quad
from scipy.interpolate import interp1d 
from CosmologicalParameterArchive import CosmologicalParameterArchive
from AstronomicalParameterArchive import AstronomicalParameterArchive 
import scipy
from scipy.optimize import minimize_scalar
from scipy.optimize import minimize
from cantrips import insertListElements
from computeMuForDifferentCosmologies import computeMuForCosmology 
#import sympy as sp

class ResidualMuCalculatorForArbitraryXofT: 

    def BoatFishSystem(state, t):
        fish, boat = state
        d_fish = fish * (2 - boat - fish)
        d_boat = -boat * (1 - 1.5 * fish)
        return [d_fish, d_boat]

    def getBoatFishState(init_state, t):
        return odeint(BoatFishSystem, init_state, t)

    def getH(self, t, tau, z):
        #print 'a = ' + str(a)
        #print 'OmL = ' + str(OmL)
        X = self.getX(t, tau, z)
        #print 'type(X) = ' + str(type(X)) 
        #print 'X = ' + str(X)
        #print 'type(X) in [list]) = ' + str(type(X) in [list]) 
        if (type(X) in [list]): X = np.array(X) 
        #print 'type(z) = ' + str(type(z)) 
        #print 'z = ' + str(z)
        #print '[self.OmM0, self.OmR0, self.OmL0] = ' + str([self.OmM0, self.OmR0, self.OmL0]) 
        return np.sqrt((1.0 + z) ** 3.0 * self.OmM0 + (1.0 + z) ** 4.0 * self.OmR0 + self.OmL0 * X)
    
    def getX(self, t, tau, z):
        #print 'tau = ' + str(tau)
        #print 'type(tau) = ' + str(type(tau))
        #print "type(tau) in ['list',np.ndarray] = " + str(type(tau) in ['list',np.ndarray])
        if not(type(tau) in ['list',np.ndarray]): tau = float(tau) 
        return self.XOfFunction(t, tau, z)

    def tOmegaLumDistSystem(self, state, z):
        t, tau, dL = state
        self.states = self.states + [state]
        #t, OmL = state
        #print '[z, t, OmL] = ' + str([z,t,OmL])
        H = self.getH(t, tau, z)
        #print 'H = ' + str(H)
        d_t = -1.0 * (-1.0 / (1.0 + z) * 1.0 / H)
        d_tau = 1.0 / H 
        d_dL = (dL / (1.0 + z)) + ((1.0 + z) / H)
        derivs = [d_t, d_tau, d_dL]
        self.derivs = self.derivs + [derivs]
        return derivs 
        #print '[d_t, d_OmL] = ' + str([d_t, d_OmL]) 
        #return [d_t, d_OmL]

    def gettOmegaLumDistSystem(self, init_state, z):
        return odeint(self.tOmegaLumDistSystem, init_state, z)

    def runCalculationWithGivenzs(self, zs):
        #print 'zs = ' + str(zs)
        if zs[0] > 0.0: zs = [0.0] + zs #Initial values are known at z = 0.0 
        self.zs = zs 
        self.states = self.gettOmegaLumDistSystem([self.t0, self.tau0, self.dL0], self.zs)
        #self.states = self.gettOmegaLumDistSystem([self.t0, self.OmL0], self.zs)
        self.tH0_of_z = [state[0] for state in self.states]
        self.tauH0_of_z = [state[1] for state in self.states]
        self.dL_of_z = [state[2] for state in self.states]
        self.X_of_z = self.getX(np.array(self.tH0_of_z), np.array(self.tauH0_of_z), np.array(self.zs))
        self.H_of_z = self.getH(np.array(self.tH0_of_z), np.array(self.tauH0_of_z), np.array(self.zs)) 
        self.a_of_z = 1.0 / (1.0 + np.array(self.zs))

    def getUnitfullTs(self):
        #print 'Units on ts are ' + self.H0_units
        return [tH0 / self.H0 for tH0 in self.tH0_of_z]
    def getUnitfullTaus(self):
        #print 'Units on ts are ' + self.H0_units
        return [tauH0 / self.H0 for tauH0 in self.tauH0_of_z]
    
    def getUnitfulldLs(self):
        #print 'Units on dL are ' + self.dL_units
        return [dL * self.dL_scaling for dL in self.dL_of_z]

    def getMus(self):
        #print 'dL_units are in ' + self.dL_units + ".  They should be mega_parsec."
        return 25.0 + 5.0 * np.log10(self.getUnitfulldLs()[1:]) #ignore first dL, as it is 0 (today) 
                                    
    #You should give the X function as a function of t, tau, and z (in that order).
    # If you only want it to be a function of only one or two of them, just make
    # the output independent of the variables that are not of interest.  
    def __init__(self, XOfFunction = lambda ts, taus, zs: 1.0 if type(ts) in [float, int] else -1.0 + np.zeros(np.shape(ts)), 
                 dL0 = 0.0, t0 = 0.0, tau0 = 0.0, 
                 canon_mus = None, a0 = None, OmM0 = None, OmL0 = None, OmR0 = None, Om0 = None, H0 = None, initial_zs = None, H0_units = 'mega_years', dL_units = 'mega_parsec'):
        astro_archive = AstronomicalParameterArchive() 
        cosmo_archive = CosmologicalParameterArchive ()
        if OmM0 is None:
            OmM0 = cosmo_archive.getOmegaM()[0]
        self.OmM0 = OmM0
        if Om0 is None:
            Om0 = cosmo_archive.getOmega0()[0]
        self.Om0 = Om0
        if OmR0 is None:
            OmR0 = cosmo_archive.getOmegaR()[0]
        self.OmR0 = OmR0
        #self.OmR0 = 0.0 
        if a0 is None:
            a0 = 1.0
        self.a0 = a0
        self.t0 = t0
        self.tau0 = tau0
        self.dL0 = dL0
        if H0 is None:
            if '_' in H0_units:
                H0_prefix, H0_base_units = H0_units.split('_')
            else:
                H0_base_units = H0_units
                H0_prefix = 'unity'
            H0 = cosmo_archive.getH0(units = H0_base_units)[0]
            H0 = H0 * 1.0 * astro_archive.getPrefix(H0_prefix)

            
        self.H0 = H0
        self.H0_units = H0_units
        if self.H0_units.lower() in [ 'mega_years', 'megayears','myears','myrs','myr','my','megayrs','megayr','megay']:
            self.H0_km_s_Mpc = (H0 * astro_archive.getSecondToYear() / astro_archive.getPrefix('mega')
                                * astro_archive.getPrefix('mega') / astro_archive.getPrefix('kilo') * 1.0 * astro_archive.getParsecToM())
        #print 'self.H0_km_s_Mpc = ' + str(self.H0_km_s_Mpc) 

        #Default units are chosen such that c / H0 is given in Mpc
        self.dL_scaling = cosmo_archive.getc() / self.H0_km_s_Mpc
        self.dL_units = dL_units 
        #print 'self.dL_scaling = ' + str(self.dL_scaling) 

        self.OmL0 = self.Om0 - self.OmR0 - self.OmM0
        #print '[self.H0, self.OmM0, self.OmR0, self.OmL0] = ' + str([self.H0, self.OmM0, self.OmR0, self.OmL0]) 
        self.X0 = XOfFunction(self.t0, self.tau0, 0.0) 
        
        self.XOfFunction = XOfFunction

        self.states = []
        self.derivs = []
            

        if not(initial_zs) is None:
            self.runCalculationWithGivenzs(initial_zs)
        else:
            self.states = []
            self.a_of_t = []
            self.OmL_of_t = []
            self.H_of_t = []
            self.X_of_t = [] 
            self.mu_of_t = []
            
 

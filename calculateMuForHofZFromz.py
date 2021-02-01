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

#This one is pretty easy, since once we have H, we are golden
def computeMusForHofZFromz(scaling_of_H_of_z, zs_to_comp, H0 = None):
    cosmo_archive = CosmologicalParameterArchive() 
    OmM0 = cosmo_archive.getOmegaM()[0]
    Om0 = cosmo_archive.getOmega0()[0]
    OmR0 = cosmo_archive.getOmegaR()[0]
    OmL0 = Om0 - OmR0 - OmM0
    
    H_of_z = lambda z: scaling_of_H_of_z (z) * np.sqrt(OmM0 * (1.0 + z) ** 3.0 + OmL0 + OmR0 * (1.0 + z) ** 4.0)
    cosmo_archive = CosmologicalParameterArchive() 
    astro_archive = AstronomicalParameterArchive()
    H0_km_s_Mpc = cosmo_archive.getH0()[0]
    print 'H0_km_s_Mpc = ' + str(H0_km_s_Mpc) 
    if H0 is None:
        H0_prefix, H0_base_units = 'mega_years'.split('_')
        H0 = cosmo_archive.getH0(units = H0_base_units)[0]
        H0 = H0 * 1.0 * astro_archive.getPrefix(H0_prefix)
    print 'H0 = ' + str(H0) 
    taus = [quad(lambda z_prime: 1.0 / H_of_z(z_prime), 0, z)[0] for z in zs_to_comp]
    dLs = [(1.0 + z) * quad(lambda z_prime: 1.0 / H_of_z(z_prime), 0, z)[0] for z in zs_to_comp]
    dL_scaling = cosmo_archive.getc() / H0_km_s_Mpc
    unitfullTaus = np.array(taus) * 1.0 / H0
    unitfulldLs = np.array(dLs) * dL_scaling 
    mus = 25.0 + 5.0 * np.log10(unitfulldLs)

    return [zs_to_comp, unitfullTaus, unitfulldLs, mus]

class ResidualMuCalculatorForArbitraryWofT: 

    def BoatFishSystem(state, t):
        fish, boat = state
        d_fish = fish * (2 - boat - fish)
        d_boat = -boat * (1 - 1.5 * fish)
        return [d_fish, d_boat]

    def getBoatFishState(init_state, t):
        return odeint(BoatFishSystem, init_state, t)

    def getH(self, z):
        #print 'a = ' + str(a)
        #print 'OmL = ' + str(OmL)
        return self.HOfFunction(z)

    def tOmegaLumDistSystem(self, state, z):
        t, tau, dL = state
        self.states = self.states + [state]
        #t, OmL = state
        #print '[z, t, OmL] = ' + str([z,t,OmL])
        H = self.getH(z)
        #print 'H = ' + str(H)
        #print 'wLambda = ' + str(wLambda)
        d_t = -1.0 * (-1.0 / (1.0 + z) * 1.0 / H)
        d_tau = 1.0 / H 
        d_OmL = 3.0 * OmL * (1 + wLambda) * 1.0 / (1.0 + z)
        d_dL = (dL / (1.0 + z)) + ((1.0 + z) / H)
        d_phi = d_t * theta
        d_theta = d_t * (- 3.0 * H * theta - self.getDVOfPhi(phi))
        derivs = [d_t, d_tau, d_OmL, d_dL, d_phi, d_theta]
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
        self.states = self.gettOmegaLumDistSystem([self.t0, self.tau0, self.OmL0, self.dL0, self.phi0, self.theta0], self.zs)
        #self.states = self.gettOmegaLumDistSystem([self.t0, self.OmL0], self.zs)
        self.tH0_of_z = [state[0] for state in self.states]
        self.tauH0_of_z = [state[1] for state in self.states]
        self.OmL_of_z = [state[2] for state in self.states]
        self.dL_of_z = [state[3] for state in self.states]
        self.phi_of_z = [state[4] for state in self.states]
        self.theta_of_z = [state[5] for state in self.states]
        self.H_of_z = self.getH(np.array(self.zs), np.array(self.OmL_of_z))
        self.a_of_z = 1.0 / (1.0 + np.array(self.zs))
        self.ws_of_z = [self.getwLambda(self.tH0_of_z[i], self.tauH0_of_z[i], self.zs[i], self.phi_of_z[i], self.theta_of_z[i]) for i in range(len(self.zs))]

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
                                    
    #You should give the w function as a function of t, tau, and z (in that order).
    # If you only want it to be a function of only one or two of them, just make
    # the output independent of the variables that are not of interest.  
    def __init__(self, HOfFunction = lambda ts, taus, zs: 70.0 if type(ts) in [float, int] else 70.0 + np.zeros(np.shape(ts)), 
                 dL0 = 0.0, t0 = 0.0, tau0 = 0.0, phi0 = 1.0, theta0 = 0.0,
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
        self.OmR0 = 0.0 
        if a0 is None:
            a0 = 1.0
        self.a0 = a0
        self.t0 = t0
        self.tau0 = tau0
        self.dL0 = dL0
        self.theta0 = theta0
        self.phi0 = phi0 
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

        #self.wOfFunction = wOfFunction
        #self.VOfPhiFunction = VOfPhiFunction 
        #self.dVOfPhiFunction = dVOfPhiFunction
        #self.wFromPhi = wFromPhi
        self.HOfFunction = HOfFunction

        self.states = []
        self.derivs = []
            

        if not(initial_zs) is None:
            self.runCalculationWithGivenzs(initial_zs)
        else:
            self.states = []
            self.a_of_t = []
            self.OmL_of_t = []
            self.H_of_t = []
            self.mu_of_t = []
            self.phi_of_t = [] 
            self.theta_of_t = []
            
 

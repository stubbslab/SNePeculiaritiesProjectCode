#Solves the set of coupled ODEs acquired when w_lambda is allowed to be an arbitrary function of time.
# The coupled ODEs are of a(t) and rho_Lambda(t)
# H(t) can then be calculated numerically from this result. 


import numpy as np
import math
from scipy.integrate import odeint
from scipy.optimize import fsolve
from scipy.optimize import minimize_scalar
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

def functToMinimizePrintLoop(funct, val):
    funct_val = funct(val)
    print 'For val ' + str(val) + ' funct(val) = ' + str(funct_val)
    return funct_val
    
def getOscillatingGCosmologyWithProperDeltaGOverG(A, nu, alpha, initial_a_s, C_power_bounds, 
                                          t0 = None, tau0 = None, 
                                          canon_mus = None, OmM0 = None, OmL0 = None, OmR0 = None, Om0 = None, H0 = None,
                                          H0_units = 'mega_years', dL_units = 'mega_parsec'):
    astro_archive = AstronomicalParameterArchive() 
    cosmo_archive = CosmologicalParameterArchive ()
    if OmM0 is None:
        OmM0 = cosmo_archive.getOmegaM()[0]
    if Om0 is None:
        Om0 = cosmo_archive.getOmega0()[0]
    if OmR0 is None:
        OmR0 = cosmo_archive.getOmegaR()[0]

    target_OmL0 = Om0 - OmR0 - OmM0
    #print ('[OmM0, Om0, OmR0, target_OmL0] = ' + str([OmM0, Om0, OmR0, target_OmL0]))

    funct_to_minimize = lambda C_power: abs(ResidualMuCalculatorForMonodromicDEType2(A, nu, alpha, initial_a_s, C = 10.0 ** C_power, t0 = t0, tau0 = tau0, canon_mus = canon_mus, OmM0 = OmM0, OmL0 = OmL0, OmR0 = OmR0, Om0 = Om0, H0 = H0, H0_units = H0_units, dL_units = dL_units).rhoDE_of_a[-1] - target_OmL0)

    best_fit_C_power_res = minimize_scalar(lambda C_power: functToMinimizePrintLoop(funct_to_minimize, C_power), bounds = C_power_bounds, method='bounded')
    best_fit_C_power = best_fit_C_power_res['x']
    print 'best_fit_C_power = ' + str(best_fit_C_power)

    return ResidualMuCalculatorForMonodromicDEType2(A, nu, alpha, initial_a_s, C = 10.0 ** best_fit_C_power, t0 = t0, tau0 = tau0, canon_mus = canon_mus, OmM0 = OmM0, OmL0 = OmL0, OmR0 = OmR0, Om0 = Om0, H0 = H0, H0_units = H0_units, dL_units = dL_units)


#Differential equations in z, calculated backwards in universe time 
class ResidualMuCalculatorForOscillatingG: 

    def BoatFishSystem(state, t):
        fish, boat = state
        d_fish = fish * (2 - boat - fish)
        d_boat = -boat * (1 - 1.5 * fish)
        return [d_fish, d_boat]

    def getBoatFishState(init_state, t):
        return odeint(BoatFishSystem, init_state, t)

    def getH(self, a, theta, phi):
        #print 'a = ' + str(a)
        #print 'OmL = ' + str(OmL)
        #OmL = self.getOmL(theta, phi)
        kappa = self.getKappa(phi) 
        #print 'type(X) = ' + str(type(X)) 
        #print 'X = ' + str(X)
        #print 'type(X) in [list]) = ' + str(type(X) in [list]) 
        if (type(kappa) in [list]): kappa = np.array(kappa) 
        #print 'type(z) = ' + str(type(z)) 
        #print 'z = ' + str(z)
        x = theta / phi
        cosmo_omega = a ** (-3.0) * self.OmM0 + a ** (-4.0) * self.OmR0 + self.OmL0 
        #print '[self.OmM0, self.OmR0, self.OmL0] = ' + str([self.OmM0, self.OmR0, self.OmL0])
        H = - x / 2.0 + np.sqrt((x / 2.0) ** 2.0 + 4.0 * (1.0 / (kappa * 6.0) * (x) ** 2.0 + cosmo_omega / phi ) )
        return H 
    
    def getKappa(self, phi):
        #print 'tau = ' + str(tau)
        #print 'type(tau) = ' + str(type(tau))
        #print "type(tau) in ['list',np.ndarray] = " + str(type(tau) in ['list',np.ndarray])
        if not(type(phi) in [list, np.ndarray]): phi= float(phi)
        kappa =  self.kappaOfFunction(phi)
        #print 'kappa = ' + str(kappa)
        if not (abs(kappa) > 0.0):
            new_kappa = self.prev_kappas[1] + self.prev_kappas[1] - self.prev_kappas[0]
            print 'bad kappa = ' + str(kappa) + ' ; self.prev_kappas = ' + str(self.prev_kappas) + ' ; new kappa = ' + str(new_kappa) 
            kappa = new_kappa 
        self.prev_kappas[0] = self.prev_kappas[1]
        self.prev_kappas[1] = kappa
        return kappa 

    def getDKappaDPhi(self, phi): 
        if not(type(phi) in [list, np.ndarray]): phi= float(phi)
        return self.dKappaDPhiOfFunction(phi)

    def getG(self, phi):
        kappa = self.getKappa(phi)
        if not(type(phi) in [list, np.ndarray]): phi= float(phi)
        G = 1.0 / phi * (4.0 * kappa + 2.0) / (3.0 * kappa + 2.0) if abs(phi) > 0.0 else 0.0 
        return G
    
    def getDG(self, phi, kappa, d_kappa_d_phi, d_phi):
        d_G = -d_phi / (phi ** 2.0) * (4.0 * kappa + 2.0) / (3.0 * kappa + 2.0) + 1.0 / phi * d_kappa_d_phi * d_phi * (4.0 * (3.0 * kappa + 2.0) - (4.0 * kappa + 2.0) * 3.0) / ((3.0 * kappa + 2.0) ** 2.0)
        return d_G
    
    def getDThetaDT(self, z, H, kappa, d_kappa_d_phi, theta):
        #if H == 1.0: print 'At H == 1.0, [z, H, omega, phi] =' + str([z, H, omega, phi])
        #if H == 1.0: print 'At H == 1.0, [z, H, omega, phi] = ' + str([z, H, omega, phi])
        #print '((3.0 * phi * H / omega) ** 2.0 + 6.0 / omega * (H ** 2.0 * phi ** 2.0 - (self.OmM0 * (1.0 + z) ** 3.0 + self.OmR0 * (1.0 + z) ** 4.0 + self.OmL0) * phi )) = ' + str(((3.0 * phi * H / omega) ** 2.0 + 6.0 / omega * (H ** 2.0 * phi ** 2.0 - (self.OmM0 * (1.0 + z) ** 3.0 + self.OmR0 * (1.0 + z) ** 4.0 + self.OmL0) * phi )))

        if kappa > 0.0: 
            d_theta_d_t = (3.0 * kappa) / (2.0 + 3.0 * kappa) * (self.OmM0 * (1.0 + z) ** 3.0 + 4.0 * self.OmL0) - theta * (3.0 * H - d_kappa_d_phi / kappa * theta / (2.0 + 3.0 * kappa))
        else:
            d_theta_d_t = 0.0 
        #if H == 1.0: print 'At H == 1.0, [z, H, omega, phi, d_phi_d_t] =' + str([z, H, omega, phi, d_phi_d_t])
        return d_theta_d_t

    def getDHDT(self, H, phi, z, kappa, d_kappa_d_phi, theta): 
        #d_t_to_d_z_scaling = - (1.0 + z) * H
        if kappa > 0.0:
            d_H_term1 = - (2.0 * self.OmR0 * (1.0 + z) ** 4.0 + self.OmM0 * (1.0 + z) ** 3.0 - 2.0 * self.OmL0) / (phi * (2.0 + 3.0 * kappa)) 
            d_H_term2 = - 3.0 * kappa * (self.OmR0 * (1.0 + z) ** 4.0 + self.OmM0 * (1.0 + z) ** 3.0 + self.OmL0) / (phi * (2.0 + 3.0 * kappa)) 
            d_H_term3 = -1.0 / 2.0 * d_kappa_d_phi / kappa * theta / (2.0 + 3.0 * kappa) * theta / phi
            d_H_term4 = - H ** 2.0
            d_H_term5 = - 1.0 / (3.0 * kappa) * theta ** 2.0 / phi ** 2.0
            d_H_term6 = H * theta / phi
            d_H_terms = [d_H_term1, d_H_term2, d_H_term3, d_H_term4, d_H_term5, d_H_term5, d_H_term6, self.d_t_to_d_z_scaling(z, H)] 
            #self.d_H_terms = self.d_H_terms + [[d_H_term1, d_H_term2, d_H_term3, d_H_term4, d_H_term5, d_H_term5, d_H_term6, self.d_t_to_d_z_scaling(z, H)]]
            d_H_d_t = (d_H_term1 + d_H_term2 + d_H_term3 + d_H_term4 + d_H_term5 + d_H_term5 + d_H_term6)
            #if H == 1.0: print 'At H == 1.0, [d_H_d_t, d_H, self.d_t_to_d_z_scaling(z, H), phi, d_phi_d_t]  = ' + str([d_H_d_t, d_H, d_H_terms, phi, d_phi_d_t])
        else:
            d_H_term1 = -1.0 / 2.0 * (4.0 * self.OmR0 * (1.0 + z) ** 4.0 + 3.0 * self.OmM0 * (1.0 + z) ** 3.0) 
            d_H_terms = [d_H_term1] + [0.0 for i in range(6-1)] + [self.d_t_to_d_z_scaling(z, H)]
            d_H_d_t = d_H_term1
        return d_H_d_t, d_H_terms 
        
    def tPhiToGSystem(self, state, z):
        #t, tau, dL, phi, H = state
        t, tau, phi, theta, H, dL = state
        self.states = self.states + [state]
        #t, OmL = state
        #print '[z, t, OmL] = ' + str([z,t,OmL])
        kappa = self.getKappa(phi) 
        #H = self.getH(a, theta, phi)
        #derivative in z 
        d_t = -1.0 * (-1.0 / (1.0 + z) * 1.0 / H)
        d_tau = 1.0 / H
        d_phi_d_t = theta
        d_phi = d_phi_d_t * self.d_t_to_d_z_scaling(z, H)
        d_kappa_d_phi = self.getDKappaDPhi(phi) 
        d_theta_d_t = self.getDThetaDT(z, H, kappa, d_kappa_d_phi, theta)
        d_theta = d_theta_d_t * self.d_t_to_d_z_scaling(z, H)
        G = self.getG(phi)
        d_G = self.getDG(phi, kappa, d_kappa_d_phi, d_phi) 
        if self.change_SN_phys:
            #d_dL = - (self.alpha / 2.0) * self.getDLRatio(z, G, dG) / self.getLRatio(z, G) * dL + (dL / (1.0 + z)) + (1.0 + z) / (H * G ** (self.alpha / 2.0))
            d_dL = - (self.alpha_SN / 2.0) * d_G / G * dL + (dL / (1.0 + z)) + (1.0 + z) / (H * G ** (self.alpha_SN / 2.0))
        else:
            d_dL = (dL / (1.0 + z)) + ((1.0 + z) / H)
        #if H == 1.0:
        #    print '[phi, d_phi_d_t, H, z, self.getDOmegaDT(phi, d_phi_d_t)] = ' + str([phi, d_phi_d_t, H, z, self.getDOmegaDT(phi, d_phi_d_t)])
        #    print '[H, phi, z, omega, d_omega_d_t, d_phi_d_t] = ' + str([H, phi, z, omega, d_omega_d_t, d_phi_d_t])
        d_H_d_t, d_H_terms = self.getDHDT(H, phi, z, kappa, d_kappa_d_phi, theta)
        d_H = d_H_d_t * self.d_t_to_d_z_scaling(z, H)  
        self.d_H_terms = self.d_H_terms + [d_H_terms]
        
        #d_dL = (dL / (1.0 + z)) + ((1.0 + z) / H)
        derivs = [d_t, d_tau, d_phi, d_theta, d_H, d_dL]
        self.used_zs = self.used_zs + [z] 
        self.derivs = self.derivs + [derivs]
        return derivs 
        #print '[d_t, d_OmL] = ' + str([d_t, d_OmL]) 
        #return [d_t, d_OmL]

    def gettPhiToGSystem(self, init_state, z):
        return odeint(self.tPhiToGSystem, init_state, z)

    def runCalculationWithGivenzs(self, zs):
        #print 'zs = ' + str(zs)
        #if zs[0] > 0.0: zs = [0.0] + zs #Initial values are known at z = 0.0
        #self.zs = zs
        self.zs = zs 
        self.states = self.gettPhiToGSystem([self.t0, self.tau0, self.phi0, self.theta0, self.H_init, self.dL_init], self.zs)
        #self.states = self.gettOmegaLumDistSystem([self.t0, self.OmL0], self.zs)
        self.tH0_of_z = [state[0] for state in self.states]
        self.tauH0_of_z = [state[1] for state in self.states]
        self.phi_of_z = [state[2] for state in self.states]
        self.d_phi_of_z = [state[3] for state in self.states]
        self.H_of_z = [state[4] for state in self.states]
        self.dL_of_z = [state[5] for state in self.states]
        self.G_of_z = [self.getG(phi) for phi in self.phi_of_z]
        self.kappa_of_z = [self.kappaOfFunction(phi) for phi in self.phi_of_z]
        self.d_kappa_d_phi_of_z = [self.getDKappaDPhi(phi) for phi in self.phi_of_z ]
        self.d_G_of_z = [self.getDG(self.phi_of_z[i], self.kappa_of_z[i], self.d_kappa_d_phi_of_z[i], self.d_phi_of_z[i]) for i in range(len(self.zs))]
        self.d_G_d_t_of_z = [self.d_G_of_z[i] * 1.0 / self.d_t_to_d_z_scaling(self.zs[i], self.H_of_z[i]) for i in range(len(self.zs))]
        #self.dV_of_a = [self.dVOfFunction(phi) for phi in self.phi_of_a]
        #self.H_of_z = self.getH(np.array(self.zs), np.array(self.theta_of_z), np.array(self.phi_of_z))
        #self.H_of_a = self.getH(np.array(self.tH0_of_a), np.array(self.tauH0_of_a), np.array(self.zs))
        self.a_of_z = 1.0 / ((np.array(self.zs)) + 1.0) 
        #self.X_of_a = [self.getX(theta) for theta in self.theta_of_a] 
        #self.w_of_a = [(1.0 - X) / (1.0 - 3.0 * X) for X in self.X_of_a]
        #self.rhoDE_of_a = [self.getOmL(self.theta_of_a[i], self.phi_of_a[i]) for i in range(len(self.zs))]
        #self.pDE_of_a = [self.VOfFunction(self.phi_of_a[i]) * (-self.X_of_a[i] + self.X_of_a[i] ** 2.0) for i in range(len(self.zs))]
        #self.w_of_a = [self.pDE_of_a[i] / self.rhoDE_of_a[i] for i in range(len(self.zs))]
        self.H_of_z_interp = interp1d (self.zs, self.H_of_z)
        #self.dL_of_z = [(1.0 + z) * quad(lambda z_int_var: 1.0 / (self.H_of_z_interp(z_int_var)), 0, z)[0] for z in self.z_of_a] 
        #self.a_of_z = 1.0 / (1.0 + np.array(self.zs))

    def getUnitfullTs(self):
        #print 'Units on ts are ' + self.H0_units
        return [tH0 / self.H0 for tH0 in self.tH0_of_a]
    
    def getUnitfullTaus(self):
        #print 'Units on ts are ' + self.H0_units
        return [tauH0 / self.H0 for tauH0 in self.tauH0_of_a]

    def getUnitlessdLs(self, target_zs = None):
        if target_zs is None:
            target_zs = self.z_of_a
        self.dL_of_z = [(1.0 + z) * quad(lambda z_int_var: 1.0 / (self.H_of_z_interp(z_int_var)), 0, z)[0] for z in target_zs]
        return self.dL_of_z 
    
    def getUnitfulldLs(self):
        #print 'Units on dL are ' + self.dL_units
        #if not(target_zs is None) or not(self.dL_of_z is None):
        #    unitless_dLs = self.getUnitlessdLs(target_zs = target_zs)
        return [dL * self.dL_scaling for dL in self.dL_of_z]

    def getMus(self, target_zs = None):
        #print 'dL_units are in ' + self.dL_units + ".  They should be mega_parsec."
        return np.array([0.0] + (25.0 + 5.0 * np.log10(self.getUnitfulldLs()[1:])).tolist()) #ignore first dL, as it is 0 (today) 
                                    
    #You should give the X function as a function of t, tau, and z (in that order).
    # If you only want it to be a function of only one or two of them, just make
    # the output independent of the variables that are not of interest.  
    def __init__(self, A, nu, alpha, psi, initial_zs, 
                       t0 = None, tau0 = None, C = None, phi0 = None, eps = None, 
                       canon_mus = None, OmM0 = None, OmL0 = None, OmR0 = None, Om0 = None, H0 = None, 
                       H0_units = 'mega_years', dL_units = 'mega_parsec', change_SN_phys = 1, alpha_SN = -3.0 / 2.0 ):
        #Monodromic DE has a specific potential.  So that is hard coded here; not something that a use gives
        phi_var_term = lambda phi, C, eps, phi0, alpha, nu, A, psi: C * abs(1.0 + eps - phi /phi0) ** (alpha) * np.sin(nu * phi / phi0 + psi)
        kappaOfSinglePhiFunction = lambda phi, C, eps, phi0, alpha, nu, A, psi: (2.0 * phi_var_term (phi, C, eps, phi0, alpha, nu, A, psi)) / (1.0 - 3.0 * phi_var_term (phi, C, eps, phi0, alpha, nu, A, psi)) 
        #omegaOfSinglePhiFunction = lambda phi, B, phi0, alpha, nu, A, psi: (2.0 * B * (1 + eps - phi /phi0) ** (-alpha) * (1 - A * np.sin(nu * phi / phi0 + psi)) - 2.0 ) / 3.0
        kappaOfFunction = lambda phis, C, eps, phi0, alpha, nu, A, psi: [kappaOfSinglePhiFunction(phi, C, eps, phi0, alpha, nu, A, psi) for phi in phis] if type(phis) in [list, np.ndarray] else kappaOfSinglePhiFunction (phis, C, eps, phi0, alpha, nu, A, psi)

        self.d_t_to_d_z_scaling = lambda z, H: 1.0 / (- (1.0 + z) * H ) # d_d_z = self.d_t_to_d_z_scaling * d_d_t
        d_phi_var_term = lambda phi, C, eps, phi0, alpha, nu, A, psi: -C * alpha * abs(1.0 + eps - phi /phi0) ** (alpha - 1.0) * np.sin(nu * phi / phi0 + psi) + C * abs(1.0 + eps - phi /phi0) ** (alpha) * nu / phi0 * np.cos(nu * phi / phi0 + psi)
        dKappaDPhiOfSinglePhiFunction = lambda phi, C, eps, phi0, alpha, nu, A, psi: ((2.0 * d_phi_var_term(phi, C, eps, phi0, alpha, nu, A, psi) ) * (1.0 - 3.0 * phi_var_term (phi, C, eps, phi0, alpha, nu, A, psi)) - (2.0 * phi_var_term (phi, C, eps, phi0, alpha, nu, A, psi)) * (-3.0 * d_phi_var_term(phi, C, eps, phi0, alpha, nu, A, psi)) ) / ((1.0 - 3.0 * phi_var_term (phi, C, eps, phi0, alpha, nu, A, psi)) ** 2.0)
        dKappaDPhiOfFunction = lambda phis, C, eps, phi0, alpha, nu, A, psi: ([dKappaDPhiOfSinglePhiFunction(phis[i], C, eps, phi0, alpha, nu, A, psi) for i in range(len(phis))] if type(phis) in [list,np.ndarray]
                                                                                      else dKappaDPhiOfSinglePhiFunction (phis, C, eps, phi0, alpha, nu, A, psi) ) 
    
        #VOfFunction = lambda phis, C, phi0, alpha, nu, A: [C * (phi / phi0) ** (-alpha) * (1.0 - A * np.sin(nu * phi)) if phi > 0.0 else 0.0 for phi in phis] if type(phis) in ['list',np.ndarray] else C * (phis / phi0) ** (-alpha) * (1.0 - A * np.sin(nu * phis)) if phis > 0.0 else 0.0
        #dVOfFunction = lambda phis, C, phi0, alpha, nu, A: [C * (-alpha) / phi0 * (phi / phi0) ** (-alpha - 1.0) * (1.0 - A * np.sin(nu * phi)) + C * (phi / phi0) ** (-alpha) * (- A * nu * np.cos(nu * phi)) if phi > 0.0 else 0.0 for phi in phis] if type(phis) in ['list',np.ndarray] else C * (-alpha) / phi0 * (phis / phi0) ** (-alpha - 1.0) * (1.0 - A * np.sin(nu * phis)) + C * (phis / phi0) ** (-alpha) * (- A * nu * np.cos(nu * phis)) if phis > 0.0 else 0.0
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

        self.change_SN_phys = change_SN_phys
        self.alpha_SN = alpha_SN 

        self.zs = initial_zs
        self.a0 = self.zs[0] 
        if t0 is None:
            #print self.OmR0 
            #print '[a0, t0, self.OmM0, self.OmR0] = ' + str([self.a0, t0, self.OmM0, self.OmR0])
            t0 = quad(lambda a: 1.0 / (a * (self.OmM0 * a ** -3.0 + self.OmR0 * a ** -4.0) ** 0.5), 0.0, self.a0)[0] - 0.0
            #print '[a0, t0] = ' + str([self.a0, t0])
            t0 = 0.0 
        self.t0 = t0
        if tau0 is None:
            tau0 = quad(lambda a: 1.0 / (a ** 2.0 * (self.OmM0 * a ** -3.0 + self.OmR0 * a ** -4.0) ** 0.5), 0.0, self.a0)[0] - 0.0
            tau0 = 0.0
        self.tau0 = tau0
        #self.dL0 = dL0
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
        self.alpha = alpha
        self.A = A
        self.nu = nu
        self.psi = psi
        
        if phi0 is None:
            phi0 = 1.0
        self.phi0 = phi0
        if C is None:
            C = 1.0
        self.C = C
        if eps is None:
            eps = 10.0 ** -4.0
        self.eps = eps 
        H_init = 1.0
        self.H_init = H_init
        dL_init = 0.0
        self.dL_init = dL_init 

        print '[self.C, self.phi0, self.alpha, self.nu, self.A, self.psi] = ' + str([self.C, self.phi0, self.alpha, self.nu, self.A, self.psi])
        #print 'omegaOfFunction([0.9, 1.0], self.B, self.phi0, self.alpha, self.nu, self.A, self.psi) = ' + str(omegaOfFunction([0.9, 1.0], self.B, self.phi0, self.alpha, self.nu, self.A, self.psi))
        self.kappaOfFunction = lambda phis: kappaOfFunction(phis, self.C, self.eps, self.phi0, self.alpha, self.nu, self.A, self.psi)
        self.dKappaDPhiOfFunction = lambda phis: dKappaDPhiOfFunction(phis, self.C, self.eps, self.phi0, self.alpha, self.nu, self.A, self.psi) 
        #self.dVOfFunction = lambda phis: dVOfFunction(phis, self.B, self.phi0, self.alpha, self.nu, self.A)
        z_init = initial_zs[0] 
        
        #V_init = self.VOfFunction(phi0)
        kappa0 = self.kappaOfFunction(phi0)
        #theta0 = self.getDPhiDT(initial_zs[0], H_init, omega0, phi0)
        theta0 = 6.0 * kappa0
        self.theta0 = theta0 
        d_kappa_d_phi0 = self.dKappaDPhiOfFunction(phi0)
        d_theta_d_t0 = self.getDThetaDT(initial_zs[0], H_init, kappa0, d_kappa_d_phi0, theta0) 
        d_phi0 = theta0 * self.d_t_to_d_z_scaling(z_init, H_init)
        d_theta0 = d_theta_d_t0 * self.d_t_to_d_z_scaling(z_init, H_init) 
        
        #print '[phi0, d_phi_d_t0, H_init, z_init, self.getDOmegaDT(phi0, d_phi_d_t0, H_init, z_init), self.dOmegaDTOfFunction(phi0, d_phi_d_t0, H_init, z_init)] = ' + str([phi0, d_phi_d_t0, H_init, z_init, self.getDOmegaDT(phi0, d_phi_d_t0), self.dOmegaDTOfFunction(phi0, d_phi_d_t0)]) 
        #print '[H_init, phi0, initial_zs[0], omega0, d_omega0, d_phi_d_t0] = ' + str([H_init, phi0, initial_zs[0], omega0, d_omega_d_t0, d_phi_d_t0])
        d_H_d_t_init, d_H_init_terms = self.getDHDT(H_init, phi0, initial_zs[0], kappa0, d_kappa_d_phi0, theta0)
        d_H_init = d_H_d_t_init * self.d_t_to_d_z_scaling(z_init, H_init) 
        
        #d_theta0 = -3.0 * theta0 / self.a0 - 1.0 / (self.a0 * H_init) * self.dVOfFunction(phi0)
        #d_theta0 = - (( theta0 * self.dVOfFunction(phi0) * (-X0 + 3.0 * X0 ** 2.0) + 6.0 * H_init * self.VOfFunction(phi0) * (-X0 + 2.0 * X0 ** 2.0) )
        #             / ( self.VOfFunction(phi0) * theta0 * (-1.0 + 6.0 * X0) * self.a0 * H_init) )
        print '[A, nu, alpha, psi, eps, kappa0, phi0, theta0, d_phi0, d_theta_d_t0, d_theta0, H_init, d_H_d_t_init, d_H_init, z_init, self.kappaOfFunction(phi0),  self.dKappaDPhiOfFunction(phi0)] = ' + str([A, nu, alpha, psi, eps, kappa0, phi0, theta0, d_phi0, d_theta_d_t0, d_theta0, H_init, d_H_d_t_init, d_H_init, z_init, self.kappaOfFunction(phi0),  self.dKappaDPhiOfFunction(phi0)])

        self.states = []
        self.derivs = []
        self.used_zs = []
        self.d_H_terms = [] 

        self.prev_kappas = [kappa0 - d_kappa_d_phi0 * 0.00001, kappa0] 
        self.runCalculationWithGivenzs(self.zs)
            
 

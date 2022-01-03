#Solves the set of coupled ODEs acquired when G derives from an evolving scalar field. 
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
import matplotlib.pyplot as plt
from matplotlib import ticker, cm, colors

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
        #OmL = self.getOmL(theta, phi)
        kappa = self.getKappa(phi) 
        if (type(kappa) in [list]): kappa = np.array(kappa) 
        x = theta / phi
        cosmo_omega = a ** (-3.0) * self.OmM0 + a ** (-4.0) * self.OmR0 + self.OmL0 
        H = - x / 2.0 + np.sqrt((x / 2.0) ** 2.0 + 4.0 * (1.0 / (kappa * 6.0) * (x) ** 2.0 + cosmo_omega / phi ) )
        return H 
    
    def getKappa(self, phi):
        if not(type(phi) in [list, np.ndarray]): phi= float(phi)
        return self.kappaOfFunction(phi)

    def getDKappaDPhi(self, phi): 
        if not(type(phi) in [list, np.ndarray]): phi= float(phi)
        return self.dKappaDPhiOfFunction(phi)

    def getG(self, phi):
        kappa = self.getKappa(phi)
        if not(type(phi) in [list, np.ndarray]): phi= float(phi)
        inv_phi = 1.0 / phi
        kappa_funct = (4.0 * kappa + 2.0) / (3.0 * kappa + 2.0)
        G = inv_phi * kappa_funct 
        return G
    
    def getDG(self, phi, kappa, d_kappa_d_phi, d_phi):
        d_G = -d_phi / (phi ** 2.0) * (4.0 * kappa + 2.0) / (3.0 * kappa + 2.0) + 1.0 / phi * d_kappa_d_phi * d_phi * (2.0) / ((3.0 * kappa + 2.0) ** 2.0)
        return d_G

    def getDDG_DT2(self, phi, theta, kappa, G, dG, z, H, dKappaDPhi):
        term1 = (theta / phi - theta ** 2.0 / phi ** 2.0) * (-G + 2.0 / ((3.0 * kappa + 2.0) ** 2.0) )
        term2 = theta / phi * (dG * (1.0 + z) * H - 8.0 / ((3.0 * kappa + 2.0) ** 2.0) * dKappaDPhi * theta)

        ddG_dt2 = term1 + term2
        return ddG_dt2 
    
    def getDThetaDT(self, z, H, kappa, d_kappa_d_phi, theta):
        if kappa > 0.0:
            #print '[d_kappa_d_phi, kappa] = ' + str([d_kappa_d_phi, kappa])
            d_theta_d_t = (3.0 * kappa) / (2.0 + 3.0 * kappa) * (self.OmM0 * (1.0 + z) ** 3.0 + 4.0 * self.OmL0) - theta * (3.0 * H - d_kappa_d_phi / kappa * theta / (2.0 + 3.0 * kappa))
        else:
            d_theta_d_t = 0.0 
        #if H == 1.0: print 'At H == 1.0, [z, H, omega, phi, d_phi_d_t] =' + str([z, H, omega, phi, d_phi_d_t])
        return d_theta_d_t

    def getDHDT(self, H, phi, z, kappa, d_kappa_d_phi, theta): 
        #d_t_to_d_z_scaling = - (1.0 + z) * H
        if kappa > 0.0:
            d_H_term1 = - H ** 2.0
            d_H_term2 = - 1.0 / (3.0 * kappa) * theta ** 2.0 / phi ** 2.0
            d_H_term3 = H * theta / phi
            d_H_term4 = - (2.0 * self.OmR0 * (1.0 + z) ** 4.0 + self.OmM0 * (1.0 + z) ** 3.0 - 2.0 * self.OmL0) / (phi * (2.0 + 3.0 * kappa)) 
            d_H_term5 = - 3.0 * kappa * (self.OmR0 * (1.0 + z) ** 4.0 + self.OmM0 * (1.0 + z) ** 3.0 + self.OmL0) / (phi * (2.0 + 3.0 * kappa)) 
            d_H_term6 = -1.0 / 2.0 * d_kappa_d_phi / kappa * theta / (2.0 + 3.0 * kappa) * theta / phi
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
        d_H_d_t, d_H_terms = self.getDHDT(H, phi, z, kappa, d_kappa_d_phi, theta)
        d_H = d_H_d_t * self.d_t_to_d_z_scaling(z, H)  
        self.d_H_terms = self.d_H_terms + [d_H_terms]
        
        #d_dL = (dL / (1.0 + z)) + ((1.0 + z) / H)
        derivs = [d_t, d_tau, d_phi, d_theta, d_H, d_dL]
        self.used_zs = self.used_zs + [z] 
        self.derivs = self.derivs + [derivs]
        return derivs 

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
        self.G_of_z = [self.getG(phi) if abs(phi) > 0.0 else np.inf for phi in self.phi_of_z ]
        self.kappa_of_z = [self.kappaOfFunction(phi) for phi in self.phi_of_z]
        self.d_kappa_d_phi_of_z = [self.getDKappaDPhi(phi) for phi in self.phi_of_z ]
        self.d_G_of_z = [self.getDG(self.phi_of_z[i], self.kappa_of_z[i], self.d_kappa_d_phi_of_z[i], self.d_phi_of_z[i]) if abs(self.phi_of_z[i]) > 0.0 else 0.0 for i in range(len(self.zs)) ]
        self.d_G_d_t_of_z = [self.d_G_of_z[i] * 1.0 / self.d_t_to_d_z_scaling(self.zs[i], self.H_of_z[i]) for i in range(len(self.zs)) ]
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
        return [dL * self.dL_scaling for dL in self.dL_of_z]

    def getMus(self, target_zs = None):
        #print 'dL_units are in ' + self.dL_units + ".  They should be mega_parsec."
        return np.array(25.0 + 5.0 * np.log10(self.getUnitfulldLs())) #ignore first dL, as it is 0 (today) 
                                    
    #You should give the X function as a function of t, tau, and z (in that order).
    # If you only want it to be a function of only one or two of them, just make
    # the output independent of the variables that are not of interest.  
    def __init__(self, A, nu, alpha, psi, C, initial_zs, 
                       t0 = None, tau0 = None, eps = None, 
                       canon_mus = None, OmM0 = None, OmL0 = None, OmR0 = None, Om0 = None, H0 = None, d_kappa_cap = 10000.0, 
                       H0_units = 'mega_years', dL_units = 'mega_parsec', change_SN_phys = 1, alpha_SN = -3.0 / 2.0, phi_frac_for_absolute_value = 0.01 ):
        #Monodromic DE has a specific potential.  So that is hard coded here; not something that a use gives
        #def correction_
        #phi_var_term = lambda phi, C, eps, phi0, alpha, nu, A, psi: ( (C if abs(1.0 + eps - phi /phi0) <= (1.0 + eps) * (phi_frac_for_absolute_value) else 1.0)) * abs(1.0 + eps - phi /phi0) ** (alpha if abs(1.0 + eps - phi /phi0) <= (1.0 + eps) * (phi_frac_for_absolute_value) else 1.0)
        #                                                                * (1.0 + A * np.sin(nu * phi / phi0 + psi) ** 2.0) / (2.0 * math.pi / math.sqrt(1.0 + A) )
        phi_var_term = lambda phi, C, eps, phi0, alpha, nu, A, psi: ( C * abs(1.0 + eps - phi /phi0) ** alpha
                                                                        * (1.0 + A * np.sin(nu * phi / phi0 + psi) ** 2.0) / (2.0 * math.pi / math.sqrt(1.0 + A) ) )
                                                                        
        #kappaOfSinglePhiFunction = lambda phi, C, eps, phi0, alpha, nu, A, psi: (2.0 * phi_var_term (phi, C, eps, phi0, alpha, nu, A, psi)) / (2.0 - 3.0 * phi_var_term (phi, C, eps, phi0, alpha, nu, A, psi))
        kappaOfSinglePhiFunction = lambda phi, C, eps, phi0, alpha, nu, A, psi: phi_var_term (phi, C, eps, phi0, alpha, nu, A, psi) 
        #omegaOfSinglePhiFunction = lambda phi, B, phi0, alpha, nu, A, psi: (2.0 * B * (1 + eps - phi /phi0) ** (-alpha) * (1 - A * np.sin(nu * phi / phi0 + psi)) - 2.0 ) / 3.0
        kappaOfFunction = lambda phis, C, eps, phi0, alpha, nu, A, psi: ( [kappaOfSinglePhiFunction(phi, C, eps, phi0, alpha, nu, A, psi) for phi in phis] if type(phis) in [list, np.ndarray]
                                                                           else kappaOfSinglePhiFunction (phis, C, eps, phi0, alpha, nu, A, psi) )

        self.d_t_to_d_z_scaling = lambda z, H: 1.0 / (- (1.0 + z) * H ) # d_d_z = self.d_t_to_d_z_scaling * d_d_t
        d_phi_var_term = lambda phi, C, eps, phi0, alpha, nu, A, psi: ((-1 if 1.0 + eps - phi /phi0 > 0.0 else 1)
                                                                       * min(C * alpha / phi0 * abs(1.0 + eps - phi /phi0) ** (alpha - 1.0) * (1.0 + A * np.sin(nu * phi / phi0 + psi) ** 2.0) / (1.0 + A / 2.0), d_kappa_cap) + C * abs(1.0 + eps - phi /phi0) ** (alpha) * (2.0 * A * nu / phi0 * np.sin(nu * phi / phi0 + psi) * np.cos(nu * phi / phi0 + psi)) / (2.0 * math.pi / math.sqrt(1.0 + A) ) )
        #dKappaDPhiOfSinglePhiFunction = lambda phi, C, eps, phi0, alpha, nu, A, psi: ((2.0 * d_phi_var_term(phi, C, eps, phi0, alpha, nu, A, psi) ) * (2.0 - 3.0 * phi_var_term (phi, C, eps, phi0, alpha, nu, A, psi)) - (2.0 * phi_var_term (phi, C, eps, phi0, alpha, nu, A, psi)) * (-3.0 * d_phi_var_term(phi, C, eps, phi0, alpha, nu, A, psi)) ) / ((2.0 - 3.0 * phi_var_term (phi, C, eps, phi0, alpha, nu, A, psi)) ** 2.0)
        dKappaDPhiOfSinglePhiFunction = lambda phi, C, eps, phi0, alpha, nu, A, psi: d_phi_var_term(phi, C, eps, phi0, alpha, nu, A, psi)  
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
        
        if C is None:
            C = 1.0
        self.C = C
        H_init = 1.0
        self.H_init = H_init
        dL_init = 0.0
        self.dL_init = dL_init
        
        z_init = initial_zs[0] 

        kappa0_of_eps = lambda eps: kappaOfSinglePhiFunction(1.0, self.C, eps, 1.0, self.alpha, self.nu, self.A, self.psi)
        d_kappa0_times_phi0_of_eps = lambda eps: dKappaDPhiOfSinglePhiFunction(1.0, self.C, eps, 1.0, self.alpha, self.nu, self.A, self.psi)
        if eps is None:
            eps = 10.0 ** -4.0
            dG_part_to_minimize_of_eps = lambda eps: - (4.0 * kappa0_of_eps(eps) + 2.0) / (3.0 * kappa0_of_eps(eps) + 2.0) + d_kappa0_times_phi0_of_eps(eps) * 2.0 / ((3.0 * kappa0_of_eps(eps) + 2.0) ** 2.0)
            eps_log_search_lims = [-10.0, 1.0]
            test_epses = np.logspace(eps_log_search_lims[0], eps_log_search_lims[1], 1001)
            #plt.scatter(test_epses, [abs(dG_part_to_minimize_of_eps (test_eps)) for test_eps in test_epses])
            #plt.xscale('log')
            #plt.xlim(min(test_epses), max(test_epses))
            #plt.show() 
            eps_min = minimize_scalar(lambda log_eps: abs(dG_part_to_minimize_of_eps (10.0 ** log_eps)), method = 'bounded', bounds = eps_log_search_lims )
            #print 'eps_min = ' + str(eps_min)
            eps = 10.0 ** eps_min['x']

        G_of_phi0 = lambda phi0: 1.0 / phi0 * (4.0 * kappa0_of_eps(eps) + 2.0) / (3.0 * kappa0_of_eps(eps) + 2.0)
        phi0_min = minimize_scalar(lambda phi0_val: abs(G_of_phi0(phi0_val) - 1.0), method = 'bounded', bounds = (0.5, 1.5))
        phi0 = phi0_min['x']
        #kappa0_funct = lambda = kappaOfFunction(1.0, self.C, self.eps, 1.0, self.alpha, self.nu, self.A, self.psi)
        #min_phi0 = .minimize_scalar()
        #if phi0 <= 4.0 / 5.0:
        kappa0 = kappaOfFunction(1.0, self.C, eps, 1.0, self.alpha, self.nu, self.A, self.psi)
        #If phi0 is less than this value, theta0 comes out complex.  So we now change our criterion for calculating phi0 and eps
        value_in_radical = ((6.0 * phi0 * kappa0) ** 2.0 - 4.0 * (6.0 * phi0 * kappa0 - 6.0 * phi0 * kappa0 * phi0)) 
        if value_in_radical < 0.0:
            print ('!!!!! kappa0 = ' + str(kappa0) + ', phi0 = ' + str(phi0) + ', which means the value in the radical is ' + str(value_in_radical) + '   Commencing a hard minimization !!!!!') 
            #Now I know my best value is when phi0 is equal to 2.0 / (3.0 * kappa0 + 2.0)
            kappa0_of_eps = lambda eps: kappaOfSinglePhiFunction(1.0, self.C, eps, minimize_scalar(lambda new_phi0: np.abs(kappaOfSinglePhiFunction(1.0, self.C, eps, new_phi0, self.alpha, self.nu, self.A, self.psi) + (2.0 * (new_phi0 - 1.0)) / (3.0 * new_phi0)))['x'], self.alpha, self.nu, self.psi) 
            phi0 = 4.0 / 5.0 
            G_of_eps = lambda new_eps: 1.0 / phi0  * (4.0 * kappa0_of_eps(new_eps) + 2.0) / (3.0 * kappa0_of_eps(new_eps) + 2.0)
            new_eps_min = minimize_scalar(lambda log_eps: abs(G_of_eps (10.0 ** log_eps) - 1.0), method = 'bounded', bounds = eps_log_search_lims )
            new_eps = new_eps_min['x']
            new_kappa0 = kappa0_of_eps(new_eps)
            new_phi = 2.0 / (3.0 * new_kapp0 + 2.0) 
            eps = new_eps
            kappa0 = new_kappa0
            phi0 = new_phi0
        self.eps = eps
        self.phi0 = phi0

        #kappa0 = kappaOfFunction(1.0, self.C, self.eps, 1.0, self.alpha, self.nu, self.A, self.psi)
        theta0 = 6.0 * kappa0 * self.phi0 / 2.0 * (1 + np.sqrt(5.0 - 4.0 / phi0))
        theta0_old = theta0
        theta0 = 6.0 * self.phi0 * kappa0 / 2.0 + np.sqrt(( 6.0 * self.phi0 * kappa0 ) ** 2.0 - 4.0 * ( 6.0 * self.phi0 * kappa0 ) * (1.0 - self.phi0)) / 2.0
        print ('[theta0, theta0_old] = ' + str([theta0, theta0_old]))  
        self.theta0 = theta0
        d_phi0 = theta0 * self.d_t_to_d_z_scaling(z_init, H_init)


        self.kappaOfFunction = lambda phis: kappaOfFunction(phis, self.C, eps, phi0, self.alpha, self.nu, self.A, self.psi)
        self.dKappaDPhiOfFunction = lambda phis:dKappaDPhiOfFunction(phis, self.C, eps, phi0, self.alpha, self.nu, self.A, self.psi) 
        
        #theta0 = self.getDPhiDT(initial_zs[0], H_init, omega0, phi0)
        d_kappa_d_phi0 = self.dKappaDPhiOfFunction(phi0)
        d_theta_d_t0 = self.getDThetaDT(initial_zs[0], H_init, kappa0, d_kappa_d_phi0, theta0) 
        d_theta0 = d_theta_d_t0 * self.d_t_to_d_z_scaling(z_init, H_init)
        G0 = self.getG(phi0) 
        d_G0 = self.getDG(phi0, kappa0, d_kappa_d_phi0, d_phi0)
        dd_G_dt2_0 = self.getDDG_DT2(self.phi0, theta0, kappa0, 1.0, d_G0, 0.0, 1.0, d_kappa_d_phi0)
        self.d_G0 = d_G0
        self.dd_G_dt2_0 = dd_G_dt2_0 
        #print '[phi0, d_phi_d_t0, H_init, z_init, self.getDOmegaDT(phi0, d_phi_d_t0, H_init, z_init), self.dOmegaDTOfFunction(phi0, d_phi_d_t0, H_init, z_init)] = ' + str([phi0, d_phi_d_t0, H_init, z_init, self.getDOmegaDT(phi0, d_phi_d_t0), self.dOmegaDTOfFunction(phi0, d_phi_d_t0)]) 
        #print '[H_init, phi0, initial_zs[0], omega0, d_omega0, d_phi_d_t0] = ' + str([H_init, phi0, initial_zs[0], omega0, d_omega_d_t0, d_phi_d_t0])
        d_H_d_t_init, d_H_init_terms = self.getDHDT(H_init, phi0, initial_zs[0], kappa0, d_kappa_d_phi0, theta0)
        d_H_init = d_H_d_t_init * self.d_t_to_d_z_scaling(z_init, H_init) 
        
        #d_theta0 = -3.0 * theta0 / self.a0 - 1.0 / (self.a0 * H_init) * self.dVOfFunction(phi0)
        #d_theta0 = - (( theta0 * self.dVOfFunction(phi0) * (-X0 + 3.0 * X0 ** 2.0) + 6.0 * H_init * self.VOfFunction(phi0) * (-X0 + 2.0 * X0 ** 2.0) )
        #             / ( self.VOfFunction(phi0) * theta0 * (-1.0 + 6.0 * X0) * self.a0 * H_init) )
        print '[A, nu, alpha, psi, eps, kappa0, phi0, theta0, d_phi0, d_theta_d_t0, d_theta0, H_init, d_H_d_t_init, d_H_init, G0, d_G0, dd_G_dt2_0, z_init, self.kappaOfFunction(phi0),  self.dKappaDPhiOfFunction(phi0)] = ' + str([A, nu, alpha, psi, eps, kappa0, phi0, theta0, d_phi0, d_theta_d_t0, d_theta0, H_init, d_H_d_t_init, d_H_init, G0, d_G0, dd_G_dt2_0, z_init, self.kappaOfFunction(phi0),  self.dKappaDPhiOfFunction(phi0)])

        self.states = []
        self.derivs = []
        self.used_zs = []
        self.d_H_terms = [] 
        
        self.runCalculationWithGivenzs(self.zs)
            
 

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
        omega = self.getOmega(phi) 
        #print 'type(X) = ' + str(type(X)) 
        #print 'X = ' + str(X)
        #print 'type(X) in [list]) = ' + str(type(X) in [list]) 
        if (type(omega) in [list]): omega = np.array(omega) 
        #print 'type(z) = ' + str(type(z)) 
        #print 'z = ' + str(z)
        x = theta / phi
        cosmo_omega = a ** (-3.0) * self.OmM0 + a ** (-4.0) * self.OmR0 + self.OmL0 
        #print '[self.OmM0, self.OmR0, self.OmL0] = ' + str([self.OmM0, self.OmR0, self.OmL0])
        H = - x / 2.0 + np.sqrt((x / 2.0) ** 2.0 + 4.0 * (omega / 6.0 * (x) ** 2.0 + cosmo_omega / phi ) )
        return H 
    
    def getOmega(self, phi):
        #print 'tau = ' + str(tau)
        #print 'type(tau) = ' + str(type(tau))
        #print "type(tau) in ['list',np.ndarray] = " + str(type(tau) in ['list',np.ndarray])
        if not(type(phi) in [list, np.ndarray]): phi= float(phi)
        return self.omegaOfFunction(phi)

    def getDOmegaDT(self, phi, d_phi_d_t):
        if not(type(phi) in [list, np.ndarray]): phi= float(phi)
        if not(type(d_phi_d_t) in [list, np.ndarray]): d_phi= float(d_phi_d_t)
        return self.dOmegaDTOfFunction(phi, d_phi_d_t)

    def getDPhiDT(self, z, H, omega, phi):
        #if H == 1.0: print 'At H == 1.0, [z, H, omega, phi] =' + str([z, H, omega, phi])
        #if H == 1.0: print 'At H == 1.0, [z, H, omega, phi] = ' + str([z, H, omega, phi])
        #print '((3.0 * phi * H / omega) ** 2.0 + 6.0 / omega * (H ** 2.0 * phi ** 2.0 - (self.OmM0 * (1.0 + z) ** 3.0 + self.OmR0 * (1.0 + z) ** 4.0 + self.OmL0) * phi )) = ' + str(((3.0 * phi * H / omega) ** 2.0 + 6.0 / omega * (H ** 2.0 * phi ** 2.0 - (self.OmM0 * (1.0 + z) ** 3.0 + self.OmR0 * (1.0 + z) ** 4.0 + self.OmL0) * phi )))

        d_phi_d_t = 3.0 * H * phi / omega +  np.sqrt((3.0 * phi * H / omega) ** 2.0 + 6.0 / omega * (H ** 2.0 * phi ** 2.0 - (self.OmM0 * (1.0 + z) ** 3.0 + self.OmR0 * (1.0 + z) ** 4.0 + self.OmL0) * phi ))
        #if H == 1.0: print 'At H == 1.0, [z, H, omega, phi, d_phi_d_t] =' + str([z, H, omega, phi, d_phi_d_t])
        return d_phi_d_t

    def getDH(self, H, phi, z, omega, d_omega_d_t, d_phi_d_t): 
        #d_t_to_d_z_scaling = - (1.0 + z) * H 
        d_H_term1 = - H ** 2.0
        d_H_term2 = - omega / 3.0 * d_phi_d_t ** 2.0 / phi ** 2.0
        d_H_term3 = H * d_phi_d_t / phi
        d_H_term4 = -( (4.0 * self.OmR0 * (1.0 + z) ** 4.0 + self.OmM0 * (1.0 + z) ** 3.0 - 2.0 * self.OmL0) * omega + 3.0 ) / ((2.0 * omega + 3.0) * phi)
        d_H_term5 = 0.5 * d_omega_d_t / (2.0 * omega + 3.0) * d_phi_d_t / phi 
        d_H_terms = [d_H_term1, d_H_term2, d_H_term3, d_H_term4, d_H_term5, d_H_term5, self.d_t_to_d_z_scaling(z, H)] 
        #self.d_H_terms = self.d_H_terms + [[d_H_term1, d_H_term2, d_H_term3, d_H_term4, d_H_term5, d_H_term5, d_H_term6, self.d_t_to_d_z_scaling(z, H)]]
        d_H_d_t = (d_H_term1 + d_H_term2 + d_H_term3 + d_H_term4 + d_H_term5 + d_H_term5)
        d_H = self.d_t_to_d_z_scaling(z, H) * d_H_d_t
        #if H == 1.0: print 'At H == 1.0, [d_H_d_t, d_H, self.d_t_to_d_z_scaling(z, H), phi, d_phi_d_t]  = ' + str([d_H_d_t, d_H, d_H_terms, phi, d_phi_d_t]) 
        return d_H, d_H_terms 
        
    def tPhiToGSystem(self, state, z):
        #t, tau, dL, phi, H = state
        t, tau, phi, H = state
        self.states = self.states + [state]
        #t, OmL = state
        #print '[z, t, OmL] = ' + str([z,t,OmL])
        omega = self.getOmega(phi) 
        #H = self.getH(a, theta, phi)
        #derivative in z 
        d_t = -1.0 * (-1.0 / (1.0 + z) * 1.0 / H)
        d_tau = 1.0 / H
        d_phi_d_t = self.getDPhiDT(z, H, omega, phi)
        d_phi = d_phi_d_t * self.d_t_to_d_z_scaling(z, H)
        d_omega_d_t = self.getDOmegaDT(phi, d_phi_d_t)
        #if H == 1.0:
        #    print '[phi, d_phi_d_t, H, z, self.getDOmegaDT(phi, d_phi_d_t)] = ' + str([phi, d_phi_d_t, H, z, self.getDOmegaDT(phi, d_phi_d_t)])
        #    print '[H, phi, z, omega, d_omega_d_t, d_phi_d_t] = ' + str([H, phi, z, omega, d_omega_d_t, d_phi_d_t])
        d_H, d_H_terms = self.getDH(H, phi, z, omega, d_omega_d_t, d_phi_d_t)
        self.d_H_terms = self.d_H_terms + [d_H_terms]
        
        #d_dL = (dL / (1.0 + z)) + ((1.0 + z) / H)
        derivs = [d_t, d_tau, d_phi, d_H]
        self.used_zs = self.used_zs + [z] 
        self.derivs = self.derivs + [derivs]
        return derivs 
        #print '[d_t, d_OmL] = ' + str([d_t, d_OmL]) 
        #return [d_t, d_OmL]

    def gettPhiToGSystem(self, init_state, a):
        return odeint(self.tPhiToGSystem, init_state, a)

    def runCalculationWithGivenzs(self, zs):
        #print 'zs = ' + str(zs)
        #if zs[0] > 0.0: zs = [0.0] + zs #Initial values are known at z = 0.0
        #self.zs = zs
        self.zs = zs 
        self.states = self.gettPhiToGSystem([self.t0, self.tau0, self.phi0, self.H_init], self.zs)
        #self.states = self.gettOmegaLumDistSystem([self.t0, self.OmL0], self.zs)
        self.tH0_of_z = [state[0] for state in self.states]
        self.tauH0_of_z = [state[1] for state in self.states]
        self.phi_of_z = [state[2] for state in self.states]
        self.H_of_z = [state[3] for state in self.states] 
        self.G_of_a = [1.0 / phi for phi in self.phi_of_z]
        self.theta_of_z = [state[3] for state in self.states]
        self.omega_of_z = [self.omegaOfFunction(phi) for phi in self.phi_of_z] 
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
    
    def getUnitfulldLs(self, target_zs = None):
        #print 'Units on dL are ' + self.dL_units
        if not(target_zs is None) or not(self.dL_of_z is None):
            unitless_dLs = self.getUnitlessdLs(target_zs = target_zs) 
        return [dL * self.dL_scaling for dL in self.dL_of_z]

    def getMus(self, target_zs = None):
        #print 'dL_units are in ' + self.dL_units + ".  They should be mega_parsec."
        return 25.0 + 5.0 * np.log10(self.getUnitfulldLs(target_zs = target_zs)) #ignore first dL, as it is 0 (today) 
                                    
    #You should give the X function as a function of t, tau, and z (in that order).
    # If you only want it to be a function of only one or two of them, just make
    # the output independent of the variables that are not of interest.  
    def __init__(self, A, nu, alpha, psi, initial_zs, 
                       t0 = None, tau0 = None, B = None, phi0 = None, eps = None, 
                       canon_mus = None, OmM0 = None, OmL0 = None, OmR0 = None, Om0 = None, H0 = None, 
                       H0_units = 'mega_years', dL_units = 'mega_parsec'):
        #Monodromic DE has a specific potential.  So that is hard coded here; not something that a use gives
        omegaOfSinglePhiFunction = lambda phi, B, phi0, alpha, nu, A, psi: (2.0 * B * abs(1 + eps - phi /phi0) ** (-alpha) * (1 - A * np.sin(nu * phi / phi0 + psi)) - 2.0 ) / 3.0
        #omegaOfSinglePhiFunction = lambda phi, B, phi0, alpha, nu, A, psi: (2.0 * B * (1 + eps - phi /phi0) ** (-alpha) * (1 - A * np.sin(nu * phi / phi0 + psi)) - 2.0 ) / 3.0
        omegaOfFunction = lambda phis, B, phi0, alpha, nu, A, psi: [omegaOfSinglePhiFunction(phi, B, phi0, alpha, nu, A, psi) for phi in phis] if type(phis) in [list,np.ndarray] else omegaOfSinglePhiFunction (phis, B, phi0, alpha, nu, A, psi)

        self.d_t_to_d_z_scaling = lambda z, H: - (1.0 + z) * H 
        dOmegaDTOfSinglePhiFunction = lambda phi, d_phi_d_t, B, phi0, alpha, nu, A, psi: ((2.0 * B * (-alpha) * (-1.0 / phi0) * abs(1 + eps - phi /phi0) ** (-alpha - 1.0) * (1 - A * np.sin(nu * phi / phi0 + psi))) / 3.0 +  (2.0 * B * abs(1 + eps - phi /phi0) ** (-alpha) * (- A * nu / phi0 * np.cos(nu * phi / phi0 + psi))) / 3.0 ) * d_phi_d_t # * self.d_t_to_d_z_scaling(z, H)
        dOmegaDTOfFunction = lambda phis, d_phi_d_ts, B, phi0, alpha, nu, A, psi: ([dOmegaDTOfSinglePhiFunction(phis[i], d_phi_d_ts[i], B, phi0, alpha, nu, A, psi) for i in range(len(phis))] if (type(phis) in [list,np.ndarray] and type(d_phi_d_ts) in [list,np.ndarray])
                                                                             else [dOmegaDTOfSinglePhiFunction(phis[i], d_phi_d_ts, B, phi0, alpha, nu, A, psi) for i in range(len(phis))] if (type(phis) in [list,np.ndarray])
                                                                             else [dOmegaDTOfSinglePhiFunction(phis, d_phi_d_ts[i], B, phi0, alpha, nu, A, psi) for i in range(len(d_phi_d_ts))] if (type(d_phi_d_ts) in [list,np.ndarray])
                                                                             else dOmegaDTOfSinglePhiFunction (phis, d_phi_d_ts, B, phi0, alpha, nu, A, psi) ) 
    
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
        if B is None:
            B = 1.0
        self.B = B
        if eps is None:
            eps = 10.0 ** -4.0
        H_init = 1.0
        self.H_init = H_init 

        print '[self.B, self.phi0, self.alpha, self.nu, self.A, self.psi] = ' + str([self.B, self.phi0, self.alpha, self.nu, self.A, self.psi])
        #print 'omegaOfFunction([0.9, 1.0], self.B, self.phi0, self.alpha, self.nu, self.A, self.psi) = ' + str(omegaOfFunction([0.9, 1.0], self.B, self.phi0, self.alpha, self.nu, self.A, self.psi))
        self.omegaOfFunction = lambda phis: omegaOfFunction(phis, self.B, self.phi0, self.alpha, self.nu, self.A, self.psi)
        self.dOmegaDTOfFunction = lambda phis, d_phi_d_ts: dOmegaDTOfFunction(phis, d_phi_d_ts, self.B, self.phi0, self.alpha, self.nu, self.A, self.psi) 
        #self.dVOfFunction = lambda phis: dVOfFunction(phis, self.B, self.phi0, self.alpha, self.nu, self.A)
        z_init = initial_zs[0] 
        
        #V_init = self.VOfFunction(phi0)
        omega_init = self.omegaOfFunction(phi0)
        omega0 = self.omegaOfFunction(phi0)
        d_phi_d_t0 = self.getDPhiDT(initial_zs[0], H_init, omega0, phi0)
        d_phi0 = d_phi_d_t0 * self.d_t_to_d_z_scaling(z_init, H_init)
        
        d_omega_d_t0 = self.dOmegaDTOfFunction(phi0, d_phi_d_t0)
        #print '[phi0, d_phi_d_t0, H_init, z_init, self.getDOmegaDT(phi0, d_phi_d_t0, H_init, z_init), self.dOmegaDTOfFunction(phi0, d_phi_d_t0, H_init, z_init)] = ' + str([phi0, d_phi_d_t0, H_init, z_init, self.getDOmegaDT(phi0, d_phi_d_t0), self.dOmegaDTOfFunction(phi0, d_phi_d_t0)]) 
        #print '[H_init, phi0, initial_zs[0], omega0, d_omega0, d_phi_d_t0] = ' + str([H_init, phi0, initial_zs[0], omega0, d_omega_d_t0, d_phi_d_t0])
        d_H_init, d_H_init_terms = self.getDH(H_init, phi0, initial_zs[0], omega0, d_omega_d_t0, d_phi_d_t0)
        
        #d_theta0 = -3.0 * theta0 / self.a0 - 1.0 / (self.a0 * H_init) * self.dVOfFunction(phi0)
        #d_theta0 = - (( theta0 * self.dVOfFunction(phi0) * (-X0 + 3.0 * X0 ** 2.0) + 6.0 * H_init * self.VOfFunction(phi0) * (-X0 + 2.0 * X0 ** 2.0) )
        #             / ( self.VOfFunction(phi0) * theta0 * (-1.0 + 6.0 * X0) * self.a0 * H_init) )
        print '[A, nu, alpha, psi, eps, phi0, d_phi_d_t0, d_phi0, H_init, d_H_init, z_init, self.omegaOfFunction(phi0),  self.dOmegaDTOfFunction(phi0, d_phi0)] = ' + str([A, nu, alpha, psi, eps, phi0, d_phi_d_t0, d_phi0, H_init, d_H_init, z_init, self.omegaOfFunction(phi0),  self.dOmegaDTOfFunction(phi0, d_phi0)])

        self.states = []
        self.derivs = []
        self.used_zs = []
        self.d_H_terms = [] 
        
        self.runCalculationWithGivenzs(self.zs)
        self.dL_of_z = None 
            
 

import numpy as np
from scipy.integrate import odeint
from scipy.interpolate import interp1d

class RingermacherModel:

    def getPhiDot(self, theta):
        phi_dot = theta

        return phi_dot

    def getDThetaDT(self, theta, phi, z):
        Omphi = self.Omphi_funct(theta, phi)
        theta_dot = -3.0 * self.getH(z, theta, phi) * theta - self.m * 2.0 * phi

        return theta_dot

    def getH(self, z, theta, phi):
        Omphi = self.Omphi_funct(theta, phi)
        H = np.sqrt(self.OmL + self.OmB * (1+z) ** 3.0 + Omphi)

        return H

    def ringermacherSystem(self, state, z):
        #t, tau, dL, phi, H = state
        t, phi, theta, dL = state
        H = self.getH(z, theta, phi)
        self.states = self.states + [state]
        d_t = -1.0 * (-1.0 / (1.0 + z) * 1.0 / H)
        d_phi_d_t = theta
        d_phi = d_phi_d_t * self.d_t_to_d_z_scaling(z, H)
        d_theta_d_t = self.getDThetaDT(theta, phi, z)
        d_theta = d_theta_d_t * self.d_t_to_d_z_scaling(z, H)
        #d_kappa_d_phi = self.getDKappaDPhi(phi)
        #d_theta_d_t = self.getDThetaDT(z, H, kappa, d_kappa_d_phi, theta)
        #d_theta = d_theta_d_t * self.d_t_to_d_z_scaling(z, H)

        d_dL = (dL / (1.0 + z)) + ((1.0 + z) / H)
        derivs = [d_t, d_phi, d_theta, d_dL]
        self.used_zs = self.used_zs + [z]
        self.derivs = self.derivs + [derivs]
        return derivs

    def getRingermacherSystem(self, init_state, z):
        return odeint(self.ringermacherSystem, init_state, z)

    def runCalculationWithGivenzs(self, zs):
       #print 'zs = ' + str(zs)
       #if zs[0] > 0.0: zs = [0.0] + zs #Initial values are known at z = 0.0
       #self.zs = zs
       self.zs = zs
       self.states = self.getRingermacherSystem([self.t0, self.phi0, self.theta0, self.dL_init], self.zs)
       #self.states = self.gettOmegaLumDistSystem([self.t0, self.OmL0], self.zs)
       self.tH0_of_z = [state[0] for state in self.states]
       self.phi_of_z = [state[1] for state in self.states]
       self.d_phi_of_z = [state[2] for state in self.states]
       self.dL_of_z = [state[3] for state in self.states]

       self.a_of_z = 1.0 / ((np.array(self.zs)) + 1.0)
       self.H_of_z = [self.getH(self.zs[i], self.d_phi_of_z[i], self.phi_of_z[i] ) for i in range(len(self.zs))]
       #print ('self.H_of_z = ' + str(self.H_of_z))
       self.H_of_z_interp = interp1d (self.zs, self.H_of_z)

    def __init__(self, OmL = 0.735, OmB = 0.043, H0 = 100.0, fphi = 3.5, t0 = None, phi0 = None, dL_init = None):

        self.OmL = OmL
        self.OmB = OmB
        self.H0 = H0
        self.fphi = fphi
        self.m = 2.0 * np.pi * self.fphi

        if t0 is None:
            t0 = 0.0
        self.t0 = t0
        if phi0 is None:
            phi0 = 0.0
        self.phi0 = phi0
        Omphi_0 = 1.0 - self.OmL - self.OmB
        print ('Omphi_0 = ' + str(Omphi_0))
        self.theta0 = np.sqrt(Omphi_0 * 2.0 - self.m ** 2.0 * self.phi0 ** 2.0)
        print ('self.theta0 = ' + str(self.theta0))
        if dL_init is None:
            dL_init = 0.0
        self.dL_init = 0.0

        #self.HFunt = lambda a_s, Omphi_local: self.H0 * np.sqrt(self.OmL + self.Omb * a_s ** -3.0 + Omphi_local)
        self.Omphi_funct = lambda theta, phi: 0.5 * theta ** 2.0 + 0.5 * self.m ** 2 * phi ** 2.0
        print ('self.Omphi_funct(self.theta0, self.phi0) = ' + str(self.Omphi_funct(self.theta0, self.phi0)))
        self.d_t_to_d_z_scaling = lambda z, H: 1.0 / (- (1.0 + z) * H ) # d_d_z = self.d_t_to_d_z_scaling * d_d_t

        self.states = []
        self.derivs = []
        self.used_zs = []
        self.d_H_terms = []

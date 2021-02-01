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

class ResidualMuCalculatorForOffCenterVoid:

    def VoidObservableSystem(init_state, z):
        return odeint(VoidObservableSystem, self.init_state, z)

    def VoidObservableSystem(self, state, z):
        t, r, theta, p = state
        init_conditions = [r0, t0, 0.0, p0]
        A = self.AOfFunction(t, r)
        Ar = self.ArOfFunction(t, r)
        At = self.AtOfFunction(t, r)
        Arr = self.ArrOfFunction(t, r)
        Atr = self.AtrOfFunction(t, r)
        k = self.kOfFunction(r)
        kr = self.krOfFunction(r)
        q = (Ar * Atr / (1.0 - k) * p ** 2.0 + At * self.single_J ** 2.0 / A ** 3.0)
        dt_dz = -(1.0 + z) / q
        dr_dz = p / q
        dtheta_dz = J / (q * A ** 2.0)
        dp_dz = 1.0 / q * ((1.0 - k)/Ar * J ** 2.0 / A ** 3.0 + 2.0 * Atr / Ar * p * (1.0 + z) - p ** 2.0 * (Arr / Ar + kr / (2.0 - 2.0 * k)))

        return [dt_dz, dr_dz, dtheta_dz, dp_dz ]

    def runCalculationWithGivenzsAndAngle(self, zs, angle):
        if zs[0] > 0.0: zs = [0.0] + zs #Initial values are known at z = 0.0
        self.zs = zs
        self.single_J = self.JOfFunction(angle)
        p0 = self.p0OfFunction(angle)
        self.states = self.gettOmegaLumDistSystem([self.t0, self.r0, self.theta0, p0], self.zs)

    def getMus(self):
        #print 'dL_units are in ' + self.dL_units + ".  They should be mega_parsec."
        return 25.0 + 5.0 * np.log10(self.getUnitfulldLs()[1:]) #ignore first dL, as it is 0 (today)


    #Referencing this paper: https://iopscience-iop-org.ezp-prod1.hul.harvard.edu/article/10.1088/1475-7516/2010/05/006/pdf
    #We assume Om(r) = Om_bg + \Delta Om * exp(-r/rs) ** 2.0
    #We assume H0(r) = 3 * H0 / 2.0 * (1/OmK(r) - OmM(r) / sqrt(OmK(r) ** 3.0 ) sinh-1(sqrt(OmK(r)  / OmM(r))))
    # and we know OmM(r) + OmK(r) = 1

    def getOmM(self, r):
        OmM = self.OmMOut + (self.DeltaOmM) * np.exp(-(r / self.rs) ** 2.0)
        return OmM

    def getH0(self, OmM):
        H0 = 3.0 * self.H0 / 2.0 * (1.0 / (1.0 - OmM) - OmM / (1.0 - OmM) ** (3.0 / 2.0) * np.arcsinh(np.sqrt((1.0 - OmM) / OmM)))
        return H0

    def getEta(self, t, r, OmM = None):
        if OmM is None: OmM = selfOmM(r)
        sinhEtaMEta = self.getH0(OmM) * t * (2.0 * (1.0 - OmM) ** (3.0 / 2.0)) / OmM
        eta = self.sinhEtaMEtaToEtaInterp(sinhEtaMEta)
        return eta

    def getA0(self, r):
        return 1

    def getA(self, t, r):
        OmM = self.OmM(r)
        A0 = self.getA0(r)
        eta = getEta(t, r, OmM = OmM)
        A = OmM / (2.0 * (1.0 - OmM)) * (np.cosh(eta) - 1.0) *  A0
        return A



    def getA(O)

    def __init__(self,
                 OmOut = 0.3, DeltaOmM = 0.0,  rs = 1.0,
                 etas_to_solve = np.linspace(0.0, 100.0, 1001),
                 self.OmMFunct = lambda ts, rs: 0.3 + ts * 0.0,
                 AOfFunction = lambda ts, rs: 1.0 if type(ts) in [float, int] else -1.0 + np.zeros(np.shape(ts)),
                 ArOfFunction = lambda ts, rs: 1.0 if type(ts) in [float, int] else -1.0 + np.zeros(np.shape(ts)),
                 ArrOfFunction = lambda ts, rs: 1.0 if type(ts) in [float, int] else -1.0 + np.zeros(np.shape(ts)),
                 AtOfFunction = lambda ts, rs: 1.0 if type(ts) in [float, int] else -1.0 + np.zeros(np.shape(ts)),
                 AtrOfFunction = lambda ts, rs: 1.0 if type(ts) in [float, int] else -1.0 + np.zeros(np.shape(ts)),
                 kOfFunction = lambda rs:
                 krOfFunction = lambda rs:
                 r0 = 0.0, t0 = 0.0, incident_angles = 0.0,
                 canon_mus = None, a0 = None, OmM0 = None, OmL0 = None, OmR0 = None, Om0 = None, H0 = None, initial_zs = None, H0_units = 'mega_years', dL_units = 'mega_parsec'):

    self.sinhEtaMEtaToEtaInterp = scipy.interp1d([np.sinh(eta) - eta for eta in etas_to_solve], etas_to_solve])
    self.OmOut = OmMOut
    self.DeltaOmM = DeltaOmM
    self.rs = rs
    self.AOfFunction = lambda ts, rs, eta: alphaFunction(ts, rs) / 2.0 * betaFunction(ts, rs) * (np.cosh(eta) - 1.0) +
    self.AOfFunction = AOfFunction
    self.ArOfFunction = ArOfFunction
    self.AtOfFunction = AtOfFunction
    self.AtrOfFunction = ArtOfFunction
    self.kOfFunction = kOfFunction
    self.incident_angle = incident_angle
    self.r0 = r0
    self.t0 = t0
    self.theta0 = 0.0
    #p0s = [np.sqrt(1.0 - self.kOfFunction(r0)) / (self.ArOfFunction(r0, t0)) * np.cos(incident_angle) for incident_angle in incident_angles]
    #self.p0 = p0
    self.p0OfFunction =  lambda angle: np.sqrt(1.0 - self.kOfFunction(r0)) / (self.ArOfFunction(t0, r0)) * np.cos(angle)
    self.JOfFunction = lambda angle: self.AOfFunction(self.t0, self.r0) * np.sin(angle)
    #self.Js = [self.AOfFunction(t0, r0) * np.sin(incident_angle) for incident_angle in incident_angles]

    astro_archive = AstronomicalParameterArchive()
    cosmo_archive = CosmologicalParameterArchive()
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
    self.A0 = ROfFunction(self.t0, self.r0)

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

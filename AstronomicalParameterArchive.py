#A set of general astronomical parameters

import math
import numpy as np

class AstronomicalParameterArchive:

    def getSteridanToSquareArcsec(self):
        return 4.255 * 10.0 ** 10.0

    def getKmPerSToPcPerYr(self):
         km_per_s_to_pc_per_yr = 1.0 / ( self.getParsecToM() / 1000.0  * self.getSecondToYear() )
         return km_per_s_to_pc_per_yr

    #gives Newton's constant in m^3 / (kg s^2)
    def getGravitationalConstant(self):
        return self.gravitational_constant

    #gives you mass of the sun, in kg
    def getSolarMass(self):
        return self.solar_mass

    #gives you mass of Milky Way
    def getMWMass(self):
        return self.MWMass

    #return intensity at earth of a stars of magnitude mags in units of kW/m^2
    def convertMagnitudesToIntensities(self,mags):
        return 10**(np.array(mags)/(-2.5))*(self.Vega_solar_lums*self.solar_lum_to_kw)/(4.0*math.pi*(self.Vega_dist_pars*self.parsec_to_m)**2)

    def getDegToRad(self):
        return self.deg_to_rad

    def getDegToAmin(self):
        return self.deg_to_amin

    def getAminToAsec(self):
        return self.amin_to_asec

    def getDegToAsec(self):
        return self.deg_to_amin * self.amin_to_asec

    def getVegaSolarLums(self):
        return self.Vega_solar_lums

    def getSolarLumToKw(self):
        return self.solar_lum_to_kw

    def getVegaDistPars(self):
        return self.Vega_dist_pars

    def getParsecToM(self):
        return self.parsec_to_m

    def getSecondToYear(self):
        return self.sec_to_yr

    def getGamma(self):
        return self.gamma

    def getElectronCharge(self):
        return self.e_charge


    def getPrefix(self, prefix):
        return self.metric_prefixes[prefix.lower()]

    def getAngularMotionConversionFactor(self):
        return self.angleMotConv

    def getPlancksConstant(self):
        return self.planck

    def getc(self):
        return self.c

    #Return RA, Dec values that cover whole sky in form of [min(RA), max(RA), min(Dec), max(Dec)] in degrees
    def getFullSky(self):
        return self.skyDegreeRange

    def __init__(self):
        self.gravitational_constant = 6.674 * 10.0 ** -11.0 #Newton's constant in m^3 / (kg s^2)
        self.solar_mass = 1.98847 * 10.0 ** 30.0 #solar mass in kg
        self.skyDegreeRange = [0.0, 360.0, -90.0, 90.0]
        self.MWMass = 1.0 * 10**12 #Solar mass of Milky Way (Wikipedia 0.8 - 1.5 * 10**12)
        self.c = 299792458.0 * 10**-3 #speed of light in km / s
        self.planck = 6.62607004 * 10. ** -34.0 #planck's constant (h, not hbar) in J s
        self.Vega_solar_lums = 40.12 #solar luminosity of Vega
        self.solar_lum_to_kw = 3.828 * 10**23 #solar luminosity in kilowatts
        self.Vega_dist_pars = 7.68 #distance to Vega in parsecs
        self.parsec_to_m = 3.086 * 10**16 #number of meters in a parsec
        self.deg_to_rad = math.pi/180.0 #convert degrees to radians
        self.deg_to_amin = 60.0 #convert degress to arcmin
        self.amin_to_asec = 60.0 #convert arcmin to arcsecond
        self.gamma=4.307 * 10**(-3) #constant scaling between (G M)/ (rs sigsqr) when M, rs, sigsqr are expressed in particular units (here, )
        self.angleMotConv = 0.0829 / 1000 # constant converting angular proper motion of galaxy into correction to los velocity, provided distance to galaxy is given in pc, position of star on sky is in degrees, and measured proper motion is in mas/yr
        self.e_charge = 1.60217662 * 10 ** -19 #electron charge in Coloumbs

        #Do calculation in Julian years
        self.sec_to_yr = 1.0 / 31557600.0

        self.metric_prefixes = {'exa'   : 10.0 ** 18.0, 'peta' : 10.0 ** 15.0, 'tera'  : 10.0 ** 12.0,  'giga'  : 10.0 ** 9.0,
                                'mega'  : 10.0 ** 6.0,  'kilo' : 10.0 ** 3.0,  'hecto' : 10.0 ** 2.0,   'deca'  : 10.0 ** 1.0,
                                'unity' : 1.0,          'deci' : 10.0 ** -1.0, 'centi' : 10.0 ** -2.0,  'milli' : 10.0 ** -3.0,
                                'micro' : 10.0 ** -6.0, 'nano' : 10.0 ** -9.0, 'pico'  : 10.0 ** -12.0, 'femto' : 10.0 ** -15.0,
                                'atto'  : 10.0 ** -18.0 }

#Stores canonical cosmological values for other programs to reference
#Allows for user to specify other cosmological values to testing how, eg, results could change if Planck data is incomplete
from AstronomicalParameterArchive import AstronomicalParameterArchive
import math

class CosmologicalParameterArchive:

    def gett0(self):
        return [self.t0, self.t0Err]

    def getH0(self, units = 'normal'):
        units = units.lower()
        if units in ['normal','standard', 'canon','usual','stand','can','norm']:
            return [self.H0, self.H0Err]
        else:
            km_to_Mpc = 1.0 / self.astro_arch.getParsecToM() * self.astro_arch.getPrefix('kilo') * 1.0 / self.astro_arch.getPrefix('mega')
            if units in ['s','second','seconds','sec','secs']:
                return [self.H0 * km_to_Mpc, self.H0Err * km_to_Mpc]
            elif units in ['yr','years','yrs','year']:
                s_to_yr = self.astro_arch.getSecondToYear()
                return [self.H0 * km_to_Mpc * 1.0 / s_to_yr, self.H0Err * km_to_Mpc * 1.0 / s_to_yr]
            else:
                'Units for H0 not recognized.  Returning default of 1/s. '
                return [self.H0 * km_to_Mpc, self.H0Err * km_to_Mpc]

    #return H0 in inverse seconds
    def getH0_invSec(self):
        return [self.H0 * (1.0 / (3.0857 * 10 ** 19)), self.H0Err * (1.0 / (3.0857 * 10 ** 19))]

    def getc(self):
        return self.c

    def geth(self):
        return [self.h,self.herr]

    def getOmegaM(self):
        return [self.OmegaM,self.OmegaMErr]

    def getOmegaLambda(self):
        return [self.OmegaLambda,self.OmegaLambdaErr]

    def getOmega0(self):
        return [self.Omega0,self.Omega0Err]

    def getOmegaR(self):
        return [self.OmegaR,self.OmegaRErr]

    def getAgeOfUniverse(self, units = 's'):
        if units.lower() in ['s', 'sec', 'seconds']:
            return [self.t0, self.t0Err]
        elif units.lower() in ['y','yr','yrs','years','year']:
            return [self.t0 * 1.0 / self.yrsToSec, self.t0Err * 1.0 / self.yrsToSec]
        else:
            return [self.t0, self.t0Err]

    #Unless otherwise specified, measured cosmological values are taken from Planck (https://arxiv.org/pdf/1502.01589.pdf)
    def __init__(self, H0 = None, H0Err = None, OmegaM = None, OmegaMErr = None, OmegaLambda = None, OmegaLambdaErr= None, Omega0 = None, Omega0Err = None, OmegaR = None, OmegaRErr = None, params_source = 'planck'):

        self.astro_arch = AstronomicalParameterArchive()
        self.params_source = params_source
        param_index_array = {'H0':0, 'H0E':1, 'OmM':2, 'OmME':3, 'OmL':4, 'OmLE':5, 'Om0':6, 'Om0E':7, 'OmR':8, 'OmRE':9}
        specified_params = [H0, H0Err, OmegaM, OmegaMErr, OmegaLambda, OmegaLambdaErr, Omega0, Omega0Err, OmegaR, OmegaRErr]

        planck_zeq = 3393.0
        planck_zeqE = 49.0

        planck_params =     [67.31,  0.96,     0.315,  0.013,    0.685,   0.013,    1.0,     0.0,      0.0, 0.0]

        planck_params[param_index_array['OmR']] = planck_params[param_index_array['OmM']] / (1.0 + planck_zeq)
        planck_params[param_index_array['OmRE']] = math.sqrt( (planck_params[param_index_array['OmME']] / (1.0 + planck_zeq)) ** 2.0
                                                              + ((planck_params[param_index_array['OmM']] * planck_zeqE) / ((1.0 + planck_zeq) ** 2.0) ) ** 2.0 )
        #In Planck, universe is assumed totally flat with (basically) off of energy density in matter and lambda.  So they are determined by each other
        planck_params[param_index_array['OmL']] = planck_params[param_index_array['Om0']] - planck_params[param_index_array['OmM']] - planck_params[param_index_array['OmR'] ]

        planck_params[param_index_array['OmLE']] = math.sqrt(planck_params[param_index_array['OmME']] ** 2.0 + planck_params[param_index_array['OmRE']] ** 2.0)


        #Pantheon offers no constraints on H0
        pantheon_params  = [70.0, 0.0, 0.284, 0.013, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]
        pantheon_params[param_index_array['OmR']] = planck_params[param_index_array['OmR']]
        pantheon_params[param_index_array['OmRE']] = planck_params[param_index_array['OmRE']]
        pantheon_params[param_index_array['OmL']] = pantheon_params[param_index_array['Om0']] - pantheon_params[param_index_array['OmM']] - pantheon_params[param_index_array['OmR'] ]
        pantheon_params[param_index_array['OmLE']] = math.sqrt(pantheon_params[param_index_array['OmME']] ** 2.0 + pantheon_params[param_index_array['OmRE']] ** 2.0)

        if self.params_source.lower() == 'pantheon':
            params_to_use = pantheon_params
        elif self.params_source.lower()  == 'planck':
            params_to_use = planck_params
        else:
            params_to_use = planck_params

        params_to_use = [specified_params[i] if specified_params[i] != None else params_to_use[i] for i in range(len(param_index_array)) ]

        #H0 units: km/s/Mpc
        self.H0 = params_to_use[param_index_array['H0']]
        self.H0Err = params_to_use[param_index_array['H0E']]

        #OmegaM unitless
        self.OmegaM = params_to_use[param_index_array['OmM']]
        self.OmegaMErr = params_to_use[param_index_array['OmME']]

        #OmegaLambda unitless
        self.OmegaLambda = params_to_use[param_index_array['OmL']]
        self.OmegaLambdaErr = params_to_use[param_index_array['OmLE']]

        #Omega0 unitless
        self.Omega0 = params_to_use[param_index_array['Om0']]
        self.Omega0Err = params_to_use[param_index_array['Om0E']]

        #OmegaR unitless
        self.OmegaR = params_to_use[param_index_array['OmR']]
        self.OmegaRErr = params_to_use[param_index_array['OmRE']]

        # h = H0 / 100.0 (km/s/Mpc)
        self.h_to_H0_scaling = 100.0
        self.h = self.H0 / self.h_to_H0_scaling
        self.herr = self.H0Err / self.h_to_H0_scaling

        #speed of light in km / s
        self.c = 299792458.0 * 10**-3

        #Number of seconds in Julian Year
        self.yrsToSec = 31557600.0

        #Age of Universe in s
        self.t0 = 13.799 * 10.0 ** 9.0 * self.yrsToSec
        self.t0Err = 0.021 * 10.0 ** 9.0 * self.yrsToSec

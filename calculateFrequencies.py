from CosmologicalParameterArchive import CosmologicalParameterArchive
from AstronomicalParameterArchive import AstronomicalParameterArchive
import numpy as np
from loadSN import loadSN
import scipy.optimize as optimize
import scipy.integrate as integrate
import math
from CosmologicalCalculator import CosmologicalCalculator 

#min_length_scale given in pc 
def calcFrequenciesFromMinSpatialDistance(data_set = 1, data_type = 'real', freq_density_scaling = 1.0, start_freq_scaling = 1.0, min_length_scale = 1.0 * 10.0 ** 6.0, use_current_scale = 1, t_or_tau = 't', max_oscillations = None, use_max_oscillations = False):
    cosmo_arch = CosmologicalParameterArchive()
    astro_arch = AstronomicalParameterArchive()

    s_to_yr = astro_arch.getSecondToYear()
    age_of_univ = cosmo_arch.getAgeOfUniverse(units= 'yr')[0]
    
    all_sn = loadSN(data_set, data_type = data_type)
    all_zs = [sn['z'] for sn in all_sn]
    if t_or_tau.lower() in ['tau', 'taus']:
        all_taus = [sn['tau'] for sn in all_sn]
        all_lookBacks = [ s_to_yr * 10.0 ** (-6.0) * tau for tau in all_taus]
    else: 
        all_ts = [sn['t'] for sn in all_sn]
        all_lookBacks = [ s_to_yr * 10.0 ** (-6.0) * t  for t in all_ts]
    
    interval = max(all_lookBacks) - min(all_lookBacks)
    fundamental_min_freq = 2.0 * math.pi / interval
    min_freq = fundamental_min_freq * start_freq_scaling 

    max_z = max(all_zs)
    pc_to_m = astro_arch.getParsecToM()
    pc_to_km = pc_to_m * 10 ** -3.0
    c = cosmo_arch.getc() #c in km / s
    H0 = cosmo_arch.getH0_invSec()[0]
    OmM = cosmo_arch.getOmegaM()[0]
    OmL = cosmo_arch.getOmegaLambda()[0]
    OmR = cosmo_arch.getOmegaR()[0]
    H_of_z = lambda zs: np.sqrt( OmM * (1.0 + zs) ** 3.0 + OmL + OmR * (1.0 + zs) ** 4.0 )
    scaling = H0/c * 10 ** -3.0 #units are 1/m
    adjacent_max_z = optimize.brenth(lambda z: min_length_scale * pc_to_m * scaling
                                     - integrate.quad(lambda zs: 1.0 / H_of_z(zs), z, max_z)[0], 0.0, max_z)
    cosmo_calc= CosmologicalCalculator(zs = [adjacent_max_z, max_z])
    min_sep_ts = cosmo_calc.ts
    if use_current_scale:
        min_t_sep = (min_length_scale * pc_to_km)/ c * s_to_yr * 10.0 ** (-6.0)
    else:
        min_t_sep = abs(min_sep_ts[1] - min_sep_ts[0]) * s_to_yr * 10 ** (-6.0) #minimum observable interval, in Myrs

    if (use_max_oscillations and not (max_oscillations is None)):
        max_freq = max_oscillations * fundamental_min_freq
    else:
        max_freq = 2.0 * math.pi / min_t_sep 

    freq_step = fundamental_min_freq / freq_density_scaling
   
    frequencies = np.arange(min_freq, max_freq, freq_step).tolist()
    print 'min_freq = ' + str(min_freq) + '1/Myr'
    print 'max_freq = ' + str(max_freq) + '1/Myr'
    print 'len(frequencies) = ' + str(len(frequencies)) 

    return frequencies 
    
    

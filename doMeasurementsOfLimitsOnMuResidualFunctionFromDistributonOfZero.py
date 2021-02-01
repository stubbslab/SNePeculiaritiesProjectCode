import math
import time 
import numpy as np
from measureLimitsOfFunctionGivenDistributionOfZero import measureLimitsOfFunctionGivenDistributionOfZero
from loadSN import loadSN 
from AstronomicalParameterArchive import AstronomicalParameterArchive
from CosmologicalParameterArchive import CosmologicalParameterArchive

if __name__ == '__main__':
    #Number of sigma
    limiting_prob = 0.999999999652544
    
    #Must be a two parameter function
    # (ideally one is like a frequency and one is like an amplitude)

    #Start with a simple oscillitory function
    omega0 = 0.05
    freq_index_range_for_finding_peak = [-300, 300]
    bounds = [0.0, 10.0]
    muOfT = lambda ts, A, omega, phi: A * omega0 / omega * np.sin(np.array(ts)*omega + phi)

    #Function to limit should be of three parameters: lookback time, param that will be minimized over, extra param to vary
    function_to_limit = lambda ts, Amp_like, freq_like: muOfT(ts, Amp_like, freq_like, 0.0)
    
    file_of_zero_limits = '/Users/sasha/Documents/Harvard/physics/stubbs/SNIsotropyProject/randomBestFitParams/params_from_300_frequency_scaling50_10_sampling_of_fitted_taualt_gaussian_to_zero.csv'

    results_for_minimization = np.genfromtxt(file_of_zero_limits, delimiter = ',').transpose()
    frequencies = results_for_minimization[0]

    freq_like_params_on_which_to_minimize = [frequencies[5 + i * 500] for i in range(len(frequencies) / 500)]
    freq_indeces_for_finding_peak = [[max(0, 5 + i * 500 + freq_index_range_for_finding_peak[0]),
                                      min(len(frequencies), 5 + i * 500 + freq_index_range_for_finding_peak[1])]
                                     for i in range(len(frequencies) / 500)]
    #freq_indeces_for_finding_peak = [None for i in range(len(frequencies) / 500)]
    init_guesses = [0.5 * freq_like_param / omega0 for freq_like_param in freq_like_params_on_which_to_minimize]

    astro_arch = AstronomicalParameterArchive()
    s_to_yr = astro_arch.getSecondToYear() 
    cosmo_archive = CosmologicalParameterArchive()
    age_of_universe = cosmo_archive.getAgeOfUniverse( units = 'yr')[0] 
    all_sn = loadSN(1, pull_extinctions = 0)
    
    all_surveys = np.unique([sn['survey'] for sn in all_sn])
    all_ts = [sn['t'] for sn in all_sn]
    all_lookBacks = [(-6.0) * (t * s_to_yr) for t in all_ts]
    all_muErrs = [sn['muErr'] for sn in all_sn]
    all_muResids = [sn['muDiff'] for sn in all_sn]
    all_zerodMuResids = [sn['muDiff'] - sn['muDiffWMean'] for sn in all_sn]

    minimized_values = []
    for i in range(len(freq_like_params_on_which_to_minimize)):
        start = time.time() 
        print 'Locating limiting case for ' + str(i) + 'th frequency of ' + str(len(freq_like_params_on_which_to_minimize))  
        freq_like_param = freq_like_params_on_which_to_minimize[i]
        init_guess = init_guesses[i]
        function_to_limit_at_specific_freq_val = lambda ts, Amp_like: function_to_limit (ts, Amp_like, freq_like_param)
        min_results = measureLimitsOfFunctionGivenDistributionOfZero(all_lookBacks, all_muErrs, function_to_limit_at_specific_freq_val, init_guess, limiting_prob,
                                                                     file_of_zero_limits, 'no saving ability yet', bounds, 
                                                                     freq_width_to_calc = 100, fitting_funct_type = 'alt_norm', 
                                                                     freq_bounds_over_which_to_find_peak = freq_indeces_for_finding_peak[i],
                                                                     freq_bounds_over_which_to_check_peak = freq_indeces_for_finding_peak[i])
        minimized_values = minimized_values + [min_results['x']]

        end = time.time() 
        print 'For freq_like_param = ' + str(freq_like_param) + ', found limit amp_like_param = ' + str(min_results['x'])  +' and it took me ' + str(end - start) + 's to do so.' 

    print 'minimized_values = '
    print minimized_values
        

    

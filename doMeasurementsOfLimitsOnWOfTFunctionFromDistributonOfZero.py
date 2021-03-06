import math
import numpy as np
from measureLimitsOfFunctionGivenDistributionOfZero import measureLimitsOfFunctionGivenDistributionOfZero
from loadSN import loadSN 
from AstronomicalParameterArchive import AstronomicalParameterArchive
from CosmologicalParameterArchive import CosmologicalParameterArchive
from calculateMuForwOftFromz import ResidualMuCalculatorForArbitraryWofT
import matplotlib.pyplot as plt
from StandardPeriodigram import StandardPeriodigram 

if __name__ == '__main__':
    #Number of sigma
    n_sigma = 7.0 
    
    #Must be a two parameter function
    # (ideally one is like a frequency and one is like an amplitude)

    #Start with a simple oscillitory function
    omega0 = 0.05
    bounds = [0.0, 100.0]
    t_or_tau = 'tau'
    if t_or_tau.lower() in ['tau','taus']: 
        wOfTTauZ = lambda ts, taus, zs, A, omega, phi: -1.0 + A * np.sin(np.array(taus)*omega + phi)
    else:
        wOfTTauZ = lambda ts, taus, zs, A, omega, phi: -1.0 + A * np.sin(np.array(ts)*omega + phi)

    #Function to limit should be of three parameters: lookback time, param that will be minimized over, extra param to vary
    #function_to_limit = lambda ts, Amp_like, freq_like: muOfT(ts, Amp_like, freq_like, 0.0)
    
    file_of_zero_limits = '/Users/sasha/Documents/Harvard/physics/stubbs/SNIsotropyProject/randomBestFitParams/params_from_300_frequency_scaling50_10_sampling_of_fitted_alt_gaussian_to_zero.csv'

    results_for_minimization = np.genfromtxt(file_of_zero_limits, delimiter = ',').transpose()
    frequencies = results_for_minimization[0]
    
    H0_in_inv_Mpc = 6.883156370706416 * 10 ** (-5.0)
    #freq_like_params_on_which_to_minimize = [frequencies[5 + i * 100] / H0_in_inv_Mpc for i in range(len(frequencies[0:1000]) / 100)]
    freq_like_params_on_which_to_minimize = [0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 200.0, 800.0]
    #freq_indeces_for_finding_peak = [[max(0, 5 + i * 50 + freq_index_range_for_finding_peak[0]),
    #                                  min(len(frequencies), 5 + i * 50 + freq_index_range_for_finding_peak[1])]
    #                                 for i in range(len(frequencies) / 50)]
    
    init_guesses = [1.0 * freq_like_param / omega0 for freq_like_param in freq_like_params_on_which_to_minimize]

    astro_arch = AstronomicalParameterArchive()
    s_to_yr = astro_arch.getSecondToYear() 
    cosmo_archive = CosmologicalParameterArchive()
    age_of_universe = cosmo_archive.getAgeOfUniverse( units = 'yr')[0]
    
    all_sn = loadSN(1, pull_extinctions = 0)
    all_surveys = np.unique([sn['survey'] for sn in all_sn])
    all_zs = [sn['z'] for sn in all_sn]
    all_ts = [sn['t'] for sn in all_sn]
    all_taus = [sn['tau'] for sn in all_sn]
    all_muErrs = [sn['muErr'] for sn in all_sn]
    all_muResids = [sn['muDiff'] for sn in all_sn]
    all_zerodMuResids = [sn['muDiff'] - sn['muDiffWMean'] for sn in all_sn]
    sorted_zs = sorted(all_zs)
    if t_or_tau.lower() in ['tau','taus']:
        all_lookBacks = [10 ** (-6.0) * tau * s_to_yr for tau in all_taus]
    else:
        all_lookBacks = [10 ** (-6.0) * t * s_to_yr for t in all_ts]
    sorted_lookBacks = [t_or_tau for _,t_or_tau in sorted(zip(all_zs, all_lookBacks))]
    sorted_muErrs = [err for _,err in sorted(zip(all_zs, all_muErrs))]
    sorted_zResids = [res for _,res in sorted(zip(all_zs, all_zerodMuResids))]

    canon_mus = ResidualMuCalculatorForArbitraryWofT(wOfFunction = lambda ts, taus, zs: -1.0 if type(ts) in [int, float] else -1.0 + np.zeros(np.shape(ts)), initial_zs = sorted_zs).getMus()

    #I think that we can define the function to just take in sorted ts, 
    # which I know maps one to one with sorted zs, and then return the calculation with the zs.
    function_to_limit = lambda sorted_lookBacks, Amp_like, freq_like: (np.array(ResidualMuCalculatorForArbitraryWofT(wOfFunction = lambda ts, taus, zs:
                                                                                                                     wOfTTauZ(ts, taus, zs, Amp_like, freq_like, 0.0), initial_zs = sorted_zs).getMus())
                                                                       -  np.array(canon_mus))

    minimized_values = []
                                             
    for i in range(len(freq_like_params_on_which_to_minimize)):
        freq_like_param = freq_like_params_on_which_to_minimize[i]
        init_guess = init_guesses[i]
        function_to_limit_at_specific_freq_val = lambda sorted_lookBacks, Amp_like: function_to_limit (sorted_lookBacks, Amp_like, freq_like_param)
        print 'Working on freq_like_param =' + str(freq_like_param)
        
        #print 'Plotting funct...'
        #plt.scatter(sorted_lookBacks, function_to_limit_at_specific_freq_val(sorted_lookBacks, init_guess))
        #plt.show()
        
        #print 'making periodigram...'
        #full_periodigram = StandardPeriodigram(sorted_lookBacks, function_to_limit_at_specific_freq_val(sorted_lookBacks, init_guess), sorted_muErrs, 
        #                                              frequencies = frequencies, apply_normalization = 0, compute_convolution = 0)
        #print 'Plotting periodigram...'
        #plt.scatter(full_periodigram.frequencies, full_periodigram.normalized_coef_mags)
        #plt.show() 
        
        min_results = measureLimitsOfFunctionGivenDistributionOfZero(sorted_lookBacks, sorted_muErrs, function_to_limit_at_specific_freq_val, init_guess, n_sigma, file_of_zero_limits, 'no saving ability yet', bounds, freq_width_to_calc = 100, freq_bounds_over_which_to_find_peak = [10, len(frequencies) / 10], freq_bounds_over_which_to_check_peak = [5, len(frequencies) / 10])
        minimized_values = minimized_values + [min_results['x']]

        print 'For freq_like_param = ' + str(freq_like_param) + ', found limit amp_like_param = ' + str(min_results['x']) 

    print 'minimized_values = '
    print minimized_values 
        

    

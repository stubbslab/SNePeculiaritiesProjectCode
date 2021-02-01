import numpy as np
import math
from loadSN import loadSN
from CosmologicalParameterArchive import CosmologicalParameterArchive
from AstronomicalParameterArchive import AstronomicalParameterArchive
from calculateFrequencies import calcFrequenciesFromMinSpatialDistance
from computeDistributionOfFourierModesOfFunction import  batchComputeDistributionOfFourierModesGivenFunction
from DirectoryArchive import DirectoryArchive
import time
import os.path

if __name__ == '__main__':
    astro_arch = AstronomicalParameterArchive()
    cosmo_arch = CosmologicalParameterArchive()

    s_to_yr = astro_arch.getSecondToYear()
    age_of_universe = cosmo_arch.getAgeOfUniverse(units = 'year')[0]

    t_or_tau = 'tau' 
    all_sn = loadSN(1)
    all_ts = [sn['t'] for sn in all_sn]
    all_taus = [sn['tau'] for sn in all_sn]
    all_errs = [sn['muErr'] for sn in all_sn]
    all_zResids = [sn['muDiff'] - sn['muDiffWMean'] for sn in all_sn]
    if t_or_tau.lower() in ['tau', 'taus']: 
        all_lookBacks = [10 ** -6.0 * (tau * s_to_yr) for tau in all_taus]
    else:
        all_lookBacks = [10 ** -6.0 * (tau * s_to_yr) for t in all_ts]

    model_type = 'alt_gaussian'
    #model_type = 'exponential'

    freq_density_scaling = 5.0
    freqs_to_calculate = calcFrequenciesFromMinSpatialDistance(freq_density_scaling = freq_density_scaling, t_or_tau = t_or_tau )

    zero_funct = lambda xs: 0.0 if type(xs) in [int, float] else np.zeros(np.shape(xs))
    n_draws = 280
    print 'About to start batch computation.'
    start = time.time()
    batch_length = 20
    frequency_batches = [ freqs_to_calculate[i * batch_length:(i+1) * batch_length] for i in range(len(freqs_to_calculate) / batch_length) ]
    if (len(freqs_to_calculate) / batch_length) * batch_length < len(freqs_to_calculate) - 1: 
        frequency_batches = frequency_batches + [freqs_to_calculate[(len(freqs_to_calculate) / batch_length) * batch_length:]]
    dir_archive = DirectoryArchive()
    dir_to_save = dir_archive.getRandomizedBestFitParamsDir()
    file_name = 'params_from_' + str(n_draws) + '_frequency_scaling' + str(int(freq_density_scaling) * 10) + '_10_sampling_of_fitted_' + t_or_tau + model_type + '_to_zero.csv'
    for i in range(len(frequency_batches)):
        print 'Working on MACRO batch ' + str(i+1) + ' of ' + str(len(frequency_batches)) + ' which contains frequencies: '
        frequency_batch = frequency_batches[i]
        print frequency_batch 
        test_results = batchComputeDistributionOfFourierModesGivenFunction(all_lookBacks, all_zResids, all_errs, zero_funct, frequency_batch, n_redraw = n_draws, show_single_batch_fit = 0, model_CPB_funct = 'zero_alt_err_funct', best_guess_minimization_params = 'alt_normal')
        new_results = [test_results[0], test_results[1], [elem[0] for elem in test_results[2]], [elem[1] for elem in test_results[2]]]

        if os.path.isfile(dir_to_save + file_name):
            old_results = np.genfromtxt(dir_to_save + file_name, delimiter=",").transpose().tolist()
        else: 
            old_results = [[] for elem in new_results]
        results_to_save = [old_results[i] + new_results[i] for i in range(len(new_results))]
        print 'Saving partial array to ' + str(dir_archive.getRandomizedBestFitParamsDir() + file_name)
        np.savetxt(dir_archive.getRandomizedBestFitParamsDir() + file_name, np.array(results_to_save).transpose(), delimiter = ',')

    print 'Done saving all arrays. '
                                            
                                             
    end = time.time()
    print 'Finished with batch computation. Took ' + str(end - start) + 's to do so.' 


    

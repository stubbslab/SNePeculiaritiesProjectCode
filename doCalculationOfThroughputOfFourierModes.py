#We want to fit a spline of a specified number of points to the SN residual data

import numpy as np
from loadSN import loadSN
import math
import random
import time 
import matplotlib.pyplot as plt 
from ApproximateFourierTransform import ApproximateFourierTransform
from calculateFrequencies import calcFrequenciesFromMinSpatialDistance 
from StandardPeriodigram import StandardPeriodigram 
from AstronomicalParameterArchive import AstronomicalParameterArchive 
from CosmologicalParameterArchive import CosmologicalParameterArchive 
from computeDistributionOfFourierModesOfFunction import showDistributionOfFourierModesOfFunction

if __name__ == '__main__':
    simul_real_data = 0
    take_medians_by_freq = 1
    use_real_errs = 0
    data_set = 1
    sin_amplitude = 1.0
    n_redraw_to_determine_mag_distr = 5
    #n_sine_waves = 10
    phases = np.linspace(0.0, 1.0 * math.pi, 4)
    min_standard_fourier_f = 0.000583636085742 * 10.0
    max_standard_fourier_f = 3.0
    in_freqs = np.linspace(min_standard_fourier_f, max_standard_fourier_f, 10)
    windowing_funct = 'rect'
    astro_arch = AstronomicalParameterArchive()
    s_to_yr = astro_arch.getSecondToYear() 
    cosmo_archive = CosmologicalParameterArchive()
    age_of_universe = cosmo_archive.getAgeOfUniverse( units = 'yr')[0] 
    
    all_sn = loadSN(1, pull_extinctions = 0)
    all_surveys = np.unique([sn['survey'] for sn in all_sn])
    all_ts = [sn['t'] for sn in all_sn]
    all_lookBacks = [-1.0 * 10 ** (-6.0) * (t * s_to_yr - age_of_universe) for t in all_ts]
    all_muErrs = [sn['muErr'] for sn in all_sn]
    all_muResids = [sn['muDiff'] for sn in all_sn]
    all_zerodMuResids = [sn['muDiff'] - sn['muDiffWMean'] for sn in all_sn]
    #for survey in all_surveys:
    #    survey_ts = [sn['t'] for sn in all_sn if sn['survey'] == survey]
    #    ordered_ts = ordered_ts + survey_ts 
    #    survey_muErrs = [sn['muErr'] for sn in all_sn if sn['survey'] == survey]
    #    ordered_muErrs = ordered_muErrs + survey_muErrs
    #    survey_muResids = [sn['muDiff'] for sn in all_sn if sn['survey'] == survey]
    #    ordered_muResids = ordered_muResids + survey_muResids
    #    weights = [1.0 / err ** 2.0 for err in survey_muErrs]
    #    wMean = sum([weights[i] * survey_muResids[i] for i in range(len(weights))]) / sum(weights)
    #    ordered_zerodMuResids = ordered_zerodMuResids + [resid - wMean for resid in survey_muResids]

    min_x = min(all_lookBacks)
    max_x = max(all_lookBacks)
    x_interval = (max_x - min_x)
    n_sn = len(all_lookBacks) 
    fourier_freq_density_scaling = 2.0 
    fourier_freqs = calcFrequenciesFromMinSpatialDistance( freq_density_scaling = fourier_freq_density_scaling, min_length_scale = 2.0 * 10.0 ** 6.0 * 20)
    #fourier_freqs = (  [small_val * standard_freq_step for small_val in np.linspace(0.1, 1.0, 10)]
    #                 + [intermediate_val * standard_freq_step for intermediate_val in range(1, n_sn / 2)]
    #                 + [large_val * standard_freq_step for large_val in range(n_sn / 2, 2 * n_sn)]
    #                )

    if use_real_errs:
        errs_to_use = all_muErrs
    else:
        errs_to_use = [1.0 for err in all_muErrs]

    coef_mag_results = []
    Myr_interval = max(all_lookBacks) - min(all_lookBacks)
    freq_range = [2.0 * math.pi / Myr_interval, 2.0 * math.pi / Myr_interval * len(all_lookBacks) / 2.0]

    #phase_inFreq_tuples = zip(phases, in_freqs)

    phase_mesh, in_freq_mesh = np.meshgrid(phases, in_freqs)
    phase_inFreq_tuples = zip(phase_mesh.flatten(), in_freq_mesh.flatten())
    
    fourier_freqs_at_points = []
    fourier_mags_at_points = [] 

    for phase_inFreq_tuple in phase_inFreq_tuples:
        start = time.time() 
        phase = phase_inFreq_tuple[0]
        freq = phase_inFreq_tuple[1] 
        print 'Working on phase ' + str(phase) + ' and frequency = ' + str(freq)
        sin_funct = lambda xs: sin_amplitude * np.sin(freq * np.array(xs) + phase)
        if simul_real_data:  
            fourier_freqs, fourier_mags = showDistributionOfFourierModesOfFunction(all_lookBacks, all_zerodMuResids, err_to_use, sin_funct,
                                                                              n_redraw = n_redraw_to_determine_mag_distr , windowing_funct = windowing_funct, 
                                                                              show = 0, save = 0, return_mags = 1, compute_total_coef_at_freq = 1
                                                                             )
        else:
            sin_vals = sin_funct(all_lookBacks)
            direct_transform = StandardPeriodigram(all_lookBacks, sin_vals, all_muErrs, frequencies = fourier_freqs, windowing_funct = windowing_funct, apply_normalization = 1)
            #direct_transform = ApproximateFourierTransform(all_lookBacks, sin_vals, errs_to_use, frequencies = fourier_freqs, windowing_funct = windowing_funct, apply_normalization = 1)
            fourier_mags = direct_transform.normalized_coef_mags
        fourier_freqs_at_points = fourier_freqs_at_points + [fourier_freqs]
        fourier_mags_at_points = fourier_mags_at_points + [fourier_mags]
        end = time.time()
        print 'Took ' + str(end - start) + 's' 
    fourier_freqs_at_points = np.array(fourier_freqs_at_points)
    fourier_mags_at_points = np.array(fourier_mags_at_points)
    if take_medians_by_freq:
        median_mags = np.median(fourier_mags_at_points, axis = 0)
        zerod_mags_at_points = ((fourier_mags_at_points) - median_mags)
    else:
        median_mags = np.median(fourier_mags_at_points, axis = 1)
        zerod_mags_at_points = ((fourier_mags_at_points).transpose() - median_mags).transpose()
    sum_of_raw_mags_at_points = np.sum(fourier_mags_at_points, axis = 1)
    sum_of_zerod_mags_at_points = np.sum(zerod_mags_at_points, axis = 1)
    max_mag_at_points = np.max(zerod_mags_at_points, axis = 1)
    indeces_of_max_mag_at_points = np.argmax(zerod_mags_at_points, axis = 1)
    peak_freq_at_points = np.array([ fourier_freqs_at_points[i][indeces_of_max_mag_at_points[i]] for i in range(len(indeces_of_max_mag_at_points)) ])

    print 'peak_freq_at_points = ' + str(peak_freq_at_points)
    print 'max_mag_at_points = ' + str(max_mag_at_points)
    print 'sum_of_raw_mags_at_points = ' + str(sum_of_raw_mags_at_points) 
    print 'sum_of_zerod_mags_at_points = ' + str(sum_of_zerod_mags_at_points) 
    print 'phase_inFreq_tuples = ' + str(phase_inFreq_tuples)

    peak_freq_mesh = np.reshape(peak_freq_at_points, np.shape(phase_mesh))
    max_mag_mesh = np.reshape(max_mag_at_points, np.shape(phase_mesh))
    sum_of_raw_mags_mesh = np.reshape(sum_of_raw_mags_at_points, np.shape(phase_mesh)) 
    sum_of_zerod_mags_mesh = np.reshape(sum_of_zerod_mags_at_points, np.shape(phase_mesh)) 

    f, axarr = plt.subplots(2, 4)
    #axarr[0, 0].contour(phase_mesh, in_freq_mesh, peak_freq_mesh)
    #axarr[0, 0].set_title('frequency of maximum mode contours')
    #axarr[0, 1].contour(phase_mesh, in_freq_mesh, max_mag_mesh)
    #axarr[0, 1].set_title('maximum magnitude contours')
    #axarr[0, 2].contour(phase_mesh, in_freq_mesh, sum_of_mags_mesh)
    #axarr[0, 2].set_title('sum of all magnitude contours')

    for i in range(3):
        axarr[0, i].set_xlabel('phase')
        axarr[1, i].set_xlabel('input frequency (1 / Myr) ')
    for i in range(0, 2): 
        axarr[i, 0].set_ylabel('frequency of maximum mode (1 / Myr)')
        axarr[i, 1].set_ylabel('magnitude of maximum mode')
        axarr[i, 2].set_ylabel('sum of magnitude of all modes')
        axarr[i, 3].set_ylabel('median-subtracted sum of magnitudes')  
    axarr[0, 0].scatter(phases, np.mean(peak_freq_mesh, axis = 0))
    axarr[0, 1].scatter(phases, np.mean(max_mag_mesh, axis = 0))
    axarr[0, 2].scatter(phases, np.mean(sum_of_raw_mags_mesh, axis = 0))
    axarr[0, 3].scatter(phases, np.mean(sum_of_zerod_mags_mesh, axis = 0))
    axarr[1, 0].scatter(in_freqs, np.mean(peak_freq_mesh, axis = 1))
    axarr[1, 1].scatter(in_freqs, np.mean(max_mag_mesh, axis = 1))
    axarr[1, 2].scatter(in_freqs, np.mean(sum_of_raw_mags_mesh, axis = 1))
    axarr[1, 3].scatter(in_freqs, np.mean(sum_of_zerod_mags_mesh, axis = 1))
    plt.suptitle('Throughput Function') 
    plt.show() 
    
    
    #        computed_mags = [pair[1] for pair in fourier_freq_mag_pairs]
    #        computed_freqs = [pair[0] for pair in fourier_freq_mag_pairs]
    #    sum_of_mags_pre_med_removal = sum(computed_mags)
    #    peak_mag_pre_med_removal = np.max(computed_mags)
    #    peak_freq_pre_med_removal = computed_freqs[np.argmax(computed_mags)]
    #    coef_mag_results = coef_mag_results + [[(freq, phase), peak_freq_pre_med_removal, peak_mag_pre_med_removal, sum_of_mags_pre_med_removal]]

    #median_mags = [ np.median([freq_mag_pairs[i][1] for freq_mag_pairs in fourier_freq_mag_pairs_at_points]) for i in range(len(fourier_freq_mag_pairs_at_points[0]))]
    #reduced_fourier_freq_mag_pairs_at_points = []
    #for i in range(n_sine_waves):
    #    fourier_freq_mag_pairs = fourier_freq_mag_pairs_at_points[i]
    #    reduced_fourier_mags = [fourier_freq_mag_pairs[j][1] - median_mags[j] for j in range(len(fourier_freq_mag_pairs))]
    #    sum_of_mags_post_med_removal = sum(reduced_fourier_mags)
    #    peak_mag_post_med_removal = np.max(reduced_fourier_mags)
    #    peak_freq_post_med_removal = computed_freqs[np.argmax(reduced_fourier_mags)]
    #    coef_mag_results[i] = coef_mag_results[i] + [peak_freq_post_med_removal, peak_mag_post_med_removal, sum_of_mags_post_med_removal]
    #print 'coef_mag_results = ' + str(coef_mag_results) 

      

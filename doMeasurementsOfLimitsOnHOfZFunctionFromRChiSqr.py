import math
import numpy as np
from measureLimitsOfFunctionGivenDistributionOfZero import ComparatorBetweenFunctionAndDistributionOfZero
from cantrips import safeSortOneListByAnother
from measureLimitsOfFunctionGivenDistributionOfZero import measureLimitsOfFunctionGivenDistributionOfZero
from loadSN import loadSN
import scipy.optimize as optimize
import scipy.integrate as integrate 
import scipy.stats as stats
from AstronomicalParameterArchive import AstronomicalParameterArchive
from CosmologicalParameterArchive import CosmologicalParameterArchive
from calculateMuForwOftFromz import ResidualMuCalculatorForArbitraryWofT
from calculateMuForHofZFromz import computeMusForHofZFromz 
import matplotlib.pyplot as plt
from StandardPeriodigram import StandardPeriodigram
from calculateMuForHOfTFromz import ResidualMuCalculatorForArbitraryHofT
import time

if __name__ == '__main__':
    #Number of sigma
    print 'Starting measurement...' 
    start = time.time() 

    t_or_tau = 'tau' 
    H0_in_inv_Myr = 6.883156370706416 * 10 ** (-5.0)
    km_s_Mpc_in_Myr = 1.022 * 10 ** (-6.0) # 1km/s/Mpc in Mega years 
    #freq_like_params_on_which_to_minimize = [frequencies[5 + i * 100] / H0_in_inv_Mpc for i in range(len(frequencies[0:1000]) / 100)]
    freq_like_params_on_which_to_minimize = [0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 200.0, 800.0]
    #freq_indeces_for_finding_peak = [[max(0, 5 + i * 50 + freq_index_range_for_finding_peak[0]),
    #                                  min(len(frequencies), 5 + i * 50 + freq_index_range_for_finding_peak[1])]
    #                                 for i in range(len(frequencies) / 50)]
    

    astro_arch = AstronomicalParameterArchive()
    cosmo_arch = CosmologicalParameterArchive()
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
    all_surveys = [sn['survey'] for sn in all_sn]
    unique_surveys = np.unique(all_surveys)
    #sorted_zs = sorted(all_zs)
    if t_or_tau.lower() in ['tau','taus']:
        all_lookBacks = [10 ** (-6.0) * tau * s_to_yr for tau in all_taus]
    else:
        all_lookBacks = [10 ** (-6.0) * t * s_to_yr for t in all_ts]
    sn_weights = [1.0 / err ** 2.0 for err in all_muErrs]
    full_mean = sum([sn_weights[i] * all_muResids[i] for i in range(len(sn_weights))]) / (sum(sn_weights))
    all_flatZerodMuResids = [resid - full_mean for resid in all_muResids]
    #sorted_lookBacks = [t_or_tau for _,t_or_tau in sorted(zip(all_zs, all_lookBacks))]
    #sorted_muErrs = [err for _,err in sorted(zip(all_zs, all_muErrs))]
    #sorted_zResids = [res for _,res in sorted(zip(all_zs, all_zerodMuResids))]
    subtract_individual_means = 1
    sorted_zs, sorted_lookBacks, sorted_muErrs, sorted_resids, sorted_resids_noBySurveyZeroing, sorted_surveys = safeSortOneListByAnother(all_zs, [all_zs, all_lookBacks, all_muErrs, all_zerodMuResids, all_flatZerodMuResids, all_surveys])     

    #canon_mus = ResidualMuCalculatorForArbitraryWofT(wOfFunction = lambda ts, taus, zs: -1.0 if type(ts) in [int, float] else -1.0 + np.zeros(np.shape(ts)), initial_zs = sorted_zs).getMus()
    canon_mus = ResidualMuCalculatorForArbitraryHofT(HScalingFunction = lambda ts, taus, zs: 1.0 if type(ts) in [float, int] else [1.0 for z in ts], sorted_zs)[3] 

    #I think that we can define the function to just take in sorted ts, 
    # which I know maps one to one with sorted zs, and then return the calculation with the zs.

    phase_pi_coef = 0.0 / 4.0
    phase_shift_int = 0.0
    phase_shift_order = 10.0
    phase = math.pi * phase_pi_coef + float(phase_shift_int) / phase_shift_order
    phase_str = r'$\pi$' + r'$\mathrm{/2}$'
    phase_str = r'$\mathrm{1/10}$'
    phase_str = r'$\mathrm{0}$'
    
    if t_or_tau.lower() in ['tau','taus']:
        HofZScaling = lambda ts, taus, zs, A, freq_like, phase: 1.0 + A * np.sin(freq_like * taus + phase)
    else:
        HofZScaling = lambda ts, taus, zs, A, freq_like, phase: 1.0 + A * np.sin(freq_like * ts + phase)
        
    function_to_limit = lambda sorted_lookBacks, A, freq_like: np.array(ResidualMuCalculatorForArbitraryHofT(HScalingFunction = lambda ts, taus, zs:
                                                                                                             HofZScaling(ts, taus, zs, A, freq_like, phase),
                                                                                                             initial_zs = sorted_zs))[3] 
    
    minimized_values = []
    #Calculate the residual over a double for loop, iterating through both the allowed omegas and allowed As
    f_min_redshift = 1.0 / 2.0
    print 'f_min_redshift = ' + str(f_min_redshift)
    fs_to_calc_in_redshift = np.linspace(0.0 + f_min_redshift, 21 * f_min_redshift, 4) # 161 81 101 #251 
    omegas_to_calc_in_redshift = fs_to_calc_in_redshift * 2.0 * math.pi
                                                                          
    #omegas_to_calc_in_H0 = omegas_to_calc_in_Myr / H0_in_inv_Myr
    print 'omegas_to_calc_in_redshift = ' + str(omegas_to_calc_in_redshift)
    freq_half_width_to_find_peak = 100
    As_to_calc = np.linspace(0.0, 0.10, 3) # 201 #101
    A_mesh, omega_mesh = np.meshgrid(omegas_to_calc_in_redshift, As_to_calc)
    #calc_chi_sqrs_yesZBS = np.zeros(np.shape(A_mesh))
    calc_chi_sqrs_noZBS = np.zeros(np.shape(A_mesh))
    #peak_CPB_vals_yesZBS = np.zeros(np.shape(A_mesh))
    peak_CPB_vals_noZBS = np.zeros(np.shape (A_mesh))
    limits_file = '/Users/sasha/Documents/Harvard/physics/stubbs/SNIsotropyProject/randomBestFitParams/params_from_1000_frequency_scaling50_10_max_oscillations_1000_sampling_of_fitted_tauexponential_to_zero_updated_SN_subtract_wMean.csv'
    n_fit_params = 3
    dof_zeroBySurvey = len(all_sn) - n_fit_params - len(unique_surveys) #must account for individual shifts from 0 as fit parameters
    #dof_zeroBySurvey = len(all_sn) - n_fit_params - 1
    dof_noZeroBySurvey = len(all_sn) - n_fit_params - 1 #must account for overall shifts from 0 as fit parameter 
    
    periodigram_comparator = ComparatorBetweenFunctionAndDistributionOfZero(limits_file, '', fitting_funct_type = 'alt_normal')
    canon_H0 = cosmo_arch.getH0()
    #min_H0_scaling = (canon_H0[0] - 5.0 * canon_H0[1]) / canon_H0[0]
    #max_H0_scaling = (canon_H0[0] + 5.0 * canon_H0[1]) / canon_H0[0] 
    for i in range(len(As_to_calc)):
        macro_element_start = time.time() 
        A = As_to_calc[i]
        print 'Working on i ' + str(i) + ' of ' + str(len(As_to_calc)) + ' at which A = ' + str(A)
        for j in range(len(omegas_to_calc_in_redshift)):
            element_start = time.time()
            omega_in_redshift = omegas_to_calc_in_redshift[j]
            print 'Working on i = ' + str(i) + ' at which A = ' + str(A) + ' and j = ' + str(j) + ' of ' + str(len(fs_to_calc_in_redshift)) + ' at which fs_in_redshift = ' + str(fs_to_calc_in_redshift) 
            new_H0_scaling = 1.0
            new_H0 = H0_in_inv_Myr * new_H0_scaling
            
            raw_theoretical_mu_resids = function_to_limit(sorted_zs, A, omega_in_redshift, phase) - np.array(canon_mus) 
            mean_theoretical_mu_resids = np.mean(raw_theoretical_mu_resids)
            theoretical_mu_resids_noZBS = [resid - mean_theoretical_mu_resids for resid in raw_theoretical_mu_resids]
            #plt.scatter(sorted_lookBacks, theoretical_mu_resids_noZBS)
            #plt.show()

            mean_theoretical_mu_resids_by_survey = {}
            minimizing_offset_by_survey = {} 
            for survey in unique_surveys:
                true_muDiffs_by_survey = [sorted_resids_noBySurveyZeroing[k] for k in range(len(sorted_resids_noBySurveyZeroing)) if sorted_surveys[k] == survey]
                theoretical_muDiffs_by_survey = [raw_theoretical_mu_resids[k] for k in range(len(raw_theoretical_mu_resids)) if sorted_surveys[k] == survey]
                muErrs_by_survey = [sorted_muErrs[k] for k in range(len(raw_theoretical_mu_resids)) if sorted_surveys[k] == survey]
                weights_by_survey = [ 1.0 / err ** 2.0 for err in muErrs_by_survey ]
                weighted_mean_for_survey = sum( [ theoretical_muDiffs_by_survey[k] * weights_by_survey[k] ** 2.0 for k in range(len(theoretical_muDiffs_by_survey)) ] ) / sum(weights_by_survey)
                minimizing_offset_by_survey[survey] = (-1.0 / sum([weight ** 2.0 for weight in weights_by_survey])
                                                       * sum ([(true_muDiffs_by_survey[k] - theoretical_muDiffs_by_survey[k]) * weights_by_survey[k] for k in range(len(true_muDiffs_by_survey))]))
                mean_theoretical_mu_resids_by_survey[survey] = weighted_mean_for_survey
            #theoretical_mu_resids_yesZBS = [ raw_theoretical_mu_resids[k] - mean_theoretical_mu_resids_by_survey[sorted_surveys[k]] for k in range(len(raw_theoretical_mu_resids)) ]
            #theoretical_mu_resids_yesZBS = [ raw_theoretical_mu_resids[k] - minimizing_offset_by_survey[sorted_surveys[k]] for k in range(len(raw_theoretical_mu_resids)) ]
            
            #new_chi_sqr_yesZBS = sum([((theoretical_mu_resids_yesZBS[k] - sorted_resids[k]) / (sorted_muErrs[k])) ** 2.0 for k in range(len(sorted_zs))]) / dof_zeroBySurvey 
            new_chi_sqr_noZBS = sum([((theoretical_mu_resids_noZBS[k] - sorted_resids_noBySurveyZeroing[k]) / (sorted_muErrs[k])) ** 2.0 for k in range(len(sorted_zs))]) / dof_noZeroBySurvey
            approx_omega_in_Planck_Myr = (sorted_zs[-1] - sorted_zs[0]) / (sorted_lookBacks[-1] - sorted_lookBacks[0]) * omega_in_redshift
            print 'approx_omega_in_Planck_Myr = ' + str(approx_omega_in_Planck_Myr )
            init_peak_guess = np.argmin(np.abs(np.array(periodigram_comparator.frequencies) - approx_omega_in_Planck_Myr))
            freq_bounds_over_which_to_find_peak = [max(0, init_peak_guess - freq_half_width_to_find_peak),
                                                   min(len(periodigram_comparator.frequencies), init_peak_guess + freq_half_width_to_find_peak)]
            #print 'init_peak_guess = ' + str(init_peak_guess)
            #print 'freq_bounds_over_which_to_find_peak = ' + str(freq_bounds_over_which_to_find_peak)
            start_freqs_of_interest = [0, min(100, freq_bounds_over_which_to_find_peak[0])] #Freqs to check since largest Fourier mode can sometimes be in early distribution
            #print 'start_freqs_of_interest  = ' + str(start_freqs_of_interest )

            #normalized_coefs_around_guess_yesZBS, CPB_vals_of_zero_distributions_around_guess_yesZBS = periodigram_comparator.computeLikelihoodsOfSignalFromZeros(sorted_lookBacks,
            #                                                                                                                                                      theoretical_mu_resids_yesZBS,
            #                                                                                                                                                      sorted_muErrs,
            #                                                                                                                                                      freq_bounds_over_which_to_find_peak = freq_bounds_over_which_to_find_peak)
            normalized_coefs_around_guess_noZBS, CPB_vals_of_zero_distributions_around_guess_noZBS = periodigram_comparator.computeLikelihoodsOfSignalFromZeros(sorted_lookBacks,
                                                                                                                                                                theoretical_mu_resids_noZBS,
                                                                                                                                                                sorted_muErrs,
                                                                                                                                                                freq_bounds_over_which_to_find_peak = freq_bounds_over_which_to_find_peak)
            if start_freqs_of_interest[1] > 0: 
                #normalized_coefs_at_start_freqs_yesZBS, CPB_vals_of_zero_distributions_at_start_freqs_yesZBS = periodigram_comparator.computeLikelihoodsOfSignalFromZeros(sorted_lookBacks,
                #                                                                                                                                                          theoretical_mu_resids_yesZBS,
                #                                                                                                                                                          sorted_muErrs,
                #                                                                                                                                                          freq_bounds_over_which_to_find_peak = start_freqs_of_interest)
                normalized_coefs_at_start_freqs_noZBS, CPB_vals_of_zero_distributions_at_start_freqs_noZBS = periodigram_comparator.computeLikelihoodsOfSignalFromZeros(sorted_lookBacks,
                                                                                                                                                                        theoretical_mu_resids_noZBS,
                                                                                                                                                                        sorted_muErrs,
                                                                                                                                                                        freq_bounds_over_which_to_find_peak = start_freqs_of_interest)
            else:
                #normalized_coefs_at_start_freqs_yesZBS = []
                #CPB_vals_of_zero_distributions_at_start_freqs_yesZBS = []
                normalized_coefs_at_start_freqs_noZBS = []
                CPB_vals_of_zero_distributions_at_start_freqs_noZBS = []
                

            #normalized_coefs_yesZBS = normalized_coefs_at_start_freqs_yesZBS + normalized_coefs_around_guess_yesZBS
            #CPB_vals_of_zero_distributions_yesZBS  = CPB_vals_of_zero_distributions_at_start_freqs_yesZBS + CPB_vals_of_zero_distributions_around_guess_yesZBS 
            normalized_coefs_noZBS = normalized_coefs_at_start_freqs_noZBS + normalized_coefs_around_guess_noZBS
            CPB_vals_of_zero_distributions_noZBS  = CPB_vals_of_zero_distributions_at_start_freqs_noZBS + CPB_vals_of_zero_distributions_around_guess_noZBS 
            
            

            peak_CPB_val_noZBS = 10.0 ** max(np.log10(CPB_vals_of_zero_distributions_noZBS))
            peak_CPB_val_noZBS_index = np.argmax(np.log10(CPB_vals_of_zero_distributions_noZBS))
            #peak_CPB_val_yesZBS = 10.0 ** max(np.log10(CPB_vals_of_zero_distributions_yesZBS))
            #peak_CPB_val_yesZBS_index = np.argmax(np.log10(CPB_vals_of_zero_distributions_yesZBS))
            #calc_chi_sqrs_yesZBS[i][j] = new_chi_sqr_yesZBS
            calc_chi_sqrs_noZBS[i][j] = new_chi_sqr_noZBS
            peak_CPB_vals_noZBS[i][j] = peak_CPB_val_noZBS
            #peak_CPB_vals_yesZBS[i][j] = peak_CPB_val_yesZBS
            element_end = time.time()
            if j % 50 == 0:
                print 'j ' + str(j) + ' of ' + str(len(omegas_to_calc_in_redshift)) + ' took ' + str(element_end - element_start) + 's'
                #print '    new_H0_scaling_first_guess * np.array((new_H0_min_bound_coefs[0], 1.0, new_H0_min_bound_coefs[1])) = ' + str(new_H0_scaling_first_guess  * np.array((new_H0_min_bound_coefs[0], 1.0, new_H0_min_bound_coefs[1]))) 
                #print '    new_H0_scaling = ' + str(new_H0_scaling)
                #print '    This quantity should be (nearly) 0.0: ' + str(function_to_constrain_at_H0(A, omega_in_Myr / (new_H0_scaling * H0_in_inv_Myr), H0_in_inv_Myr * new_H0_scaling))
                #print '    omega_in_Myr / (new_H0_scaling * H0_in_inv_Myr) = ' + str(omega_in_Myr / (new_H0_scaling * H0_in_inv_Myr))
                #print '    omegas_to_calc_in_H0[j]  = ' + str(omegas_to_calc_in_H0[j])
                #print '    peak_CPB_val = ' + str(peak_CPB_val)
                #print '    peak_CPB_val_index = ' + str(peak_CPB_val_index)
 
        macro_element_end = time.time() 
        print 'i ' + str(i) + ' of ' + str(len(As_to_calc)) + ' took  ' + str(macro_element_end - macro_element_start) + 's'
    #plt.scatter(sorted_zs, theoretical_mu_resids)
    #plt.scatter(sorted_zs, sorted_resids)
    #plt.errorbar(sorted_zs, sorted_resids, yerr = sorted_muErrs, fmt = None)
    #plt.show()
    res_limites_from_fourier = []

    file_name = 'ForHSin_phi_pix' + str(np.around(phase_pi_coef, 3)) + '+' + str(int(phase_shift_int)) + '_' + str(phase_shift_order) +  'withFourier_correctedDOF' + '_final11'
    #np.save('/Users/sasha/Desktop/RChiSqrLimits/' + 'rChiSqr_yesZBS' + file_name, calc_chi_sqrs_yesZBS)
    np.save('/Users/sasha/Desktop/RChiSqrLimits/' + 'rChiSqr_noZBS' + file_name, calc_chi_sqrs_noZBS) 
    #np.save('/Users/sasha/Desktop/RChiSqrLimits/' + 'peakCPBVals_yesZBS' + file_name, peak_CPB_vals_yesZBS)
    np.save('/Users/sasha/Desktop/RChiSqrLimits/' + 'peakCPBVals_noZBS' + file_name, peak_CPB_vals_noZBS)

    ref_r_chi_sqr = sum([((0.0 - sorted_resids_noBySurveyZeroing[k]) / (sorted_muErrs[k])) ** 2.0 for k in range(len(sorted_zs))]) / dof_noZeroBySurvey
    print 'ref_r_chi_sqr = ' + str(ref_r_chi_sqr) 
    ref_r_chi_sqr_prob = 1.0 - stats.chi2.cdf(dof_noZeroBySurvey * ref_r_chi_sqr, dof_noZeroBySurvey)
    print 'ref_r_chi_sqr_prob = ' + str(ref_r_chi_sqr_prob)

    N_peaks = 500.0 
    ref_fourier_max_CPB_value = 0.9993306871
    ref_fourier_CPB_prob = 1.0 - ref_fourier_max_CPB_value ** (N_peaks)
    #ref_fourier_CPB_prob = 1.0 
    print 'ref_fourier_CPB_prob = ' + str(ref_fourier_CPB_prob) 

    target_probs = [0.0, 0.5, 0.9, 0.99, 0.999, 0.9999]
    target_relative_r_chi_sqr_probs = [1.0 / ref_r_chi_sqr_prob] + [1.0 - prob for prob in target_probs]
    target_relative_fourier_probs =   [1.0 / ref_fourier_CPB_prob] + [1.0 - prob for prob in target_probs]
    print 'target_relative_fourier_probs = ' + str(target_relative_fourier_probs)

    threshold_r_chi_sqr_vals_noZeroBySurvey = [optimize.minimize_scalar(lambda rchisqr: abs((1.0 - stats.chi2.cdf(dof_noZeroBySurvey * rchisqr, dof_noZeroBySurvey)) / ref_r_chi_sqr_prob - prob),
                                                         bounds = [0.1, 2.0], method = 'bounded' )['x'] for prob in target_relative_r_chi_sqr_probs]

    print 'threshold_r_chi_sqr_vals_noZeroBySurvey  = ' + str(threshold_r_chi_sqr_vals_noZeroBySurvey )
    
    
    #levels = np.array([target_prob ** (1.0 / N_peaks) for target_prob in target_probs])
    #contours = plt.contourf(A_mesh, omega_mesh, (1.0 - peak_CPB_vals_yesZBS ** N_peaks) / ref_fourier_CPB_prob,
    #                        levels =list(reversed( target_relative_fourier_probs)), colors = list(reversed(['g','b','m','r','orange'])), alpha = 0.5)
    #plt.xticks(np.linspace(0.0, 0.005, 11) * 2.0 * math.pi / (km_s_Mpc_in_Myr * 100.0), np.linspace(0.0, 0.005, 11))
    #plt.xlim(min(omegas_to_calc_in_H0), max(omegas_to_calc_in_H0)) 
    #label_size = 20.0
    #plt.xlabel(r'$f$' + r'$_w / h$' + ' (1/Myr)', fontsize = label_size, labelpad = label_size / 2.0)
    #plt.ylabel('$A$' + r'$_w$', fontsize = label_size)
    #plt.text((np.max(omegas_to_calc_in_H0) + np.mean(omegas_to_calc_in_H0))/2.5, np.mean(As_to_calc), r'$\phi_w$' + ' = ' + phase_str, fontsize = label_size )
    #plt.title('',  fontsize = label_size)
    #cbar = plt.colorbar(contours)
    #cbar.ax.set_yticklabels([r'$P_{\mathrm{rej}}$ = ' + str(np.around(1.0 - prob, 5)) for prob in target_probs])
    #plt.ylim([min(As_to_calc), max(As_to_calc)])
    #plt.savefig('/Users/sasha/Desktop/RChiSqrLimits/' + 'peakCPBVals_yesZBS' + file_name + '.png')
    #plt.close('all')

    print '(1.0 - peak_CPB_vals_noZBS ** N_peaks) / ref_fourier_CPB_prob = '
    print (1.0 - peak_CPB_vals_noZBS ** N_peaks) / ref_fourier_CPB_prob

    print 'threshold_r_chi_sqr_vals_noZeroBySurvey  = ' +str(threshold_r_chi_sqr_vals_noZeroBySurvey )
    chi_sqr_noZBS_contours = plt.contourf(A_mesh, omega_mesh, calc_chi_sqrs_noZBS, levels = threshold_r_chi_sqr_vals_noZeroBySurvey , colors = ['yellow','g','b','m','r','orange'], alpha = 0.5)
    plt.xticks(np.linspace(0.0, max(omegas_to_calc_in_redshift), 11))
    plt.xlim(min(omegas_to_calc_in_redshift), max(omegas_to_calc_in_redshift)) 
    label_size = 20.0
    plt.xlabel(r'$f$' + r'$_H$' + ' ( frequencies in redshift space)', fontsize = label_size, labelpad = label_size / 2.0)
    plt.ylabel('$A$' + r'$_H$', fontsize = label_size)
    plt.text((np.max(omegas_to_calc_in_redshift) + np.mean(omegas_to_calc_in_redshift))/2.5, np.mean(As_to_calc), r'$\phi_H$' + ' = ' + phase_str, fontsize = label_size )
    plt.title('',  fontsize = label_size)
    cbar = plt.colorbar(chi_sqr_noZBS_contours) 
    cbar.ax.set_yticklabels([''] + ['$R_{\mathrm{rej}}$ = ' + str(1.0 - np.around(prob, 5)) for prob in target_probs])
    plt.ylim([min(As_to_calc), max(As_to_calc)])
    plt.tight_layout()
    plt.savefig('/Users/sasha/Desktop/RChiSqrLimits/' + 'chSqrVals_noZBS' + file_name + '.png')
    plt.close('all')

    #
    #contours = plt.contourf(A_mesh, omega_mesh, (1.0 - peak_CPB_vals_noZBS ** N_peaks) / ref_fourier_CPB_prob,
    #                        levels = list(reversed(target_relative_fourier_probs)), colors = list(reversed(['yellow', 'g','b','m','r','orange'])), alpha = 0.5)
    print 'target_relative_fourier_probs = ' + str(target_relative_fourier_probs) 
    #Note: in these contours, I have to multiply everything by -1 so that levels are increasing.  plt.contourf seems to only work if contours are increasing.  
    contours = plt.contourf(A_mesh, omega_mesh, -1.0 * (1.0 - peak_CPB_vals_noZBS ** N_peaks) / ref_fourier_CPB_prob,
                            levels = [-1.0 * prob for prob in target_relative_fourier_probs], colors = ['yellow', 'g','b','m','r','orange'], alpha = 0.5)  
    plt.xticks(np.linspace(0.0, max(omegas_to_calc_in_redshift), 11) * 2.0 * math.pi)
    plt.xlim(min(omegas_to_calc_in_redshift), max(omegas_to_calc_in_redshift))  
    label_size = 20.0
    plt.xlabel(r'$f$' + r'$_H$' + ' ( frequencies in redshift space)', fontsize = label_size, labelpad = label_size / 2.0)
    plt.ylabel('$A$' + r'$_H$', fontsize = label_size)
    plt.text((np.max(omegas_to_calc_in_redshift) + np.mean(omegas_to_calc_in_redshift))/2.5, np.mean(As_to_calc), r'$\phi_H$' + ' = ' + phase_str, fontsize = label_size )
    plt.title('',  fontsize = label_size)
    cbar = plt.colorbar(contours)
    cbar.ax.set_yticklabels([''] + [r'$R_{\mathrm{rej}}$ = ' + str(np.around(1.0 - prob, 5)) for prob in target_probs])
    plt.ylim([min(As_to_calc), max(As_to_calc)])
    plt.tight_layout()
    plt.savefig('/Users/sasha/Desktop/RChiSqrLimits/' + 'peakCPBVals_noZBS' + file_name + '.png')
    plt.close('all')
    
    
    #levels = np.linspace(min_r_chi_sqr, threshold_r_chi_sqr, 100)
    #chi_sqr_yesZBS_contours = plt.contourf(A_mesh, omega_mesh, calc_chi_sqrs_yesZBS, levels = threshold_r_chi_sqr_vals_noZeroBySurvey , colors = ['g','b','m','r','orange'], alpha = 0.5)
    #plt.xticks(np.linspace(0.0, 0.005, 11) * 2.0 * math.pi / (km_s_Mpc_in_Myr * 100.0), np.linspace(0.0, 0.005, 11))
    #plt.xlim(min(omegas_to_calc_in_H0), max(omegas_to_calc_in_H0)) 
    #label_size = 20.0
    #plt.xlabel(r'$f$' + r'$_w / h$' + ' (1/Myr)', fontsize = label_size, labelpad = label_size / 2.0)
    #plt.ylabel('$A$' + r'$_w$', fontsize = label_size)
    #plt.text((np.max(omegas_to_calc_in_H0) + np.mean(omegas_to_calc_in_H0))/2.5, np.mean(As_to_calc), r'$\phi_w$' + ' = ' + phase_str, fontsize = label_size )
    #plt.title('',  fontsize = label_size)
    #cbar = plt.colorbar(chi_sqr_yesZBS_contours)
    #cbar.ax.set_yticklabels([''] + ['$P_{\mathrm{rej}}$ = ' + str(1.0 - np.around(prob, 5)) for prob in target_probs[1:]])
    #plt.tight_layout()
    #plt.savefig('/Users/sasha/Desktop/RChiSqrLimits/' + 'chiSqrVals_yesZBS' + file_name + '.png')
    #plt.close('all')

    end = time.time()
    print 'Took ' + str(end - start)  + 's to do computation.'  
    #plt.show()
    


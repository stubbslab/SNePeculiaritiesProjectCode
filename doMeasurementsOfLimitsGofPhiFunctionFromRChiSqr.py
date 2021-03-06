
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
import matplotlib.pyplot as plt
from StandardPeriodigram import StandardPeriodigram
from calculateMuForXOfTauFromz import ResidualMuCalculatorForArbitraryXofT
import time
from scipy.interpolate import interp1d 
from calculateMuForFullOscillatingGravity_V6 import ResidualMuCalculatorForOscillatingG 

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

    nu = 200.0 
    nu_str = r'$\mathrm{200}$'
    #phase_str = r'$\mathrm{1/10}$'
    #phase_str = r'$\mathrm{0}$'
    alpha = 0.5
    alpha_str = '5_10'
    A = 7.0 
    filled_zs = np.linspace(10.0 ** -10.0, max(sorted_zs), 10001)

    canon_mus = interp1d(filled_zs, np.array(ResidualMuCalculatorForOscillatingG (A, nu, alpha, 0.0, 0.0, filled_zs).getMus()))(sorted_zs)

    #canon_mus = ResidualMuCalculatorForArbitraryWofT(wOfFunction = lambda ts, taus, zs: -1.0 if type(ts) in [int, float] else -1.0 + np.zeros(np.shape(ts)), initial_zs = sorted_zs).getMus()
    canon_mus_filled = ResidualMuCalculatorForArbitraryXofT(XOfFunction = lambda ts, taus, zs: 1.0 if type(ts) in [float, int, np.float64] else [1.0 for t in ts], initial_zs = np.array([0.0] + filled_zs.tolist())).getMus()
    canon_interp = interp1d(filled_zs, canon_mus_filled) 
    canon_mus2 = [canon_interp(z) for z in sorted_zs] 

    #I think that we can define the function to just take in sorted ts, 
    # which I know maps one to one with sorted zs, and then return the calculation with the zs.

    function_to_limit = lambda sorted_lookBacks, resid_mu_calculator: ( interp1d(filled_zs, np.array(resid_mu_calculator.getMus()))(sorted_zs) 
                                                                        - np.array(canon_mus) )
    #function_to_limit = lambda sorted_lookBacks, C, psi: ( interp1d(filled_zs, np.array(ResidualMuCalculatorForOscillatingG (A, nu, alpha, psi, C, filled_zs).getMus()))(sorted_zs) 
    #                                                       - np.array(canon_mus) )
      
    
    minimized_values = []
    #Calculate the residual over a double for loop, iterating through both the allowed omegas and allowed As
    C_min = 0.001
    C_max = 0.400 
    #f_min_H0 = 1.0 / integrate.quad(lambda z: 1.0 / np.sqrt((1.0 + z) ** 3.0 * cosmo_archive.getOmegaM()[0]+ (1.0 + z) ** 4.0 * cosmo_archive.getOmegaR()[0]+ cosmo_archive.getOmegaLambda()[0]), min(sorted_zs), max(sorted_zs))[0] #independent of h, H0, etc.  This quantity is f_min / H0
    #print 'f_min_H0 = ' + str(f_min_H0)
    psis_to_calc = np.linspace(0.0 * np.pi, 1.0 * np.pi, 41) # 41
    #omegas_to_calc_in_H0 = fs_to_calc_in_H0 * 2.0 * math.pi
    print 'psis_to_calc = ' + str(psis_to_calc)

    Cs_to_calc = np.linspace(C_min, C_max, 51) #51 np.linspace(-0.45, 0.45, 91) # 91
    C_mesh, psi_mesh = np.meshgrid(psis_to_calc, Cs_to_calc)
    #calc_chi_sqrs_yesZBS = np.zeros(np.shape(A_mesh))
    calc_chi_sqrs_noZBS = np.zeros(np.shape(C_mesh))
    #peak_CPB_vals_yesZBS = np.zeros(np.shape(A_mesh))
    peak_CPB_vals_noZBS = np.zeros(np.shape (C_mesh))
    peak_CPB_noZBS_index = np.zeros(np.shape (C_mesh))
    dG_over_Gs = np.zeros(np.shape(C_mesh))
    ddG_over_Gs = np.zeros(np.shape(C_mesh)) 

    limits_file = '/Users/sasha/Documents/Harvard/physics/stubbs/SNIsotropyProject/randomBestFitParams/params_from_1000_frequency_scaling50_10_max_oscillations_1000_sampling_of_fitted_tauexponential_to_zero_updated_SN_subtract_wMean.csv'
    n_fit_params = 5 
    dof_zeroBySurvey = len(all_sn) - n_fit_params - len(unique_surveys) #must account for individual shifts from 0 as fit parameters
    #dof_zeroBySurvey = len(all_sn) - n_fit_params - 1
    dof_noZeroBySurvey = len(all_sn) - n_fit_params - 1 #must account for overall shifts from 0 as fit parameter 
    
    periodigram_comparator = ComparatorBetweenFunctionAndDistributionOfZero(limits_file, '', fitting_funct_type = 'alt_normal')
    canon_H0 = cosmo_arch.getH0()
    #min_H0_scaling = (canon_H0[0] - 5.0 * canon_H0[1]) / canon_H0[0]
    #max_H0_scaling = (canon_H0[0] + 5.0 * canon_H0[1]) / canon_H0[0] 
    for i in range(len(Cs_to_calc)):
        macro_element_start = time.time() 
        C = Cs_to_calc[i]
        print 'Working on i ' + str(i) + ' of ' + str(len(Cs_to_calc)) + ' at which C = ' + str(C)
        for j in range(len(psis_to_calc)):
            element_start = time.time()
            psi = psis_to_calc[j]
            print 'Working on j ' + str(j) + ' of ' + str(len(psis_to_calc)) + ' at which psi = ' + str(psi)
            print 'Measuring distance modulus residuals... '
            new_resid_calc = ResidualMuCalculatorForOscillatingG (A, nu, alpha, psi, C, filled_zs)
            #raw_theoretical_mu_resids = function_to_limit(sorted_lookBacks, new_resid_calc C, psi)
            raw_theoretical_mu_resids = function_to_limit(sorted_lookBacks, new_resid_calc)
            mean_theoretical_mu_resids = np.mean(raw_theoretical_mu_resids)
            theoretical_mu_resids_noZBS = [resid - mean_theoretical_mu_resids for resid in raw_theoretical_mu_resids]
            new_dG_over_G = new_resid_calc.d_G_d_t_of_z[0]
            new_ddG_over_G = new_resid_calc.dd_G_dt2_0 
            print 'new_dG_over_G = ' + str(new_dG_over_G)
            print 'new_ddG_over_G = ' + str(new_ddG_over_G) 
            #plt.scatter(sorted_zs, theoretical_mu_resids_noZBS)
            #plt.title('Oscillating Gravity From Phi C = ' + str(C) + ' psi = ' + str(psi) )
            #plt.savefig('/Users/sasha/Desktop/muResidsOscillatingG1/' + 'OcsillatingG_C' + str(abs(C)) + '_psi_' + str(psi) + '_alpha_' + alpha_str + '_nu_' + str(int(nu)) + '.png')
            #plt.close('all') 
            #plt.show() 

            print 'Computing reduced chi square... '
            new_chi_sqr_noZBS = sum([((theoretical_mu_resids_noZBS[k] - sorted_resids_noBySurveyZeroing[k]) / (sorted_muErrs[k])) ** 2.0 for k in range(len(sorted_zs))]) / dof_noZeroBySurvey
            
            #init_peak_guess = np.argmin(np.abs(np.array(periodigram_comparator.frequencies) - omega_in_Planck_Myr))
            #freq_bounds_over_which_to_find_peak = [max(0, init_peak_guess - freq_half_width_to_find_peak),
            #                                       min(len(periodigram_comparator.frequencies), init_peak_guess + freq_half_width_to_find_peak)]
            #print 'init_peak_guess = ' + str(init_peak_guess)
            #print 'freq_bounds_over_which_to_find_peak = ' + str(freq_bounds_over_which_to_find_peak)
            #start_freqs_of_interest = [0, min(100, freq_bounds_over_which_to_find_peak[0])] #Freqs to check since largest Fourier mode can sometimes be in early distribution
            freq_bounds_over_which_to_find_peak = None 

            print 'Computing normalized periodogram...'
            normalized_coefs_noZBS, CPB_vals_of_zero_distributions_noZBS = periodigram_comparator.computeLikelihoodsOfSignalFromZeros(sorted_lookBacks,
                                                                                                                                                                theoretical_mu_resids_noZBS,
                                                                                                                                                                sorted_muErrs,
                                                                                                                                                                freq_bounds_over_which_to_find_peak = freq_bounds_over_which_to_find_peak)
                
            #normalized_coefs_noZBS = normalized_coefs_at_start_freqs_noZBS + normalized_coefs_around_guess_noZBS
            #CPB_vals_of_zero_distributions_noZBS  = CPB_vals_of_zero_distributions_at_start_freqs_noZBS + CPB_vals_of_zero_distributions_around_guess_noZBS 

            peak_CPB_val_noZBS = 10.0 ** max(np.log10(CPB_vals_of_zero_distributions_noZBS))
            peak_CPB_val_noZBS_index = np.argmax(np.log10(CPB_vals_of_zero_distributions_noZBS))
            #print 'peak_CPB_val_noZBS_index = ' + str(peak_CPB_val_noZBS_index)
            calc_chi_sqrs_noZBS[i][j] = new_chi_sqr_noZBS
            peak_CPB_vals_noZBS[i][j] = peak_CPB_val_noZBS
            peak_CPB_noZBS_index[i][j] = peak_CPB_val_noZBS_index
            dG_over_Gs[i][j] = new_dG_over_G 
            ddG_over_Gs[i][j] = new_ddG_over_G 
            #peak_CPB_vals_yesZBS[i][j] = peak_CPB_val_yesZBS
            element_end = time.time()
            if j % 1 == 0:
                print 'j ' + str(j) + ' of ' + str(len(psis_to_calc)) + ' took ' + str(element_end - element_start) + 's'
                #print '    new_H0_scaling_first_guess * np.array((new_H0_min_bound_coefs[0], 1.0, new_H0_min_bound_coefs[1])) = ' + str(new_H0_scaling_first_guess  * np.array((new_H0_min_bound_coefs[0], 1.0, new_H0_min_bound_coefs[1]))) 
                #print '    new_H0_scaling = ' + str(new_H0_scaling)
                #print '    This quantity should be (nearly) 0.0: ' + str(function_to_constrain_at_H0(A, omega_in_Myr / (new_H0_scaling * H0_in_inv_Myr), H0_in_inv_Myr * new_H0_scaling))
                #print '    omega_in_Myr / (new_H0_scaling * H0_in_inv_Myr) = ' + str(omega_in_Myr / (new_H0_scaling * H0_in_inv_Myr))
                #print '    omegas_to_calc_in_H0[j]  = ' + str(omegas_to_calc_in_H0[j])
                #print '    peak_CPB_val = ' + str(peak_CPB_val)
                #print '    peak_CPB_val_index = ' + str(peak_CPB_val_index)
 
        macro_element_end = time.time() 
        print 'i ' + str(i) + ' of ' + str(len(Cs_to_calc)) + ' took  ' + str(macro_element_end - macro_element_start) + 's'
    #plt.scatter(sorted_zs, theoretical_mu_resids)
    #plt.scatter(sorted_zs, sorted_resids)
    #plt.errorbar(sorted_zs, sorted_resids, yerr = sorted_muErrs, fmt = None)
    #plt.show()
    res_limites_from_fourier = []

    #file_name = 'ForMonodromicType1_alpha' + str(np.around(phase_pi_coef, 3)) + '+' + str(int(phase_shift_int)) + '_' + str(phase_shift_order) +  'withFourier' + '_correctedDOF_final1'
    file_name = 'ForGofPhi_alpha_' + str(int(np.around(alpha * 10.0, 0))) + '_10' + '_A_' + str(int(A)) + '_nu_' + str(int(nu)) + 'withFourier' + '_correctedDOF_corrected_theta0_final1'
    #np.save('/Users/sasha/Desktop/RChiSqrLimits/' + 'rChiSqr_yesZBS' + file_name, calc_chi_sqrs_yesZBS)
    np.save('/Users/sasha/Desktop/RChiSqrLimits/' + 'rChiSqr_noZBS' + file_name, calc_chi_sqrs_noZBS)
    #np.save('/Users/sasha/Desktop/RChiSqrLimits/' + 'peakCPBVals_yesZBS' + file_name, peak_CPB_vals_yesZBS)
    np.save('/Users/sasha/Desktop/RChiSqrLimits/' + 'peakCPBVals_noZBS' + file_name, peak_CPB_vals_noZBS)
    np.save('/Users/sasha/Desktop/RChiSqrLimits/' + 'peakCPBVals_noZBS_index' + file_name, peak_CPB_noZBS_index)
    np.save('/Users/sasha/Desktop/RChiSqrLimits/' + 'dG_over_G_vals' + file_name, dG_over_Gs)
    np.save('/Users/sasha/Desktop/RChiSqrLimits/' + 'ddG_over_G_vals' + file_name, ddG_over_Gs) 

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
    chi_sqr_noZBS_contours = plt.contourf(C_mesh, psi_mesh, calc_chi_sqrs_noZBS, levels = threshold_r_chi_sqr_vals_noZeroBySurvey , colors = ['yellow','g','b','m','r','orange'], alpha = 0.5)
    plt.xticks(np.linspace(0.0, math.pi, 5), ['0', r'$\pi$' + r'$\mathrm{/4}$', r'$\pi$' + r'$\mathrm{/2}$', r'$\mathrm{3}$' + r'$\pi$' + r'$\mathrm{/4}$', '$\pi$'] )
    plt.xlim(min(psis_to_calc), max(psis_to_calc)) 
    label_size = 20.0
    plt.xlabel(r'$\psi_{G}$' + '(unitless)', fontsize = label_size, labelpad = label_size / 2.0)
    plt.ylabel('$C$' + r'$_{G}$' + '(unitless)', fontsize = label_size)
    #plt.text((np.max(omegas_to_calc_in_H0) + np.mean(omegas_to_calc_in_H0))/2.5, np.mean(As_to_calc), r'$\phi_X$' + ' = ' + phase_str, fontsize = label_size )
    plt.text((np.max(psis_to_calc) + np.mean(psis_to_calc))/2.5, np.mean(Cs_to_calc), r'$\nu$' + r'$_G$' + ' = ' + nu_str, fontsize = label_size )
    plt.title('',  fontsize = label_size)
    cbar = plt.colorbar(chi_sqr_noZBS_contours) 
    cbar.ax.set_yticklabels([''] + ['$R_{\mathrm{rej}}$ = ' + str(1.0 - np.around(prob, 5)) for prob in target_probs])
    plt.ylim([min(Cs_to_calc), max(Cs_to_calc)])
    plt.tight_layout()
    plt.savefig('/Users/sasha/Desktop/RChiSqrLimits/' + 'chSqrVals_noZBS' + file_name + '.png')
    plt.close('all')

    #
    #contours = plt.contourf(A_mesh, omega_mesh, (1.0 - peak_CPB_vals_noZBS ** N_peaks) / ref_fourier_CPB_prob,
    #                        levels = list(reversed(target_relative_fourier_probs)), colors = list(reversed(['yellow', 'g','b','m','r','orange'])), alpha = 0.5)
    print 'target_relative_fourier_probs = ' + str(target_relative_fourier_probs) 
    #Note: in these contours, I have to multiply everything by -1 so that levels are increasing.  plt.contourf seems to only work if contours are increasing.  
    contours = plt.contourf(C_mesh, psi_mesh, -1.0 * (1.0 - peak_CPB_vals_noZBS ** N_peaks) / ref_fourier_CPB_prob,
                            levels = [-1.0 * prob for prob in target_relative_fourier_probs], colors = ['yellow', 'g','b','m','r','orange'], alpha = 0.5)  
    plt.xticks(np.linspace(0.0, math.pi, 5), ['0', r'$\pi$' + r'$\mathrm{/4}$', r'$\pi$' + r'$\mathrm{/2}$', r'$\mathrm{3}$' + r'$\pi$' + r'$\mathrm{/4}$', '$\pi$'] )
    plt.yticks([0.005, 0.01, 0.015, 0.02])
    label_size = 20.0
    plt.xlabel(r'$\psi$' + r'$_G$' + ' (unitless)', fontsize = label_size, labelpad = label_size / 2.0)
    plt.ylabel('$C$' + r'$_G$' + '(unitless)', fontsize = label_size)
    #plt.text((np.max(omegas_to_calc_in_H0) + np.mean(omegas_to_calc_in_H0))/2.5, np.mean(As_to_calc), r'$\phi_X$' + ' = ' + phase_str, fontsize = label_size )
    plt.text((np.max(psis_to_calc) + np.mean(psis_to_calc))/2.5, np.mean(Cs_to_calc), r'$\nu$' + r'$_G$' + ' = ' + nu_str, fontsize = label_size )
    plt.title('',  fontsize = label_size)
    cbar = plt.colorbar(contours)
    cbar.ax.set_yticklabels([''] + [r'$R_{\mathrm{rej}}$ = ' + str(np.around(1.0 - prob, 5)) for prob in target_probs])
    plt.ylim([min(Cs_to_calc), max(Cs_to_calc)])
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
    


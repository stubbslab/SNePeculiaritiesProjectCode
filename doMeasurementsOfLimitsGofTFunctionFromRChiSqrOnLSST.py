import math
import numpy as np
from measureLimitsOfFunctionGivenDistributionOfZero import measureLimitsOfFunctionGivenDistributionOfZero
from loadSN import loadSN
import scipy.optimize as optimize 
from AstronomicalParameterArchive import AstronomicalParameterArchive
from CosmologicalParameterArchive import CosmologicalParameterArchive
from calculateMuForGOfTFromz import ResidualMuCalculatorForArbitraryGofT
import matplotlib.pyplot as plt
from StandardPeriodigram import StandardPeriodigram
import time

if __name__ == '__main__':
    #Number of sigma
    start = time.time() 
    r_chi_sqr_target = 5.0  
    
    #Must be a two parameter function
    # (ideally one is like a frequency and one is like an amplitude)

    #Start with a simple oscillitory function
    omega0 = 0.05
    bounds = [0.0, 100.0]
    t_or_tau = 'tau'
    if t_or_tau.lower() in ['tau','taus']: 
        GOfTTauZ = lambda ts, taus, zs, A, omega, phi: 1.0 + A * np.sin(omega * np.array(taus) + phi)
    else:
        GOfTTauZ = lambda ts, taus, zs, A, omega, phi: 1.0 + A * np.sin(omega * np.array(ts) + phi)
    
    #file_of_zero_limits = '/Users/sasha/Documents/Harvard/physics/stubbs/SNIsotropyProject/randomBestFitParams/params_from_300_frequency_scaling50_10_sampling_of_fitted_alt_gaussian_to_zero.csv'

    #results_for_minimization = np.genfromtxt(file_of_zero_limits, delimiter = ',').transpose()
    #frequencies = results_for_minimization[0]
    
    H0_in_inv_Myr = 6.883156370706416 * 10 ** (-5.0)
    #freq_like_params_on_which_to_minimize = [frequencies[5 + i * 100] / H0_in_inv_Mpc for i in range(len(frequencies[0:1000]) / 100)]
    freq_like_params_on_which_to_minimize = [0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 200.0, 800.0]
    #freq_indeces_for_finding_peak = [[max(0, 5 + i * 50 + freq_index_range_for_finding_peak[0]),
    #                                  min(len(frequencies), 5 + i * 50 + freq_index_range_for_finding_peak[1])]
    #                                 for i in range(len(frequencies) / 50)]
    

    astro_arch = AstronomicalParameterArchive()
    s_to_yr = astro_arch.getSecondToYear() 
    cosmo_archive = CosmologicalParameterArchive()
    age_of_universe = cosmo_archive.getAgeOfUniverse( units = 'yr')[0]
    
    all_sn = loadSN(1, pull_extinctions = 0, data_type = 'art_lsst')
    all_surveys = np.unique([sn['survey'] for sn in all_sn])
    all_zs = [sn['z'] for sn in all_sn]
    all_ts = [sn['t'] for sn in all_sn]
    all_taus = [sn['tau'] for sn in all_sn]
    all_muErrs = [sn['muErr'] for sn in all_sn]
    all_muResids = [sn['muDiff'] for sn in all_sn]
    all_zerodMuResids = [sn['muDiff'] - sn['muDiffWMean'] for sn in all_sn]
    #sorted_zs = sorted(all_zs)
    if t_or_tau.lower() in ['tau','taus']:
        all_lookBacks = [10 ** (-6.0) * tau * s_to_yr for tau in all_taus]
    else:
        all_lookBacks = [10 ** (-6.0) * t * s_to_yr for t in all_ts]
    #sorted_lookBacks = [t_or_tau for _,t_or_tau in sorted(zip(all_zs, all_lookBacks))]
    #sorted_muErrs = [err for _,err in sorted(zip(all_zs, all_muErrs))]
    #sorted_zResids = [res for _,res in sorted(zip(all_zs, all_zerodMuResids))]
    subtract_individual_means = 1
    if subtract_individual_means:
        sorted_zs, sorted_lookBacks, sorted_muErrs, sorted_resids = safeSortOneListByAnother(all_zs, [all_zs, all_lookBacks, all_muErrs, all_zerodMuResids])
        #sorted_resids = [res for _,res in sorted(zip(all_zs, all_zerodMuResids))]
    else:
        sorted_zs, sorted_lookBacks, sorted_muErrs, sorted_resids = safeSortOneListByAnother(all_zs, [all_zs, all_lookBacks, all_muErrs, all_Resids])
        #sorted_resids = [res for _,res in sorted(zip(all_zs, all_muResids))]   
    canon_mus = ResidualMuCalculatorForArbitraryGofT(Gof = lambda ts, taus, zs: 1.0 if type(ts) in [int, float] else 1.0 + np.zeros(np.shape(ts)), initial_zs = sorted_zs).getMus()

    
    asymptotic_z = 1000.0
    phase_pi_coef = 0.5
    phase_shift_int = -1.0
    phase_shift_order = 100
    phase = math.pi * phase_pi_coef + phase_shift_int / phase_shift_order
    change_SN_phys = 1
    extended_zs = np.linspace(0.0, asymptotic_z, 2000)
    asymptotic_canon_mus = ResidualMuCalculatorForArbitraryGofT(Gof = lambda ts, taus, zs: 1.0 if type(ts) in [int, float] else 1.0 + np.zeros(np.shape(ts)),
                                                                initial_zs = extended_zs, change_SN_phys = change_SN_phys).getMus()
    

    #I think that we can define the function to just take in sorted ts, 
    # which I know maps one to one with sorted zs, and then return the calculation with the zs.
    
    function_to_limit = lambda sorted_lookBacks, Amp_like, freq_like, new_H0_scaling: (np.array(ResidualMuCalculatorForArbitraryGofT(Gof = lambda ts, taus, zs:
                                                                                                                                     GOfTTauZ(ts, taus, zs, Amp_like, freq_like, phase),
                                                                                                                                     initial_zs = sorted_zs, H0 = H0_in_inv_Myr * new_H0_scaling,
                                                                                                                                     change_SN_phys = change_SN_phys).getMus())
                                                                       -  np.array(canon_mus))
    function_to_constrain_at_H0 = lambda Amp_like, freq_like, new_H0: abs(np.array(ResidualMuCalculatorForArbitraryGofT(Gof = lambda ts, taus, zs:
                                                                                                                        GOfTTauZ(ts, taus, zs, Amp_like, freq_like, phase),
                                                                                                                        initial_zs = extended_zs, H0 = new_H0,
                                                                                                                        change_SN_phys = change_SN_phys).getMus()[-1]
                                                                                                  )
                                                                                         -  np.array(asymptotic_canon_mus)[-1])

    minimized_values = []
    #Calculate the residual over a double for loop, iterating through both the allowed omegas and allowed As
    omegas_to_calc_in_H0 = np.linspace(0.0, 200.0, 250)
    omegas_to_calc_in_Myr = omegas_to_calc_in_H0 * H0_in_inv_Myr
    As_to_calc = np.linspace(0.0, 0.5, 200)
    A_mesh, omega_mesh = np.meshgrid(omegas_to_calc_in_Myr, As_to_calc)
    calc_chi_sqrs = np.zeros(np.shape(A_mesh))
    for i in range(len(As_to_calc)):
        print 'Working on i ' + str(i) + ' of ' + str(len(As_to_calc))
        A = As_to_calc[i]
        for j in range(len(omegas_to_calc_in_Myr)):
            element_start = time.time() 
            omega_in_Myr = omegas_to_calc_in_Myr[j]
            new_H0_scaling = 1.0 
            #I now need to find an H0 for our fit that brings cosmologies into asymptotic agreement (agreement at z = 1000)
            #new_H0_scaling = optimize.minimize_scalar(lambda new_H0_scaling:  function_to_constrain_at_H0(A, omega_in_Myr / (new_H0_scaling * H0_in_inv_Myr),
            #                                                                                              H0_in_inv_Myr * new_H0_scaling),
            #                                          bounds = [0.1, 10.0], method = 'Bounded')['x']
            new_H0 = H0_in_inv_Myr * new_H0_scaling
            
            theoretical_mu_resids = function_to_limit(sorted_lookBacks, A, omega_in_Myr / (new_H0_scaling * H0_in_inv_Myr), new_H0_scaling)
            new_chi_sqr = sum([((theoretical_mu_resids[k] - sorted_resids[k]) / (sorted_muErrs[k])) ** 2.0 for k in range(len(sorted_zs))]) / len(sorted_resids)
            calc_chi_sqrs[i][j] = new_chi_sqr
            element_end = time.time()
            #if j % 50 == 0:
            #    print 'j ' + str(j) + ' of ' + str(len(omegas_to_calc_in_Myr)) + ' took ' + str(element_end - element_start) + 's'
            #    print '    new_H0_scaling = ' + str(new_H0_scaling)
            #    print '    omega_in_Myr / (new_H0_scaling * H0_in_inv_Myr) = ' + str(omega_in_Myr / (new_H0_scaling * H0_in_inv_Myr))
            #    print '    omega_in_Myr = ' + str(omega_in_Myr)
            #    print '    omegas_to_calc_in_H0[j]  = ' + str(omegas_to_calc_in_H0[j])
    #plt.scatter(sorted_zs, theoretical_mu_resids)
    #plt.scatter(sorted_zs, sorted_resids)
    #plt.errorbar(sorted_zs, sorted_resids, yerr = sorted_muErrs, fmt = None)
    #plt.show() 

    print 'calc_chi_sqrs = '
    print calc_chi_sqrs
    
    #threshold_r_chi_sqr = 1.22642
    n_fit_params = 3
    threshold_r_chi_sqr = 1.0 + 4.865 * math.sqrt(2.0 / (len(all_sn) - n_fit_params))
    min_r_chi_sqr = 0.95

    #According to Varitions of the Gravitational Constant from Lunar Laser Ranging Data 
    n_sigma_for_existing_G_lims = 5.0
    G_dot_G_lims_inv_Myr = [abs((2.0 - 7.0 * n_sigma_for_existing_G_lims ) * 10.0 ** -13.0 * 10.0 ** 6.0), abs((2.0 + 7.0 * n_sigma_for_existing_G_lims )  * 10.0 ** -13.0 * 10.0 ** 6.0)]
    G_ddot_G_lims_inv_Myr_sqr = [abs((4.0 - 5.0 * n_sigma_for_existing_G_lims ) * 10.0 ** -15.0 * 10.0 ** 12.0), abs((4.0 + 5.0 * n_sigma_for_existing_G_lims ) * 10.0 ** -15.0 * 10.0 ** 12.0) ]
    if math.cos(phase) < 0.0:
        limit_curve_from_G_dot = lambda omegas: [- G_dot_G_lims_inv_Myr[0] / (omega *  math.cos(phase) )  for omega in omegas]
    else:
        limit_curve_from_G_dot = lambda omegas: [G_dot_G_lims_inv_Myr[1] / (omega *  math.cos(phase) )  for omega in omegas]
    if math.sin(phase) < 0.0: 
        limit_curve_from_G_ddot = lambda omegas: [G_ddot_G_lims_inv_Myr_sqr[0] / (omega *  math.sin(phase) )  for omega in omegas]
    else:
        limit_curve_from_G_ddot = lambda omegas: [G_ddot_G_lims_inv_Myr_sqr[1] / (omega *  math.sin(phase) ) for omega in omegas]
    
    levels = np.linspace(min_r_chi_sqr, threshold_r_chi_sqr, 100)
    contours = plt.contourf(A_mesh, omega_mesh, calc_chi_sqrs, levels = levels)
    plt.xlabel(r'$\omega$' + ' (1/Myr)')
    plt.ylabel('A')
    plt.title('')
    file_name = 'rChiSqrForGSin_phi_pix' + str(np.around(phase_pi_coef, 3)) + '+' + str(int(phase_shift_int)) + '_' + str(phase_shift_order) + '_SNPhys' + str(change_SN_phys) + '_2.png'
    cbar = plt.colorbar(contours)

    
    G_dot_lim = plt.plot(omegas_to_calc_in_Myr, limit_curve_from_G_dot(omegas_to_calc_in_Myr), c = 'c')
    #G_ddot_lim = plt.plot(omegas_to_calc_in_Myr, limit_curve_from_G_ddot(omegas_to_calc_in_Myr), c = 'm')
    #print 'omegas_to_calc_in_Myr = ' + str(omegas_to_calc_in_Myr)
    #print 'limit_curve_from_G_dot(omegas_to_calc_in_Myr) = ' + str(limit_curve_from_G_dot(omegas_to_calc_in_Myr)) 
    #print 'limit_curve_from_G_ddot(omegas_to_calc_in_Myr) = ' + str(limit_curve_from_G_ddot(omegas_to_calc_in_Myr))
    plt.ylim([min(As_to_calc), max(As_to_calc)])
    #plt.legend([ G_dot_lim[0], G_ddot_lim[0]],
    #           [str(n_sigma_for_existing_G_lims) + r'$\sigma$' + ' limits from ' + r'$\dot{G}$' + '/' + r'$G$',str(n_sigma_for_existing_G_lims) + r'$\sigma$' + ' limits from ' + r'$\ddot{G}$' + '/' + r'$G$'])
    plt.legend([ G_dot_lim[0]], [str(n_sigma_for_existing_G_lims) + r'$\sigma$' + ' limits from ' + r'$\dot{G}$' + '/' + r'$G$'])
    
    plt.savefig('/Users/sasha/Desktop/' + file_name)
    end = time.time()
    print 'Took ' + str(end - start)  + 's to do computation.'  
    #plt.show()


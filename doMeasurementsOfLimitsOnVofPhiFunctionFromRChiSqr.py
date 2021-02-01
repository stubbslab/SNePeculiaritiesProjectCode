import math
import numpy as np
from measureLimitsOfFunctionGivenDistributionOfZero import measureLimitsOfFunctionGivenDistributionOfZero
from loadSN import loadSN
import scipy.optimize as optimize 
from AstronomicalParameterArchive import AstronomicalParameterArchive
from CosmologicalParameterArchive import CosmologicalParameterArchive
from calculateMuForwOftFromz import ResidualMuCalculatorForArbitraryWofT
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
    VofPhi = lambda phis, A, alpha, nu, C, phi_0: C * (np.array(phis)/phi_0) ** (-alpha) * (1 - A*np.sin(np.array(phis)*nu))
    dVofPhi = lambda phis, A, alpha, nu, C, phi_0: (C * (-alpha) * 1.0 / phi_0 * (np.array(phis)/phi_0) ** (-alpha-1.0) * (1 - A*np.sin(np.array(phis)*nu))
                                                            + C * (np.array(phis)/phi_0) ** (-alpha) * (- A * nu * np.cos(np.array(phis)*nu)) )
    
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
    subtract_individual_means = 1
    if subtract_individual_means:
        sorted_resids = [res for _,res in sorted(zip(all_zs, all_zerodMuResids))]
    else:
        sorted_resids = [res for _,res in sorted(zip(all_zs, all_muResids))]   

    canon_mus = ResidualMuCalculatorForArbitraryWofT(wOfFunction = lambda ts, taus, zs: -1.0 if type(ts) in [int, float] else -1.0 + np.zeros(np.shape(ts)), initial_zs = sorted_zs).getMus()

    
    asymptotic_z = 1000.0
    alpha_fixed = 0.2
    C_fixed = 1.0
    phi_0_fixed = 1.0
    theta_0_fixed = 0.0 
    extended_zs = np.linspace(0.0, asymptotic_z, int(asymptotic_z * 10))
    asymptotic_canon_mus = ResidualMuCalculatorForArbitraryWofT(wOfFunction = lambda ts, taus, zs: -1.0 if type(ts) in [int, float] else -1.0 + np.zeros(np.shape(ts)), initial_zs = extended_zs).getMus()
    

    #I think that we can define the function to just take in sorted ts, 
    # which I know maps one to one with sorted zs, and then return the calculation with the zs.
    
    function_to_limit = lambda sorted_lookBacks, Amp_like, freq_like, new_H0_scaling: (np.array(ResidualMuCalculatorForArbitraryWofT(VOfPhiFunction = lambda phis: VofPhi(phis, Amp_like,
                                                                                                                                                                  alpha_fixed, freq_like,
                                                                                                                                                                  C_fixed, phi_0_fixed),
                                                                                                                                     dVOfPhiFunction = lambda phis: dVofPhi(phis, Amp_like,
                                                                                                                                                                  alpha_fixed, freq_like,
                                                                                                                                                                  C_fixed, phi_0_fixed),
                                                                                                                                     initial_zs = sorted_zs, theta0 = theta_0_fixed, 
                                                                                                                                     H0 = H0_in_inv_Myr * new_H0_scaling, wFromPhi = 1).getMus())
                                                                                       -  np.array(canon_mus))
    function_to_constrain_at_H0 = lambda Amp_like, freq_like, new_H0: abs(np.array(ResidualMuCalculatorForArbitraryWofT(VOfPhiFunction = lambda phis: VofPhi(phis, Amp_like,
                                                                                                                                                    alpha_fixed, freq_like,
                                                                                                                                                    C_fixed, phi_0_fixed),
                                                                                                                        dVOfPhiFunction = lambda phis: dVofPhi(phis, Amp_like,
                                                                                                                                                                  alpha_fixed, freq_like,
                                                                                                                                                                  C_fixed, phi_0_fixed),
                                                                                                                        initial_zs = extended_zs, H0 = new_H0, theta0 = theta_0_fixed, wFromPhi = 1).getMus()[-1]
                                                                                                  )
                                                                                         -  np.array(asymptotic_canon_mus)[-1])


    minimized_values = []
    #Calculate the residual over a double for loop, iterating through both the allowed nus and allowed As
    nus_to_calc_in_H0 = np.linspace(0.0, 100.0 ,200)
    nus_to_calc_in_Myr = nus_to_calc_in_H0 * 6.883156370706416 * 10 ** (-5.0)
    As_to_calc = np.linspace(0.0, 0.5, 150)
    A_mesh, nu_mesh = np.meshgrid(nus_to_calc_in_Myr, As_to_calc)
    calc_chi_sqrs = np.zeros(np.shape(A_mesh))
    print 'np.shape(calc_chi_sqrs)= ' + str(np.shape(calc_chi_sqrs)) 
    for i in range(len(As_to_calc)):
        print 'Working on i ' + str(i) + ' of ' + str(len(As_to_calc))
        A = As_to_calc[i]
        for j in range(len(nus_to_calc_in_Myr)):
            element_start = time.time() 
            nu_in_Myr = nus_to_calc_in_Myr[j]
            #print 'j = ' + str(j) 
            #print '    nus_to_calc_in_H0[j]  = ' + str(nus_to_calc_in_H0[j])
            #print 'A = ' + str(A) 
            
            #I now need to find an H0 for our fit that brings cosmologies into asymptotic agreement (agreement at z = 1000)
            
            new_H0_scaling = optimize.minimize_scalar(lambda new_H0_scaling:  function_to_constrain_at_H0(A, nu_in_Myr / (new_H0_scaling * H0_in_inv_Myr),
                                                                                                          H0_in_inv_Myr * new_H0_scaling),
                                                      bounds = [0.1, 10.0], method = 'Bounded')['x']
            
            #print '    new_H0_scaling = ' + str(new_H0_scaling)
            #print '    nu_in_Myr / (new_H0_scaling * H0_in_inv_Myr) = ' + str(nu_in_Myr / (new_H0_scaling * H0_in_inv_Myr))
            #new_H0_scaling = 1.0
            new_H0 = H0_in_inv_Myr * new_H0_scaling
            
            theoretical_mu_resids = function_to_limit(sorted_lookBacks, A, nu_in_Myr / (new_H0_scaling * H0_in_inv_Myr), new_H0_scaling)
            new_chi_sqr = sum([((theoretical_mu_resids[k] - sorted_resids[k]) / (sorted_muErrs[k])) ** 2.0 for k in range(len(sorted_zs))]) / len(sorted_resids)
            calc_chi_sqrs[i][j] = new_chi_sqr
            element_end = time.time()
            if j % 10 == 1:
                print 'j ' + str(j) + ' of ' + str(len(nus_to_calc_in_Myr)) + ' took ' + str(element_end - element_start) + 's'
                print '    new_H0_scaling = ' + str(new_H0_scaling)
                print '    nu_in_Myr / (new_H0_scaling * H0_in_inv_Myr) = ' + str(nu_in_Myr / (new_H0_scaling * H0_in_inv_Myr))
                print '    nus_to_calc_in_H0[j]  = ' + str(nus_to_calc_in_H0[j])
                print '    A = ' + str(A) 
    #plt.scatter(sorted_zs, theoretical_mu_resids)
    #plt.scatter(sorted_zs, sorted_resids)
    #plt.errorbar(sorted_zs, sorted_resids, yerr = sorted_muErrs, fmt = None)
    #plt.show() 

    print 'calc_chi_sqrs = '
    print calc_chi_sqrs

    threshold_r_chi_sqr = 1.22667
    min_r_chi_sqr = 0.98
    levels = np.linspace(min_r_chi_sqr, threshold_r_chi_sqr, 100)
    contours = plt.contourf(A_mesh, nu_mesh, calc_chi_sqrs, levels = levels)
    plt.xlabel(r'$\omega$' + ' (1/Myr)')
    plt.ylabel('A')
    plt.title('')
    file_name = 'rChiSqrForVPhiSin_alpha2_10_2.png'
    cbar = plt.colorbar(contours)
    plt.savefig('/Users/sasha/Desktop/' + file_name)
    end = time.time()
    print 'Took ' + str(end - start)  + 's to do computation.'  
    #plt.show()


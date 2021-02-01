import numpy as np
import csv 
from scipy.interpolate import interp1d
from loadSN import loadSN
from AstronomicalParameterArchive import AstronomicalParameterArchive
from CosmologicalParameterArchive import CosmologicalParameterArchive 
from calculateMuForDaDtofTFromz import ResidualMuCalculatorForSpecifiedDaDtofT
from calculateMuForXOfTauFromz import ResidualMuCalculatorForArbitraryXofT
import matplotlib.pyplot as plt
from cantrips import safeSortOneListByAnother
from measureLimitsOfFunctionGivenDistributionOfZero import ComparatorBetweenFunctionAndDistributionOfZero

if __name__  == "__main__":

    astro_arch = AstronomicalParameterArchive()
    cosmo_arch = CosmologicalParameterArchive()
    s_to_yr = astro_arch.getSecondToYear() 
    cosmo_archive = CosmologicalParameterArchive()
    age_of_universe = cosmo_archive.getAgeOfUniverse( units = 'yr')[0]
    t_or_tau = 'tau' 
    H0_in_inv_Myr = 6.883156370706416 * 10 ** (-5.0)
    
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
    subtract_individual_means = 1
    sorted_zs, sorted_lookBacks, sorted_muErrs, sorted_resids, sorted_resids_noBySurveyZeroing, sorted_surveys = safeSortOneListByAnother(all_zs, [all_zs, all_lookBacks, all_muErrs, all_zerodMuResids, all_flatZerodMuResids, all_surveys])
    filled_zs = np.linspace(10.0 ** -10.0, max(sorted_zs), 10001)

    #canon_mus_filled = ResidualMuCalculatorForSpecifiedDaDtofT(DaDtPerturbationFunction = lambda ts, taus, zs: 0.0 if type(ts) in [float, int] else 0.0 + np.zeros(np.shape(ts)), initial_zs = np.array([0.0] + filled_zs.tolist())).getMus()
    canon_mus = interp1d(filled_zs, np.array(ResidualMuCalculatorForSpecifiedDaDtofT(DaDtPerturbationFunction = lambda ts, taus, zs: 0.0 if type(ts) in [float, int] else 0.0 + np.zeros(np.shape(ts)), initial_zs = filled_zs).getMus()))(sorted_zs)

    target_dir = '/Users/sasha/Documents/Harvard/physics/stubbs/SNIsotropyProject/'
    file_of_resid_da_of_t = 'my_approx_of_t_vs_aDot_of_Ringermacher_1.csv'
    rows = []
    with open(target_dir + file_of_resid_da_of_t) as csv_file:
        csv_reader = csv.DictReader(csv_file)
        rows = []
        for row in csv_reader: 
            rows = rows + [[float(row['x']), float(row['Curve1'])]] 

    interp_ts = [row[0] for row in rows]
    interp_da_resids = [row[1] for row in rows]
    print ('[interp_ts,interp_da_resids] =' + str([interp_ts,interp_da_resids]) )

    da_of_t_interp = interp1d(interp_ts, interp_da_resids)

    min_t = min(interp_ts) 
    min_interp_t = min(interp_ts) 
    max_interp_t = max(interp_ts) 
    extended_da_of_t = lambda ts: (da_of_t_interp(ts) if (ts >= min_interp_t and ts <= max_interp_t) else 0.0 ) if type(ts) in [int, float, np.float64] else [da_of_t_interp(t) if (t >= min_interp_t and t <= max_interp_t) else 0.0 for t in ts ]

    limits_file = '/Users/sasha/Documents/Harvard/physics/stubbs/SNIsotropyProject/randomBestFitParams/params_from_1000_frequency_scaling50_10_max_oscillations_1000_sampling_of_fitted_tauexponential_to_zero_updated_SN_subtract_wMean.csv'
    n_fit_params = 5 
    dof_zeroBySurvey = len(all_sn) - n_fit_params - len(unique_surveys) #must account for individual shifts from 0 as fit parameters
    #dof_zeroBySurvey = len(all_sn) - n_fit_params - 1
    dof_noZeroBySurvey = len(all_sn) - n_fit_params - 1 #must account for overall shifts from 0 as fit parameter 
    
    periodigram_comparator = ComparatorBetweenFunctionAndDistributionOfZero(limits_file, '', fitting_funct_type = 'alt_normal')
    
    function_to_limit = lambda sorted_lookBacks, resid_mu_calculator: ( interp1d(filled_zs, np.array(resid_mu_calculator.getMus()))(sorted_zs) 
                                                                        - np.array(canon_mus) )

    new_resid_calc = ResidualMuCalculatorForSpecifiedDaDtofT(DaDtPerturbationFunction = lambda ts, taus, zs: extended_da_of_t(ts), initial_zs = filled_zs)
    raw_theoretical_mu_resids = function_to_limit(sorted_lookBacks, new_resid_calc)
    mean_theoretical_mu_resids = np.mean(raw_theoretical_mu_resids)
    theoretical_mu_resids_noZBS = [resid - mean_theoretical_mu_resids for resid in raw_theoretical_mu_resids]
    #plt.scatter(sorted_zs, theoretical_mu_resids_noZBS)
    #plt.xlabel(r'$z$')
    #plt.ylabel(r'$\Delta\mu$' + ' (mean subtracted)')
    #plt.title('Predicted distance modulus residuals') 
    #plt.show() 
    new_chi_sqr_noZBS = sum([((theoretical_mu_resids_noZBS[k] - sorted_resids_noBySurveyZeroing[k]) / (sorted_muErrs[k])) ** 2.0 for k in range(len(sorted_zs))]) / dof_noZeroBySurvey

    freq_bounds_over_which_to_find_peak = None 

    print ('Computing normalized periodogram...')
    frequencies = periodigram_comparator.frequencies 
    normalized_coefs_noZBS, CPB_vals_of_zero_distributions_noZBS = periodigram_comparator.computeLikelihoodsOfSignalFromZeros(sorted_lookBacks, theoretical_mu_resids_noZBS, sorted_muErrs, freq_bounds_over_which_to_find_peak = freq_bounds_over_which_to_find_peak)
    
    plt.scatter(frequencies[0:300] / (2.0 * np.pi) * 1000.0, np.sqrt(2.0 * np.array(normalized_coefs_noZBS[0:300])) )
    xlabel = 'Scaled Periodogram Frequency, ' + r'$f_n / h_{100}$' + ' (Gyr' + r'$^{-1}$' + ')'
    ylabel = 'Approximate Fourier Amplitude, ' + r'$\sqrt{2Q_n}$' + ' (mags)'
    plt.xlim([0.0, frequencies[300] / (2.0 * np.pi) * 1000.0 * 1.01])
    plt.ylim([0.0, 0.025])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title('Partial Periodogram of Distance Modulus Residuals')
    plt.show() 

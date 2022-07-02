import numpy as np
import matplotlib.pyplot as plt
import scipy
import cantrips as can
import AstronomicalParameterArchive as apa
import CosmologicalParameterArchive as cpa
import calculateMuForwOftFromz as calcMuForw
import loadSN as lsn
import emcee
import corner
from matplotlib import ticker, cm
import random
import sys

def param_fit_funct_full(zs, H0, OmM, w, OmL, OmR, Om0):
    zs = list(zs)
    w_of_funct = lambda w_zs, w = w: w
    funct_mus = getMusForCosmology(zs, [H0, OmM, OmL, OmR, Om0], w_of_funct)
    #print ('zs = ' + str(zs))
    return funct_mus

def readInCovarianceFile(covariance_file):
    f = can.readInColumnsToList(covariance_file)[0]

    N = int(float(f[0]))
    covmat_unraveled = np.array(f[1:],dtype='float')
    covmat = covmat_unraveled.reshape((N,N)).tolist()
    return covmat

def generateCovarianceFileForSubsetOfSurveys(new_covariance_file, HF_surveys_to_include,
                                             baseline_covariance_file = 'stubbs/SNIsotropyProject/OriginalSNDataFiles/Pantheon+SH0ES_122221_1.cov',
                                             cepheid_file_name = 'calibratorset.txt', calibrator_surveys_to_include = 'all',
                                             dir_base = '', sn_data_type = 'pantheon_plus' , surveys_to_pair = { 'FOUNDATION':'PS1MD', 'SWIFTNEW':'SWIFT'}):
    orig_covariances = [float(elem) for elem in can.readInColumnsToList(dir_base + baseline_covariance_file, verbose = 0)[0]]
    orig_n_sn = int(orig_covariances[0])
    orig_covariances = orig_covariances[1:]

    print ('HF_surveys_to_include = ' + str(HF_surveys_to_include))
    HF_surveys_to_include = HF_surveys_to_include + [key for key in surveys_to_pair.keys() if surveys_to_pair[key] in HF_surveys_to_include ]
    print ('HF_surveys_to_include = ' + str(HF_surveys_to_include))

    cepheid_sn_ids = lsn.loadCepheidSN(cepheid_file_name)

    all_sn = lsn.loadSN(1, pull_extinctions = 0, verbose = 0, data_type = sn_data_type)
    orig_surveys_in_order = [sn['survey'] for sn in all_sn]
    all_surveys = np.unique(orig_surveys_in_order).tolist()
    if HF_surveys_to_include == 'all':
        HF_surveys_to_include = all_surveys
    if calibrator_surveys_to_include == 'all':
        calibrator_surveys_to_include = all_surveys
    all_baseline_sn = can.flattenListOfLists([[sn for sn in all_sn if sn['survey'] == survey] for survey in all_surveys])
    all_baseline_ids = [sn['SNID'] for sn in all_baseline_sn]
    all_baseline_surveys = [sn['survey'] for sn in all_baseline_sn]
    all_baseline_cephPaired = [True if id in cepheid_sn_ids else False for id in all_baseline_ids]
    indeces_to_keep = [i for i in range(len(all_sn)) if ((all_baseline_cephPaired[i] and all_baseline_surveys[i] in calibrator_surveys_to_include) or all_baseline_surveys[i] in HF_surveys_to_include)]

    print ('len(orig_surveys_in_order) = ' + str(len(orig_surveys_in_order)))
    #indeces_to_keep = [i for i in range(orig_n_sn) if orig_surveys_in_order[i] in surveys_to_include ]
    print ('len(indeces_to_keep) = ' + str(len(indeces_to_keep)) )
    truncated_covariances = [ orig_covariances[j] for j in range(len(orig_covariances)) if (j % orig_n_sn in indeces_to_keep and j // orig_n_sn in indeces_to_keep) ]
    truncated_covariances = orig_covariances.copy()
    n_carry_through_covariances = 0
    for j in range(len(orig_covariances)):
        if (j % orig_n_sn in indeces_to_keep and j // orig_n_sn in indeces_to_keep) :
            truncated_covariances[j] = orig_covariances[n_carry_through_covariances]
            n_carry_through_covariances = n_carry_through_covariances + 1
    truncated_covariances = truncated_covariances[0:n_carry_through_covariances]
    print ('n_carry_through_covariances = ' + str(n_carry_through_covariances))
    print ('len(truncated_covariances) = ' + str(len(truncated_covariances)))
    can.saveListsToColumns([len(indeces_to_keep)] + truncated_covariances, new_covariance_file, dir_base, )

    return 0

def generateToySNeCovarianceFile(n_extra_sn, new_covariance_file, baseline_covariance_file = 'stubbs/SNIsotropyProject/OriginalSNDataFiles/Pantheon+SH0ES_122221_1.cov', dir_base = ''):
    """
    Take in the Pantheon+ SNe covariance file, and add new lines for new
       artificial SNe.  These new lines are just 0's (assume no
       covariances with artificial SNe).  The 0's just need to be placed
       at the right spots to preserve the actual covariances between the
       real data points.
    """
    orig_covariance = [float(elem) for elem in can.readInColumnsToList(dir_base + baseline_covariance_file, verbose = 0)[0]]
    orig_n_sn = int(orig_covariance[0])
    orig_covariance = orig_covariance[1:]
    orig_lines = [orig_covariance[i * orig_n_sn:(i+1) * orig_n_sn] for i in range(orig_n_sn)]
    new_lines = [orig_line + [0.0 for i in range(n_extra_sn)] for orig_line in orig_lines]
    new_lines = new_lines + [[0.0 for i in range(n_extra_sn + orig_n_sn)] for j in range(n_extra_sn)]

    new_covariances = can.flattenListOfLists(new_lines)
    can.saveListsToColumns([int(orig_n_sn + n_extra_sn)] + new_covariances, dir_base + new_covariance_file, '')

    return 1



def generateToySNeData(randomize_orig_sn = 1,
                       n_sn_per_survey = 200, surveys = ['A', 'B'],
                       H0 = 70.0, OmegaM = 0.3, Omega0 = 1.0, OmegaR = 0.0001,  w_params = [-1],
                       sn_data_type = 'pantheon_plus', cepheid_file_name = 'calibratorset.txt',
                       HF_surveys_to_include = 'all', z_ranges = [[0.001, 0.1], [0.1, 1.0]], mu_err = 'median',
                       colors = ['r','b','g','orange','cyan','magenta','purple','darkblue','darkgreen', 'skyblue',],
                        art_data_file = 'stubbs/variableMuFits/ArtificialSurveys/ArtificialSNe_CSP.csv', dir_base = ''):
    art_data_file = dir_base  + art_data_file
    n_surveys = len(surveys)
    cepheid_sn_ids = lsn.loadCepheidSN(cepheid_file_name)
    orig_sn = lsn.loadSN(1, pull_extinctions = 0, verbose = 0, data_type = sn_data_type)
    loadedWOfFunction = lambda zs: w_params[0]
    OmegaLambda = Omega0 - OmegaM - OmegaR

    all_surveys = np.unique([sn['survey'] for sn in orig_sn]).tolist()
    if HF_surveys_to_include == 'all':
        HF_surveys_to_include = all_surveys
    all_baseline_sn = can.flattenListOfLists([[sn for sn in orig_sn if sn['survey'] == survey] for survey in all_surveys])
    all_baseline_ids = [sn['SNID'] for sn in all_baseline_sn]
    all_baseline_surveys = [sn['survey'] for sn in all_baseline_sn]
    all_baseline_cephPaired = [True if id in cepheid_sn_ids else False for id in all_baseline_ids]
    indeces_to_use = [i for i in range(len(orig_sn)) if (all_baseline_cephPaired[i] or all_baseline_surveys[i] in HF_surveys_to_include)]

    baseline_sn = [all_baseline_sn[i] for i in indeces_to_use]
    baseline_ids = [sn['SNID'] for sn in baseline_sn]
    baseline_cephPaired = [sn['SNID'] for sn in baseline_sn]
    baseline_surveys = [sn['survey'] for sn in baseline_sn]
    baseline_zs = [sn['z'] for sn in baseline_sn]
    baseline_mus = [sn['mu'] for sn in baseline_sn]
    baseline_muErrs = [sn['muErr'] for sn in baseline_sn]
    median_errs = [np.median([ baseline_muErrs[i] for i in range(len(baseline_muErrs)) if (baseline_zs[i] < z_range[1] and baseline_zs[i] > z_range[0]) ]) for z_range in z_ranges]
    mean_errs = [np.mean([ baseline_muErrs[i] for i in range(len(baseline_muErrs)) if (baseline_zs[i] < z_range[1] and baseline_zs[i] > z_range[0]) ])  for z_range in z_ranges]
    baseline_colors = [sn['color'] for sn in baseline_sn]
    baseline_zs, baseline_mus, baseline_muErrs, baseline_ids, baseline_surveys, baseline_colors, baseline_cephPaired = can.safeSortOneListByAnother(baseline_zs, [baseline_zs, baseline_mus, baseline_muErrs, baseline_ids, baseline_surveys, baseline_colors, baseline_cephPaired ])

    scalar_params = [H0, OmegaM, OmegaLambda, OmegaR , Omega0]
    #print ('baseline_zs = ' + str(baseline_zs))
    """
    calc_mus = getMusForCosmology(baseline_zs, scalar_params, loadedWOfFunction )
    if randomize_orig_sn:
        baseline_mus = np.random.normal(calc_mus, baseline_muErrs).tolist()
    wOfFunction = lambda zs: fitter_to_plot.w_of_funct(zs, *w_params)
    param_fit_funct = lambda zs, H0, OmM, w: param_fit_funct_full(zs, H0, OmM, w, OmegaLambda, OmegaR , Omega0)
    zs_to_fit = baseline_zs
    mus_to_fit = baseline_mus
    init_guess = [H0, OmegaM, -1]
    #print ('[len(zs_to_fit), len(mus_to_fit), len(canon_mus), len(start_mus)] = ' + str([len(zs_to_fit), len(mus_to_fit), len(canon_mus), len(start_mus)] ))
    best_fit_params = scipy.optimize.curve_fit(param_fit_funct, zs_to_fit, mus_to_fit, p0 = init_guess, sigma = baseline_muErrs)[0]
    param_fit_funct = lambda zs, H0, OmM, w: param_fit_funct_full(zs, H0, OmM, w, OmegaLambda, OmegaR , Omega0)
    best_fit_params = scipy.optimize.curve_fit(param_fit_funct, zs_to_fit, mus_to_fit, p0 = init_guess)[0]
    H0, OmegaM, best_fit_w = [best_fit_params[0], best_fit_params[1], best_fit_params[2]]

    best_fit_w = best_fit_params[2]
    loadedWOfFunction = lambda zs, best_fit_w = best_fit_w: best_fit_w

    #all_cephPaired = [np.any([(i % n_sn_per_survey == j) for j in range(int(n_sn_per_survey * 0.05))]) for i in range(len(all_ids)) ]]

    """

    astro_arch = apa.AstronomicalParameterArchive()
    cosmo_arch = cpa.CosmologicalParameterArchive()
    km_s_Mpc_in_Myr = astro_arch.getKmPerSToPcPerYr()
    zs = can.flattenListOfLists([10.0 ** (np.random.random(n_sn_per_survey ) * (np.log10(z_ranges[i][1]) - np.log10(z_ranges[i][0])) + np.log10(z_ranges[i][0])) for i in range(len(surveys))])
    sorted_zs = can.safeSortOneListByAnother(zs, [zs]) [0]
    adjustedWOfFunction = lambda ts, taus, zs: loadedWOfFunction(zs)

    #residual_mu_calculator = calcMuForw.ResidualMuCalculatorForArbitraryWofT(wOfFunction = adjustedWOfFunction, initial_zs = sorted_zs, astro_archive = astro_arch, cosmo_archive = cosmo_arch,H0 = H0 * km_s_Mpc_in_Myr, OmM0 = OmegaM, OmL0 = OmegaLambda, Om0 = Omega0, OmR0 = OmegaR)

    if len(sorted_zs) > 0:
        residual_mu_calculator = calcMuForw.ResidualMuCalculatorForArbitraryWofT(wOfFunction = adjustedWOfFunction, initial_zs = sorted_zs, astro_archive = astro_arch, cosmo_archive = cosmo_arch,H0 = H0 * km_s_Mpc_in_Myr, OmM0 = OmegaM, OmL0 = OmegaLambda, Om0 = Omega0, OmR0 = OmegaR)
        sorted_mus = residual_mu_calculator.getMus()
        sorted_muErrs = can.flattenListOfLists([[(median_errs[j] if mu_err == 'median' else mean_errs[j] if mu_err == 'mean' else mu_err)  for i in range(n_sn_per_survey)] for j in range(len(z_ranges))])
        sorted_mus = [np.random.normal(sorted_mus[i], sorted_muErrs[i]) for i in range(len(sorted_mus))]
    else:
        sorted_mus, sorted_muErrs = [[], []]
    all_surveys = can.flattenListOfLists([[survey for i in range(n_sn_per_survey)] for survey in surveys])
    unordered_indeces = list(range(len(sorted_zs)))
    random.shuffle(unordered_indeces )
    all_indeces, all_zs, all_mus, all_muErrs = can.safeSortOneListByAnother(unordered_indeces, [unordered_indeces, sorted_zs, sorted_mus, sorted_muErrs])
    all_surveys = baseline_surveys + all_surveys
    all_zs = baseline_zs + all_zs
    all_mus = baseline_mus + all_mus
    all_muErrs = baseline_muErrs + all_muErrs
    #plt.scatter(all_zs, all_mus, c = 'r', marker = '.')
    #plt.scatter(baseline_zs, baseline_mus, c = 'g', marker = '.')
    #forced_mus = calcMuForw.ResidualMuCalculatorForArbitraryWofT(wOfFunction = lambda ts, taus, zs: -1, initial_zs = sorted_zs, astro_archive = astro_arch, cosmo_archive = cosmo_arch,H0 = 70.0 * km_s_Mpc_in_Myr, OmM0 = 0.3, OmL0 = 0.7, Om0 = 1.0, OmR0 = 0.0).getMus()
    #plt.plot(sorted_zs, forced_mus)
    #plt.show()
    all_ids = baseline_ids + ['SN' + str(i) for i in all_indeces]
    all_cephPaired = [np.any([(i % n_sn_per_survey == j) for j in range(int(n_sn_per_survey * 0.05))]) for i in range(len(all_ids)) ]
    all_cephPaired = [True if id in cepheid_sn_ids else False for id in all_ids]
    all_colors = can.flattenListOfLists([[color for i in range(n_sn_per_survey)] for color in colors[0:len(surveys)]])
    all_colors = baseline_colors + all_colors
    header = 'SNID, Surveys, z, mu, muErr, Color, CephPaired?'
    print ('[len(all_ids), len(all_surveys), len(all_zs), len(all_mus), len(all_muErrs), len(all_colors), len(all_cephPaired)] = ' + str([len(all_ids), len(all_surveys), len(all_zs), len(all_mus), len(all_muErrs), len(all_colors), len(all_cephPaired)]))
    can.saveListsToColumns([all_ids, all_surveys, all_zs, all_mus, all_muErrs, all_colors, all_cephPaired], art_data_file, '', header = header, sep = ', ')
    return 1

def getMusForCosmology(zs_to_calc_mus, scalar_params, wOfFunction, astro_archive = None, cosmo_archive = None):
    """
    The residual mu calculator assumes that the given w is a function of t, tau, and z.
    So you need to give a function that takes in all three of those values, even if it
    only depends on one of them.
    In this code, we assume w depends only on z.
    """
    if astro_archive == None:
        astro_archive = apa.AstronomicalParameterArchive()
    if cosmo_archive == None:
        cosmo_archive = cpa.CosmologicalParameterArchive()
    km_s_Mpc_in_Myr = astro_archive.getKmPerSToPcPerYr()
    #if wOfFunction == None:
    #    wOfFunction = self.default_w_of_funct
    adjustedWOfFunction = lambda ts, taus, zs: wOfFunction(zs)
    H0, OmegaM, OmegaLambda, OmegaR , Omega0 = scalar_params
    residual_mu_calculator = calcMuForw.ResidualMuCalculatorForArbitraryWofT(wOfFunction = adjustedWOfFunction, initial_zs = zs_to_calc_mus, astro_archive = astro_archive, cosmo_archive = cosmo_archive,
                                                               H0 = H0 * km_s_Mpc_in_Myr, OmM0 = OmegaM, OmL0 = OmegaLambda, Om0 = Omega0, OmR0 = OmegaR)
    #print ('[H0,  OmegaM, OmegaLambda, Omega0, OmegaR] = ' + str([H0,  OmegaM, OmegaLambda, Omega0, OmegaR] ))
    #print ('[residual_mu_calculator.H0, residual_mu_calculator.OmM0, residual_mu_calculator.OmL0, residual_mu_calculator.OmR0, residual_mu_calculator.Om0] = ' + str([residual_mu_calculator.H0, residual_mu_calculator.OmM0, residual_mu_calculator.OmL0, residual_mu_calculator.OmR0, residual_mu_calculator.Om0]))
    calc_mus = residual_mu_calculator.getMus()
    return calc_mus

def makePlotOfPosteriorContours(mcmc_output_files_sets, params_in_sequence, xlims, ylims,
                                dir_base = '', x_scaling = 1,
                                xlabel = '$\Omega_M$', ylabel = r'$w$', col_elem_width = 2, row_elem_width = 2,
                                data_load_dir = 'stubbs/variableMuFits/mcmcResults/',
                                delimiter = ', ', n_ignore = 1, n_ignore_end = 0, labelsize = 14,
                                n_axarr_cols = 5, bins = 50, xticks = None, yticks = None, in_plot_loc = [0.5, 0.5],
                                fitter_levels = can.niceReverse([0.0] + (1.0 - np.exp(-(np.arange(1.0, 3.1, 1.0) ** 2.0 )/ 2.0)).tolist()),
                                show_fig = 0, save_fig = 0, save_fig_name = 'TestSaveMakePlotOfPosteriorContours.pdf',
                                plot_save_dir = 'stubbs/variableMuFits/plots/PosteriorWidths/',
                                text_loc = [0.05, 0.99], plot_text = None, in_plot_text_size = 8,
                                contour_area_ax = None, area_contour_num = 0, contour_area_x_vals = None, area_plot_color = 'k', area_plot_marker = 'x',
                                n_points_to_lin_fit = 0, in_plot_lin_text = ['', ''], n_burn_in = 20000):
    data_load_dir = dir_base + data_load_dir
    plot_save_dir = dir_base + plot_save_dir
    n_axarr_rows = ((len(mcmc_output_files_sets) - 1) // n_axarr_cols) + 1
    n_axarr_elems = n_axarr_rows * n_axarr_cols
    f, axarr = plt.subplots(n_axarr_rows, n_axarr_cols, squeeze = False, figsize = (col_elem_width * n_axarr_cols, row_elem_width * n_axarr_rows))
    areas_to_plot = [-1 for i in mcmc_output_files_sets ]
    for i in range(len(mcmc_output_files_sets)):
        col_index = i % n_axarr_cols
        row_index = i // n_axarr_cols
        ax = axarr[row_index, col_index]
        mcmc_output_files = mcmc_output_files_sets[i]
        mcmc_mesh_stack = [[] for j in range(len(mcmc_output_files))]
        for j in range(len(mcmc_output_files)):
            mcmc_output_file = mcmc_output_files[j]
            mcmc_data_set = can.readInColumnsToList(data_load_dir + mcmc_output_file, n_ignore = n_ignore, n_ignore_end = n_ignore_end, delimiter = delimiter)
            #print ('mcmc_data_set = ' + str(mcmc_data_set))
            x_data = [float(elem) - float(mcmc_data_set[params_in_sequence[0] - 1][-1]) for elem in mcmc_data_set[params_in_sequence[0] - 1][n_burn_in:-1]]
            x_data = np.array(x_data) - np.median(x_data)
            y_data = [float(elem) - float(mcmc_data_set[params_in_sequence[1] - 1][-1]) for elem in mcmc_data_set[params_in_sequence[1] - 1][n_burn_in:-1]]
            y_data = np.array(y_data) - np.median(y_data)
            x_mesh, y_mesh, mcmc_mesh = getContourMeshFromMCMCChainComponents([x_data, y_data], xlims, ylims, bins )
            mcmc_mesh_stack[j] = mcmc_mesh
        full_mcmc_mesh = np.mean(mcmc_mesh_stack, axis = 0)
        sorted_full_mcmc_mesh = np.flip(np.sort(full_mcmc_mesh.flatten()))
        full_levels = determineContourLevelsFromMCMCPosterior(sorted_full_mcmc_mesh, fitter_levels, samples_already_sorted = 1)
        contours = ax.contour(x_mesh, y_mesh, full_mcmc_mesh, levels =  full_levels)
        contour_vertices = [contours.allsegs[-1 - i] for i in range(1, len(full_levels)  )] #nth outer contour is -1 - n
        contour_areas = [np.sum([np.abs(can.measureAreaOfContour(vertices_set)) for vertices_set in contour]) for contour in contour_vertices]
        print ('contour_areas = ' + str(contour_areas))
        areas_to_plot[i] = contour_areas[area_contour_num]
        if i == 0:
            if xticks == None:
                xticks = ax.get_xticks()
            if yticks == None:
                yticks = ax.get_yticks()
        ax.set_yticks(yticks)
        if col_index != 0:
            ax.set_yticklabels([])
        ax.set_xticks(xticks)
        if row_index != n_axarr_rows - 1:
            ax.set_xticklabels([])
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)
        if plot_text != None:
            ax.text(*text_loc, plot_text[i], horizontalalignment='left',verticalalignment='top', transform = ax.transAxes, fontsize = in_plot_text_size)
    for i in range(len(mcmc_output_files_sets), n_axarr_elems):
        col_index = i % n_axarr_cols
        row_index = i // n_axarr_cols
        ax = axarr[row_index, col_index]
        if col_index != 0:
            ax.set_yticklabels([])
        ax.set_xticks(xticks)
        ax.set_yticks(yticks)
        if row_index != n_axarr_rows - 1:
            ax.set_xticklabels([])
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)
    plt.subplots_adjust (wspace=0.05, hspace=0.05, bottom=0.15, left = 0.1)
    f.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
    plt.xlabel(xlabel, fontsize = labelsize)
    plt.ylabel(ylabel, fontsize = labelsize, labelpad = 15)

    if contour_area_ax == None:
        f2, contour_area_ax = plt.subplots(1,1)
    if contour_area_x_vals == None:
        contour_area_x_vals = range(len(areas_to_plot))
    print ('contour_area_x_vals = ' + str(contour_area_x_vals))
    contour_area_x_vals_to_plot = [val * x_scaling for val in contour_area_x_vals]
    print ('contour_area_x_vals_to_plot = ' + str(contour_area_x_vals_to_plot))
    contour_area_ax.scatter(contour_area_x_vals_to_plot, np.array(areas_to_plot) / areas_to_plot[0] * 100 - 100, color = area_plot_color, marker = area_plot_marker )
    x_lims, y_lims = [contour_area_ax.get_xlim(), contour_area_ax.get_ylim()]
    print ('x_lims = ' + str(x_lims))
    if n_points_to_lin_fit > 1:
        x_vals_to_plot, y_vals_to_plot = [ [val * x_scaling for val in contour_area_x_vals], np.array(areas_to_plot) / areas_to_plot[0] * 100 - 100]
        print ('x_vals_to_plot = ' + str(x_vals_to_plot))
        lin_fit = np.polyfit(contour_area_x_vals[:n_points_to_lin_fit], y_vals_to_plot[:n_points_to_lin_fit], 1)
        try:

            plateau_fit = scipy.optimize.curve_fit(plateauFitFunct, np.array(contour_area_x_vals), np.array(y_vals_to_plot), p0 = (lin_fit[0], 1, 1))
        except RuntimeError:
            print ('PLATEAU FIT FAILED')
            plateau_fit = [ (lin_fit[0], 1, 1), [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0],  [0.0, 0.0, 1.0]] ]

        print ('plateau_fit = ' + str(plateau_fit))
        x_lims, y_lims = [contour_area_ax.get_xlim(), contour_area_ax.get_ylim()]
        #contour_area_ax.plot(x_lims, np.poly1d(lin_fit)(x_lims), c = 'r', linestyle = '--', alpha = 0.5)
        contour_area_ax.plot(np.linspace(0.0, x_lims[1], 51), plateauFitFunct(np.linspace(0.0, x_lims[1], 51) / x_scaling, *plateau_fit[0]), c = 'r', linestyle = '--', alpha = 0.5)
        print ('x_lims = ' + str(xlims))
        print('plateauFitFunct(np.linspace(0.0, x_lims[1], 51) / x_scaling, *plateau_fit[0]) = ' + str(plateauFitFunct(np.linspace(0.0, x_lims[1], 51) / x_scaling, *plateau_fit[0])))
        contour_area_ax.set_xlim(x_lims)
        contour_area_ax.set_ylim(y_lims)
        #ax.text(r'$\\frac{}{d \\sigma_{G,S}} = $' + str(can.round_to_n(lin_fit[0], 3)), color = 'blue')
        #contour_area_ax.text(in_plot_loc[0], in_plot_loc[1],  in_plot_lin_text[0] + str(can.round_to_n(lin_fit[0] / 1000 * 25, 3)) +  in_plot_lin_text[1], color = 'blue', transform = contour_area_ax.transAxes, fontsize = in_plot_text_size)
        contour_area_ax.text(in_plot_loc[0], in_plot_loc[1],  in_plot_lin_text[0] + str(can.round_to_n(plateau_fit[0][0] * plateau_fit[0][1], 2)) +  in_plot_lin_text[1], color = 'blue', transform = contour_area_ax.transAxes, fontsize = in_plot_text_size)
    return axarr, contour_area_ax

def determineContourLevelsFromMCMCPosterior(mcmc_samples, fraction_levels, samples_already_sorted = 0):
    if samples_already_sorted:
        sorted_samples = mcmc_samples[:]
    else:
        sorted_samples = np.sort(mcmc_samples)
    mcmc_samples = mcmc_samples.tolist()
    mcmc_sums = [ 0.0 for sample in mcmc_samples ]
    total_sample_points = np.sum(mcmc_samples)
    mcmc_sums[0] = mcmc_samples[0] / total_sample_points
    for i in range(1, len(mcmc_sums)):
        mcmc_sums[i] = mcmc_sums[i-1] + mcmc_samples[i] / total_sample_points
    indeces_of_levels = [ np.argmin(np.abs(np.array(mcmc_sums) - level) )  for level in fraction_levels ]
    scaled_levels = [ mcmc_samples[index] for index in indeces_of_levels ]
    scaled_levels = np.unique(scaled_levels)
    #print ('mcmc_samples = ' + str(mcmc_samples))
    #print ('fraction_levels  = ' + str(fraction_levels ))
    #print ('scaled_levels = ' + str(scaled_levels))
    return scaled_levels

def getContourMeshFromMCMCChainComponents(mcmc_samples_to_plot,  xlims, ylims, bins ):
    #print ('mcmc_samples_to_plot = ' + str(mcmc_samples_to_plot))
    x_bin_edges = np.linspace(xlims[0], xlims[1], bins+1).tolist()
    x_bin_centers = [(x_bin_edges[index] + x_bin_edges[index+1]) / 2.0  for index in range(0, len(x_bin_edges)-1)]
    x_bin_edges[0], x_bin_edges[-1] = [-np.inf, np.inf]
    y_bin_edges = np.linspace(ylims[0], ylims[1], bins+1).tolist()
    y_bin_centers = [(y_bin_edges[index] + y_bin_edges[index+1]) / 2.0  for index in range(0, len(y_bin_edges)-1)]
    y_bin_edges[0], y_bin_edges[-1] = [-np.inf, np.inf]
    x_mesh, y_mesh = np.meshgrid(x_bin_centers, y_bin_centers)
    mcmc_mesh = np.zeros((len(x_bin_centers), len(y_bin_centers)))
    #print ('[len(ref_mcmc_samples_to_plot[0]), len(ref_mcmc_samples_to_plot[1])] = ' + str([len(ref_mcmc_samples_to_plot[0]), len(ref_mcmc_samples_to_plot[1])]))
    for l in range(len(mcmc_samples_to_plot[0])):
        x_val, y_val = [mcmc_samples_to_plot[0][l], mcmc_samples_to_plot[1][l]]
        x_index = [index for index in range(len(x_bin_centers)) if mcmc_samples_to_plot[0][l] < x_bin_edges[index+1] and mcmc_samples_to_plot[0][l] >= x_bin_edges[index]]
        x_index = x_index[0]
        y_index = [index for index in range(len(y_bin_centers)) if mcmc_samples_to_plot[1][l] < y_bin_edges[index+1] and mcmc_samples_to_plot[1][l] >= y_bin_edges[index]]
        #print ('[[y_bin_edges[index], y_bin_edges[index + 1]] for index in range(len(y_bin_centers))] = ' + str([[y_bin_edges[index], y_bin_edges[index + 1]] for index in range(len(y_bin_centers))]))
        y_index = y_index[0]
        mcmc_mesh[y_index, x_index] = mcmc_mesh[y_index, x_index] + 1
    return [x_mesh, y_mesh, mcmc_mesh]

"""
def plateauFitFunct(xs, A, n):
    print ('[A, n] = ' + str([A, n]))
    fitted_vals = A * (xs ** n)
    return fitted_vals
"""

def plateauFitFunct(xs, A, p, alpha):
    #print ('[A, p] = ' + str([A, p]))
    fitted_vals = A * (p - 1 / (xs ** alpha + 1/p) )
    return fitted_vals

def makePlotOfPosteriorWidths(mcmc_output_files_by_x_val_dict, param_in_sequence,
                              dir_base = '', bins = 50, x_scaling = 1, show_errs = 0,
                              xlabel = '$\Delta \mu_S$ Gaussian prior width (mags)', ylabel = r'$H_0$ MCMC posterior $\sigma$ (km/s/Mpc)', title = '',
                              data_load_dir = 'stubbs/variableMuFits/mcmcResults/',
                              delimiter = ', ', n_ignore = 1, n_ignore_end = 0, ax = None, labelsize = 8, xlims = None, ylims = None,
                              zero_uncertainties_at_zero_prior = 1, plot_sigmas = 1, n_points_to_lin_fit = 0, in_plot_lin_text = ['', ''], in_plot_loc = [0.5, 0.5], in_plot_text_size = 8,
                              ref_param_vals_for_other_ax = None, plot_color = 'k', plot_marker = 'x', show_all_chain_posteriors = 0,
                              show_fig = 0, save_fig = 0, save_fig_name = 'TestSaveMakePlotOfPosteriorWidths.pdf',
                              plot_save_dir = 'stubbs/variableMuFits/plots/PosteriorWidths/',
                              plot_load_dir = 'stubbs/variableMuFits/plots/PosteriorWidths/', n_burn_in = 20000):
    #param in sequence is: is this the first, second, third param analyzed?
    data_load_dir = dir_base  + data_load_dir
    plot_save_dir = dir_base + plot_save_dir
    plot_load_dir = dir_base + plot_load_dir
    x_vals = list(mcmc_output_files_by_x_val_dict.keys())
    y_val_stds_dict = {}
    y_val_means_dict = {}
    print ('mcmc_output_files_by_x_val_dict = ' + str(mcmc_output_files_by_x_val_dict))
    for x_val in x_vals:
        mcmc_output_files = mcmc_output_files_by_x_val_dict[x_val]
        read_in_vals = [[] for i in range(len(mcmc_output_files))]
        for i in range(len(mcmc_output_files)):
            mcmc_output_file = mcmc_output_files[i]
            mcmc_data = can.readInColumnsToList(data_load_dir + mcmc_output_file, delimiter = ', ', n_ignore = n_ignore, convert_to_float = 1, n_ignore_end = 0, verbose = 0)
            loaded_data = mcmc_data[param_in_sequence][n_burn_in:]
            #print ('loaded_data = ' + str(loaded_data))
            loaded_data = [elem - loaded_data[-1] for elem in loaded_data[0:-1]]
            #print ('loaded_data = ' + str(loaded_data))
            data_std, data_mean = [np.std(loaded_data), np.mean(loaded_data)]
            read_in_vals[i] = [data_std, data_mean]
        y_val_stds_dict[x_val] = [read_in_val[0] for read_in_val in read_in_vals]
        y_val_means_dict[x_val] = [read_in_val[1] for read_in_val in read_in_vals]
        print ('y_val_stds_dict = ' + str(y_val_stds_dict))
        print ('y_val_means_dict = ' + str(y_val_means_dict))
    x_vals_to_plot = [val * x_scaling for val in x_vals]
    mean_param_stds = [np.mean(y_val_stds_dict[x_val]) for x_val in x_vals]
    mean_param_means = [np.mean(y_val_means_dict[x_val]) for x_val in x_vals]
    if zero_uncertainties_at_zero_prior:
        zerod_y_val_std = mean_param_stds[0]
    else:
        zerod_y_val_std = 0.0
    #Since these are statistical uncertainties, we do the subtraction in quadrature.
    print ('zerod_y_val_std = ' + str(zerod_y_val_std ))
    if plot_sigmas:
        y_vals_to_plot = [np.sqrt(y_val ** 2.0 - zerod_y_val_std ** 2.0) if y_val >= zerod_y_val_std else -np.sqrt(zerod_y_val_std ** 2.0 - y_val ** 2.0) for y_val in mean_param_stds]
        y_vals_to_plot_errs = [np.std(y_val_stds_dict[x_val])  for x_val in x_vals]
    else:
        print ('mean_param_means = ' + str(mean_param_means))
        y_vals_to_plot = [y_val - mean_param_means[0] for y_val in mean_param_means]
        y_vals_to_plot_errs = [np.std(y_val_means_dict[x_val])  for x_val in x_vals]
    print ('[x_vals_to_plot, y_vals_to_plot, y_vals_to_plot_errs] = ' + str([x_vals_to_plot, y_vals_to_plot, y_vals_to_plot_errs]))

    if show_all_chain_posteriors:
        [ax.scatter([x_vals_to_plot[i] for y_val in y_val_stds_dict[x_vals[i]]], np.array(y_val_stds_dict[x_vals[i]]) - zerod_y_val_std, c = 'r', marker = '.') for i in range(len(x_vals_to_plot))]
    scat = ax.scatter(x_vals_to_plot, y_vals_to_plot, c = plot_color, marker = plot_marker)
    if show_errs: ax.errorbar(x_vals_to_plot, y_vals_to_plot, yerr = y_vals_to_plot_errs, color = 'k', fmt = 'none')
    if n_points_to_lin_fit > 1:
        lin_fit = np.polyfit(x_vals[:n_points_to_lin_fit], y_vals_to_plot[:n_points_to_lin_fit], 1)
        try:
            plateau_fit = scipy.optimize.curve_fit(plateauFitFunct, np.array(x_vals), np.array(y_vals_to_plot), p0 = (lin_fit[0], 1, 1))
        except(RuntimeError):
            print ('PLATEAU FIT FAILED')
            plateau_fit = [ (lin_fit[0], 1, 1), [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0],  [0.0, 0.0, 1.0]] ]
        print ('plateau_fit = ' + str(plateau_fit))
        xlims, ylims = [ax.get_xlim(), ax.get_ylim()]
        #ax.plot(xlims, np.poly1d(lin_fit)(xlims), c = 'r', linestyle = '--', alpha = 0.5)
        print ('Text here is: ' + in_plot_lin_text[0])
        ax.plot(np.linspace(0.0, xlims[1], 51), plateauFitFunct(np.linspace(0.0, xlims[1], 51) / x_scaling, *plateau_fit[0]), c = 'r', linestyle = '--', alpha = 0.5)
        #ax.text(r'$\\frac{}{d \\sigma_{G,S}} = $' + str(can.round_to_n(lin_fit[0], 3)), color = 'blue')
        #ax.text(in_plot_loc[0], in_plot_loc[1], in_plot_lin_text[0] + str(can.round_to_n(lin_fit[0] / 1000 * 25, 3)) +  in_plot_lin_text[1], color = 'blue', transform = ax.transAxes, fontsize = in_plot_text_size)
        ax.text(in_plot_loc[0], in_plot_loc[1], in_plot_lin_text[0] + str(can.round_to_n(plateau_fit[0][0] * plateau_fit[0][1], 2)) +  in_plot_lin_text[1], color = 'blue', transform = ax.transAxes, fontsize = in_plot_text_size)
    if not(xlims == None):
        ax.set_xlim(xlims)
    if not(ylims == None):
        ax.set_ylim(ylims)
    ax.set_xlabel(xlabel, fontsize = labelsize)
    ax.set_ylabel(ylabel, fontsize  = labelsize)
    ax.set_title(title, fontsize = labelsize + 2 )

    if ref_param_vals_for_other_ax != None:
        twin_ax = ax.twinx()
        if ylims == None:
            ylims = ax.get_ylim()
        print ('!!! ylims = ' + str(ylims) + '!!!')
        meas_param_central_val = ref_param_vals_for_other_ax[0]
        ref_param_mean = ref_param_vals_for_other_ax[1]
        ref_param_sigma = ref_param_vals_for_other_ax[2]
        #paired_sigmas = [np.sqrt(average_err ** 2.0 + ref_param_sigma ** 2.0) for average_err in std_average_errs]
        #param_diffs = [val - ref_param_mean for val in average_vals]
        print ('ylims = ' + str(ylims))
        print ('[ref_param_mean, ref_param_sigma] = ' + str([ref_param_mean, ref_param_sigma] ))
        print ('np.mean([y_val_means_dict[x_val] for x_val in x_vals]) = ' + str(np.mean([y_val_means_dict[x_val] for x_val in x_vals])))
        #ylim_sigmas = [np.sqrt(ylims[0] ** 2.0 + ref_param_sigma ** 2.0), np.sqrt(ylims[1] ** 2.0 + ref_param_sigma ** 2.0)]
        #ylim_sigma_increases = [np.sqrt((ylims[0] + zerod_y_val_std) ** 2.0 + ref_param_sigma ** 2.0), np.sqrt((ylims[1] + zerod_y_val_std) ** 2.0 + ref_param_sigma ** 2.0)]
        #ylim_sigma_increases = [np.sqrt((ylims[0] ** 2.0 + zerod_y_val_std) ** 2.0 + ref_param_sigma ** 2.0), np.sqrt((ylims[1] + zerod_y_val_std) ** 2.0 + ref_param_sigma ** 2.0)]
        #ylim_Nsigmas = [abs(ref_param_mean - np.mean([y_val_means_dict[x_val] for x_val in x_vals])) / ylim_sigmas[i] for i in range(len(ylim_sigmas))]
        #ylim_Nsigmas = [abs(ref_param_mean - meas_param_central_val) / ylim_sigmas[i] for i in range(len(ylim_sigmas))]
        #ylim_sigma_percentage_increases = [(np.sqrt((y_vals_to_plot[0] + zerod_y_val_std) ** 2.0 + ref_param_sigma ** 2.0) / sig - 1) * 100 for sig in ylim_sigmas]
        #print ('ylim_sigmas = ' + str(ylim_sigmas))
        #print ('ylim_sigma_percentage_increases = ' + str(ylim_sigma_percentage_increases))
        #print ('ylim_Nsigmas = ' + str(ylim_Nsigmas))
        #twin_ax.set_ylim(ylim_Nsigmas)
        #twin_ax.set_ylim(ylim_sigma_percentage_increases)
        init_yticks = [elem for elem in np.linspace(np.min([val for val in y_vals_to_plot if not(np.isnan(val)) ]), np.max([val for val in y_vals_to_plot if not(np.isnan(val)) ]), 4)]

        #extra_sigmas = [np.sqrt(tick ** 2.0 - zerod_y_val_std ** 2.0) for tick in init_yticks]
        extra_sigmas = [tick for tick in init_yticks]
        HT_sigmas = [np.sqrt(sigma ** 2.0 + ref_param_sigma ** 2.0) if sigma >= 0.0 else np.sqrt(-(sigma ** 2.0) + ref_param_sigma ** 2.0) for sigma in extra_sigmas]
        HT_seps = [can.round_to_n(abs(ref_param_mean - meas_param_central_val) / sigma, 3) for sigma in HT_sigmas]
        corrected_yticks = [np.sqrt((abs(ref_param_mean - meas_param_central_val) / sep) ** 2.0 - ref_param_sigma ** 2.0 ) for sep in HT_seps]
        print ('init_yticks = ' + str(init_yticks))
        print ('zerod_y_val_std = ' +str(zerod_y_val_std))
        print ('extra_sigmas = ' + str(extra_sigmas))
        print ('HT_sigmas  = ' + str(HT_sigmas ))
        print ('HT_seps  = ' + str(HT_seps ))
        print ('corrected_yticks = ' + str(corrected_yticks))
        yticklabels = [str(sep) + r'$\sigma$' for sep in HT_seps]
        twin_ax.set_ylim(ylims)
        twin_ax.set_yticks(init_yticks)
        twin_ax.set_yticklabels(yticklabels)
        twin_ax.set_ylabel(ref_param_vals_for_other_ax[3], fontsize = labelsize)

    if save_fig:
        plt.savefig(plot_load_dir + save_fig_name)
    if show_fig:
        plt.show()

    return scat


    loaded_data_sets = [can.readInColumnsToList(mcmc_posterior_file, file_dir = data_load_dir, n_ignore = n_ignore, n_ignore_end = n_ignore_end, delimiter = delimiter) for mcmc_posterior_file in mcmc_posterior_files]

    x_vals_sets_to_plot = can.flattenListOfLists([[float(elem) for elem in loaded_data_set[0]] for loaded_data_set in loaded_data_sets])
    param_val_sets = can.flattenListOfLists([[ float(elem) for elem in loaded_data_set[(param_in_sequence - 1) * 3 + 1] ] for loaded_data_set in loaded_data_sets])
    param_neg_errs_sets = [[ float(elem) for elem in loaded_data_set[(param_in_sequence - 1) * 3 + 2] ] for loaded_data_set in loaded_data_sets]
    param_plus_errs_sets = [[ float(elem) for elem in loaded_data_set[(param_in_sequence - 1) * 3 + 3] ]  for loaded_data_set in loaded_data_sets]
    average_errs_sets = can.flattenListOfLists([(np.array(param_neg_errs_sets[i]) + np.array(param_plus_errs_sets[i])) / 2.0 for i in range(len(mcmc_posterior_files))])

    full_vals_dict = {x_val:[] for x_val in x_vals_sets_to_plot}
    mean_vals_dict = {x_val:[] for x_val in x_vals_sets_to_plot}
    std_vals_dict = {x_val:[] for x_val in x_vals_sets_to_plot}
    full_errs_dict = {x_val:[] for x_val in x_vals_sets_to_plot}
    mean_errs_dict = {x_val:[] for x_val in x_vals_sets_to_plot}
    std_errs_dict = {x_val:[] for x_val in x_vals_sets_to_plot}

    for i in range(len(x_vals_sets_to_plot)):
        x_val = x_vals_sets_to_plot[i]
        full_vals_dict[x_val] = full_vals_dict[x_val] + [param_val_sets[i]]
        full_errs_dict[x_val] = full_errs_dict[x_val] + [average_errs_sets[i]]

    for x_val in full_vals_dict.keys():
        mean_vals_dict[x_val] = np.mean(full_vals_dict[x_val])
        std_vals_dict[x_val] = np.std(full_vals_dict[x_val]) / np.sqrt(len(full_vals_dict[x_val]))
        mean_errs_dict[x_val] = np.mean(full_errs_dict[x_val])
        std_errs_dict[x_val] = np.std(full_errs_dict[x_val]) / np.sqrt(len(full_errs_dict[x_val]))

    if ax == None:
        f, ax = plt.subplots(1,1)

    x_vals_to_plot = full_vals_dict.keys()
    average_vals = [mean_vals_dict[key] for key in x_vals_to_plot]
    average_errs = [mean_errs_dict[key] for key in x_vals_to_plot]
    std_average_errs = [std_errs_dict[key] for key in x_vals_to_plot]

    ax.scatter(x_vals_to_plot, average_errs, c = plot_color, marker = plot_marker)
    if len(mcmc_posterior_files) > 1:
        ax.errorbar(x_vals_to_plot, average_errs, yerr = std_average_errs, color = 'k', fmt = 'none')
    ax.set_xlabel(xlabel, fontsize = labelsize)
    ax.set_ylabel(ylabel, fontsize  = labelsize)

    if ref_param_vals_for_other_ax != None:
        twin_ax = ax.twinx()
        ylims = ax.get_ylim()
        ref_param_mean = ref_param_vals_for_other_ax[0]
        ref_param_sigma = ref_param_vals_for_other_ax[1]
        #paired_sigmas = [np.sqrt(average_err ** 2.0 + ref_param_sigma ** 2.0) for average_err in std_average_errs]
        #param_diffs = [val - ref_param_mean for val in average_vals]
        ylim_sigmas = [np.sqrt(ylims[0] ** 2.0 + ref_param_sigma ** 2.0), np.sqrt(ylims[1] ** 2.0 + ref_param_sigma ** 2.0)]
        ylim_Nsigmas = [abs(ref_param_mean - np.mean(average_vals)) / ylim_sigmas[i] for i in range(len(ylim_sigmas))]
        #print ('n_sigma_seps = ' + str(n_sigma_seps))
        twin_ax.set_ylim(ylim_Nsigmas)
        twin_ax.set_ylabel(ref_param_vals_for_other_ax[2])

    if save_fig:
        plt.savefig(plot_load_dir + save_fig_name)
    if show_fig:
        plt.show()
    return 1


class CosmicFitter:
    #Do cosmic fits, allowing for survey-by-survey offsets in \mu and cosmic parameters.
    # We distinguish between supernovae that are paired with cepheids and supernovae
    # that are in the Hubble flow - these guys are affected differently by the offsets
    # in cosmic fit parameters.

    def getScaledHOfZ(self, z, OmegaM, OmegaR, Omegalambda, wOfFunctionParams):
        H = np.sqrt( OmegaM * (1.0 + z) ** 3.0 + OmegaR * (1.0 + z) ** 4.0 + OmegaLambda * (1.0 + z) ** ((1.0 + self.w_of_funct(z, wOfFunctionParams)) * 3.0) )
        return H

    def getMuFromSimpleWofZCosmology(self, redshifts, cosmic_and_w_params):
        speedol = self.astro_arch.getc()
        scalar_params, wOfFunctionParams, muOffsets = self.getAllFitParametersFromListOfParamsToVary(cosmic_and_w_params, ['H0', 'OmM', 'Om0', 'OmR'] + self.w_id_strs)
        H0, OmegaM, OmegaLambda, OmegaR , Omega0 = scalar_params
        scaling_in_Mpc = speedol / H0 #H0 is in km/s/Mpc  c is in km/s
        #dl_int_arg = lambda zInt: 1.0 / np.sqrt( OmegaM * (1.0 + zInt) ** 3.0 + OmegaR * (1.0 + zInt) ** 4.0 + OmegaLambda * (1.0 + zInt) ** ((1.0 + self.w_of_funct(zInt, wOfFunctionParams)) * 3.0) )
        dl_int_arg = lambda zInt: 1.0 / self.getScaledHofZ(zInt, OmegaM, OmegaR, Omegalambda, wOfFunctionParams)
        dLs = [scaling_in_Mpc * (1.0 + z) * scipy.integrate.quad(dl_int_arg, 0.0, z) [0] for z in redshifts]
        mus = 25.0 + 5.0 * np.log10(dLs)
        return mus

    def getVariableWCosmicFunction(self, w_of_funct_str = 'normal'):
        if w_of_funct_str.lower() == 'wa':
            default_w_params = [-1.0, 0.0]
            w_param_bounds = [[-2.0, 0.0], [-1.0, 1.0]]
            w_of_funct = lambda z, params: params[0] + params[1] * np.array(z) / (1.0 + np.array(z))
            w_id_strs = ['w0','wa']
            self.truth_vals['w0'] = -1.0
            self.truth_vals['wa'] = 0.0
        elif w_of_funct_str.lower() == 'w0':
            default_w_params = [-1.0]
            w_param_bounds = [[-2.0, 0.0]]
            w_of_funct = lambda z, params: params[0]
            w_id_strs = ['w0']
            self.truth_vals['w0'] = -1.0
        else:
            default_w_params = []
            w_param_bounds = []
            w_of_funct = lambda z, params: -1.0
            w_id_strs = []

        return w_of_funct, default_w_params, w_param_bounds, w_id_strs

    def updateUsedSN(self, z_lims = [-np.inf, np.inf], only_include_cepheid_surveys = 0, do_not_include_cepheid_surveys = 0, match_no_ceph_to_cephs = 0,
                        surveys_to_include = 'all'):
        """
        Update the set of supernovae that we use in the analysis by various cuts.
        """
        if surveys_to_include == 'all':
            surveys_to_include = self.all_surveys
        self.z_lims = z_lims
        include_indeces = [i for i in range(len(self.full_sorted_zs)) if self.full_sorted_zs[i] > z_lims[0] and self.full_sorted_zs[i] < z_lims[1]]
        include_indeces = [ i for i in include_indeces if self.full_sorted_surveys[i] in surveys_to_include ]
        if only_include_cepheid_surveys:
            include_indeces = [i for i in include_indeces if self.full_sorted_surveys[i] in self.full_cepheid_surveys]
        elif do_not_include_cepheid_surveys:
            include_indeces = [i for i in include_indeces if not(self.full_sorted_surveys[i] in self.full_cepheid_surveys)]
        elif match_no_ceph_to_cephs:
            #We need to select supernovae that are NOT from cepheid survays to try to match the redshift distribution of supernovae that ARE from cepheid surveys.
            print ('Not ready yet.')
        self.sorted_zs, self.sorted_sn, self.sorted_mus, self.sorted_muErrs, self.sorted_plotColors, self.sorted_surveys, self.sorted_ids = [[arr[index] for index in include_indeces] for arr in [self.full_sorted_zs, self.full_sorted_sn, self.full_sorted_mus, self.full_sorted_muErrs, self.full_sorted_plot_colors, self.full_sorted_surveys, self.full_sorted_ids]]
        self.sorted_cov_matrix = [[self.full_sorted_cov_matrix[include_col_index][include_row_index] for include_row_index in include_indeces] for include_col_index in include_indeces]
        self.sorted_inv_cov_matrix = np.linalg.inv(self.sorted_cov_matrix)
        self.sorted_mu_weights = 1.0 / np.array(self.sorted_muErrs) ** 2.0
        self.cepheid_indeces = np.unique( [ i for i in range(len(self.sorted_ids)) if self.sorted_ids[i] in self.cepheid_sne] )
        self.non_cepheid_indeces = [index for index in range(len(self.sorted_sn)) if not(index in  self.cepheid_indeces)]
        self.nonCeph_zs, self.nonCeph_sn, self.nonCeph_mus, self.nonCeph_muErrs, self.nonCeph_mu_weights, self.nonCeph_plotColors, self.nonCeph_surveys, self.nonCeph_ids = [[arr[index] for index in self.non_cepheid_indeces] for arr in [self.sorted_zs, self.sorted_sn, self.sorted_mus, self.sorted_muErrs, self.sorted_mu_weights, self.sorted_plotColors, self.sorted_surveys, self.sorted_ids]]
        self.nonCeph_cov_matrix = [[self.sorted_cov_matrix[col_index][row_index] for row_index in self.non_cepheid_indeces] for col_index in self.non_cepheid_indeces]
        self.nonCeph_inv_cov_matrix = np.linalg.inv(self.nonCeph_cov_matrix)
        self.cepheid_surveys = np.unique([self.sorted_surveys[index] for index in self.cepheid_indeces  ])
        unique_surveys = np.unique(self.sorted_surveys)
        mean_redshifts_by_survey = [np.mean([self.sorted_zs[i] for i in range(len(self.sorted_zs)) if self.sorted_surveys[i] == survey]) for survey in unique_surveys]
        self.surveys = can.safeSortOneListByAnother(mean_redshifts_by_survey, [unique_surveys])[0]
        self.n_surveys = len(self.surveys)
        self.canon_mus = getMusForCosmology(self.nonCeph_zs, [self.H0, self.OmegaM, self.OmegaLambda, self.OmegaR , self.Omega0], lambda zs: self.w_of_funct(zs, self.default_w_params), astro_archive = self.astro_arch, cosmo_archive = self.cosmo_arch) #calcMuForw.ResidualMuCalculatorForArbitraryWofT(wOfFunction = self.w_of_funct, initial_zs = self.sorted_zs).getMus()
        self.null_chi_square = np.sum((np.array(self.nonCeph_mus) - np.array(self.canon_mus)) ** 2.0 / np.array(self.nonCeph_muErrs) ** 2.0  )
        print ('self.null_chi_square = ' + str(self.null_chi_square ) )
        if len(self.params_to_overplot) > 0:
            print ('Running fitting for basic cosmological model...')
            print ('[self.init_guess_params_to_overplot, self.params_to_overplot] = ' + str([self.init_guess_params_to_overplot, self.params_to_overplot] ))
            self.basic_mcmc_sampler = self.doFullCosmicFit( self.init_guess_params_to_overplot, self.params_to_overplot, do_init_minim = 1, n_mcmc_steps = self.overplot_mcmc_steps, n_mcmc_chains = self.overplot_mcmc_chains, do_all_mu_offsets = 0, overwrite_mcmc = 0, show_mcmc = 0, save_mcmc = 0, overplot_basic_plot = 0, z_lims = self.z_lims)
        return 1


    def readInSNData(self, randomize_sn = 0, force_to_background_cosmology = 0, z_lims = [-np.inf, np.inf], surveys_to_pair = { 'FOUNDATION':'PS1MD', 'SWIFTNEW':'SWIFT'},
                     toy_sn_surveys_to_add = [] ):
        """
        We need to read in the supernovae that we use in the analysis.
        """

        print ('sn_data_type = ' + str(self.sn_data_type))
        if self.sn_data_type == 'toy_surveys':
            print ('Here 1')
            print ('self.sn_toy_data_file = ' + str(self.sn_toy_data_file))
            orig_ids, orig_surveys, orig_zs, orig_mus, orig_muErrs, orig_snPlotColors, cepheid_paired = can.readInColumnsToList(self.sn_toy_data_file, n_ignore = 1, delimiter = ', ')
            orig_surveys = [surveys_to_pair[survey] if survey in surveys_to_pair.keys() else survey for survey in orig_surveys ]

            orig_mus = [float(elem) for elem in orig_mus]
            orig_muErrs = [float(elem) for elem in orig_muErrs]
            orig_zs = [float(elem) for elem in orig_zs]
            cepheid_paired = [1 if elem == 'True' else 0 for elem in cepheid_paired]
            orig_sn = [{'SNID':orig_ids[i], 'survey':orig_surveys[i], 'z':orig_zs[i], 'mu':orig_mus[i], 'muErr':orig_muErrs[i], 'color':orig_snPlotColors[i]} for i in range(len(orig_ids))]
            cepheid_sn_ids = [orig_ids[i] for i in range(len(orig_ids)) if cepheid_paired[i]]
        else:
            print ('Here 2')
            orig_sn = lsn.loadSN(1, pull_extinctions = 0, verbose = 0, data_type = self.sn_data_type)
            orig_ids = [sn['SNID'] for sn in orig_sn]
            orig_surveys = [sn['survey'] for sn in orig_sn]
            orig_surveys = [surveys_to_pair[survey] if survey in surveys_to_pair.keys() else survey for survey in orig_surveys ]
            orig_zs = [sn['z'] for sn in orig_sn]
            orig_muErrs = [sn['muErr'] for sn in orig_sn]
            orig_mus = [sn['mu'] for sn in orig_sn]
            orig_snPlotColors = [sn['color'] for sn in orig_sn]
            cepheid_sn_ids = lsn.loadCepheidSN(self.cepheid_file_name)
        #print ('len(orig_sn) = ' + str(len(orig_sn)))
        #print ('[orig_ids[0], orig_surveys[0], orig_zs[0], orig_muErrs[0], orig_mus[0]] = ' + str([orig_ids[0], orig_surveys[0], orig_zs[0], orig_muErrs[0], orig_mus[0]]))
        unique_surveys = np.unique(orig_surveys)
        indeces_by_id = {}
        for i in range(len(orig_sn)):
            id = orig_ids[i]
            if id in indeces_by_id:
                indeces_by_id[id] = indeces_by_id[id]  + [i]
            else:
                indeces_by_id[id] = [i]
        """
        all_zs, all_mus, all_muErrs, all_snPlotColors, all_sn, all_ids, all_surveys = [[-1 for id in indeces_by_id.keys()] for arr in [orig_zs, orig_mus, orig_muErrs, orig_snPlotColors, orig_sn, orig_ids, orig_surveys]]
        print ('[np.max(orig_muErrs), all_ids[np.argmax(orig_muErrs)]] = ' + str([np.max(orig_muErrs), orig_ids[np.argmax(orig_muErrs)]]))
        for i in range(len(indeces_by_id.keys())):
            sn_id =  list(indeces_by_id.keys())[i]
            orig_id_indeces = indeces_by_id[sn_id]
            all_zs[i] = np.mean([orig_zs[index] for index in orig_id_indeces])
            all_mus[i] = np.average([orig_mus[index] for index in orig_id_indeces], weights = [orig_muErrs[index] ** -2.0 for index in orig_id_indeces])
            all_muErrs[i] = np.sqrt(np.sum([orig_muErrs[index] ** 2.0 for index in orig_id_indeces]) / np.sqrt(len(orig_id_indeces)))
            all_snPlotColors[i] = orig_snPlotColors[orig_id_indeces[0]]
            all_surveys[i] = orig_surveys[orig_id_indeces[0]]
            all_sn[i] = orig_sn[orig_id_indeces[0]]
            all_ids[i] = sn_id
        """
        all_zs, all_mus, all_muErrs, all_snPlotColors, all_sn, all_ids, all_surveys = [orig_zs, orig_mus, orig_muErrs, orig_snPlotColors, orig_sn, orig_ids, orig_surveys]
        print ('self.covariance_file = ' + str(self.covariance_file))
        off_diag_covariance = readInCovarianceFile(self.covariance_file)
        print ('len(all_muErrs) = ' + str(len(all_muErrs)))
        diag_covariance = [[0.0 for i in range(len(all_muErrs))] for j in range(len(all_muErrs))]
        for i in range(len(all_muErrs)): diag_covariance[i][i] = all_muErrs[i] ** 2.0
        full_covariance = (np.array(off_diag_covariance) + np.array(diag_covariance)).tolist()
        #full_covariance = np.array(diag_covariance)
        semi_sorted_cov_matrix =  can.safeSortOneListByAnother(all_zs, full_covariance)
        sorted_zs, sorted_sn, sorted_muErrs, sorted_mus, sorted_surveys, sorted_ids, sorted_plot_colors, sorted_cov_matrix = can.safeSortOneListByAnother(all_zs, [all_zs, all_sn, all_muErrs, all_mus, all_surveys, all_ids, all_snPlotColors, semi_sorted_cov_matrix])
        cepheid_indeces = np.unique( [ i for i in range(len(sorted_ids)) if sorted_ids[i] in cepheid_sn_ids ] )
        #print ('Following sn ids are repeated N times, and is (not) a cepheid sn_ext_errs: ')
        #for key in indeces_by_id:
        #    if len(indeces_by_id[key]) > 1: print(key + ' => ' + str(len(indeces_by_id[key])) + ' ' + str(key in cepheid_sn_ids))
        if randomize_sn:
            scalar_params = [self.H0, self.OmegaM, self.OmegaLambda, self.OmegaR, self.Omega0]
            wOfFunctionParams = self.default_w_params
            loadedWOfFunction = lambda zs: self.w_of_funct(zs, wOfFunctionParams[:])
            calc_mus = getMusForCosmology(sorted_zs, scalar_params, loadedWOfFunction, astro_archive = self.astro_arch, cosmo_archive = self.cosmo_arch, )
            if force_to_background_cosmology:
                sorted_mus = calc_mus[:]
            else:
                sorted_mus = np.random.normal(calc_mus, sorted_muErrs)
        #sorted_zs = sorted(all_zs)
        #s_to_yr = self.astro_arch.getSecondToYear()
        #if self.lookback_time_type.lower() in ['tau','taus']:
        #    all_lookBacks = [10 ** (-6.0) * tau * s_to_yr for tau in all_taus]
        #else:
        #    all_lookBacks = [10 ** (-6.0) * t * s_to_yr for t in all_ts]
        #sn_weights = [1.0 / err ** 2.0 for err in sorted_muErrs]
        #full_mean = sum([sn_weights[i] * sorted_muResids[i] for i in range(len(sn_weights))]) / (sum(sn_weights))
        #sorted_resids = [resid - full_mean for resid in all_muResids]
        return sorted_zs, sorted_sn, sorted_muErrs, sorted_mus, sorted_surveys, sorted_ids, sorted_plot_colors, sorted_cov_matrix, cepheid_sn_ids, cepheid_indeces
        #return sorted_sn, sorted_zs, sorted_lookBacks, sorted_mus, sorted_muErrs, sorted_resids, sorted_zerod_mus, sorted_surveys, sorted_ids, sorted_plot_colors, cepheid_sn_ids

    def getAllFitParametersFromListOfParamsToVary(self, new_params, params_to_vary):
        """
        Params to vary can be:
        ids of w params (typically w0 or w0 and wa) => parameters for w of z function
        'H0' => H0/H0_assumed (unitless)
        'L' => SNe luminosity / Assumed SN luminosity
        'OmM' => OmegaM
        'OmR' => OmegaR
        'Om0' => Omega0
        'mus' => All muOffsets
        OmegaLambda is calculated as: Om0 - OmR - OmM

        Note that those parameters that are really sets of parameters (w model parameters and
         mu offsets), this key string must appear everywhere that one of those parameters is given
        """
        # H0, OmegaM, OmegaR, Omega0, OmegaLambda, muOffsets, wOfFunctionParams,
        scalar_params = [-1 for param in self.scalar_param_strs]
        for i in range(len(self.scalar_param_strs)):
            param_str = self.scalar_param_strs[i]
            if param_str in params_to_vary:
                scalar_params[i] = new_params[params_to_vary.index(param_str)]
            else:
                scalar_params[i] = self.scalar_params[i]
        scalar_params[self.scalar_param_strs.index('OmL')] = scalar_params[self.scalar_param_strs.index('Om0')] - scalar_params[self.scalar_param_strs.index('OmM')] - scalar_params[self.scalar_param_strs.index('OmR')]
        wOfFunctionParams = [0.0 for param in self.default_w_params]
        for i in range(len(self.w_id_strs)):
            param_str = self.w_id_strs[i]
            if param_str in params_to_vary:
                wOfFunctionParams[i] = new_params[params_to_vary.index(param_str)]
            else:
                wOfFunctionParams[i] = self.default_w_params[i]
        mu_offsets_by_survey_strs = ['mu' + survey for survey in self.surveys]
        mu_offsets = [ 0.0 for offset_str in mu_offsets_by_survey_strs ]
        for i in range(len(mu_offsets_by_survey_strs)):
            mu_offset_str = mu_offsets_by_survey_strs[i]
            if mu_offset_str  in params_to_vary:
                mu_offsets[i] = new_params[params_to_vary.index(mu_offset_str)]
            else:
                mu_offsets[i] = self.default_mu_offset

        return([scalar_params, wOfFunctionParams, mu_offsets])

    def mcmc_log_prior(self, param_vals, params_to_vary):
        #includes all params, whether they were varied or not.  They always come out in a fixed order.

        scalar_params, wOfFunctionParams, muOffsets = self.getAllFitParametersFromListOfParamsToVary(param_vals, params_to_vary)
        H0, OmegaM, OmegaLambda, Omega0, OmegaR = scalar_params
        total_prior = 0.0
        for i in range(len(scalar_params)):
            scalar_param = scalar_params[i]
            cosmic_param_range = self.cosmic_param_bounds[i][1] - self.cosmic_param_bounds[i][0]
            if scalar_param > self.cosmic_param_bounds[i][1] or scalar_param < self.cosmic_param_bounds[i][0] :
                #print ('Out of Bound A: scalar_param = ' + str(scalar_param) + '; self.cosmic_param_bounds[i] = ' + str(self.cosmic_param_bounds[i]))
                #total_prior = total_prior + (-np.inf)
                dist_out_of_bounds = (scalar_param - self.cosmic_param_bounds[i][1] if scalar_param > self.cosmic_param_bounds[i][1] else self.cosmic_param_bounds[i][0] - scalar_param)
                total_prior = total_prior + (-np.exp(dist_out_of_bounds / cosmic_param_range * 100))
            else:
                total_prior = total_prior + np.log(1.0 / (cosmic_param_range))
        for i in range(len(wOfFunctionParams)):
            w_param = wOfFunctionParams[i]
            w_param_range = self.w_param_bounds[i][1] - self.w_param_bounds[i][0]
            if w_param > self.w_param_bounds[i][1] or scalar_params[i] < self.w_param_bounds[i][0] :
                #print ('Out of Bound B')
                #total_prior = total_prior + (-np.inf)
                dist_out_of_bounds = (scalar_param - self.w_param_bounds[i][1] if scalar_param > self.w_param_bounds[i][1] else self.w_param_bounds[i][0] - scalar_param)
                total_prior = total_prior + (-np.exp(dist_out_of_bounds / w_param_range * 100))
            else:
                total_prior = total_prior + np.log(1.0 / (w_param_range))
        for i in range(len(muOffsets)):
            mu_offset = muOffsets[i]
            survey = self.surveys[i]
            mu_prior_mean, mu_prior_sigma = self.muOffsetPriors[survey]
            if self.mu_prior_type in ['gaussian', 'gauss', 'normal', 'norm']:
                prior_addition = np.log( 1.0 / np.sqrt(2.0 * np.pi * mu_prior_sigma ** 2.0) * np.exp(- (mu_offset - mu_prior_mean)**2.0 / (2.0 * mu_prior_sigma ** 2.0)) )
            else:
                if mu_offset > (mu_prior_mean + mu_prior_sigma) or mu_offset < (mu_prior_mean - mu_prior_sigma):
                    #print ('Out of Bound C')
                    prior_addition = -np.inf
                else:
                    prior_addition = np.log(1.0 / (mu_prior_sigma * 2))
            total_prior = total_prior + prior_addition
            #print ('[survey, mu_offset, mu_prior_mean, mu_prior_sigma, prior_addition ] = ' + str([survey, mu_offset, mu_prior_mean, mu_prior_sigma, prior_addition ]))
        return total_prior

    #def determineLuminosityOffsetFromCepheids(self, mu_offset_uncertainties, ceph_surveys, mu_offsets_by_survey_dict, cepheid_paired_sn_indeces, extra_cepheid_offsets_by_survey = {} ):
    def determineLuminosityOffsetFromCepheids(self, ceph_inv_cov_matrix, ceph_surveys, mu_offsets_by_survey_dict, cepheid_paired_sn_indeces, extra_cepheid_offsets_by_survey = {} ):
        #cepheid_paired_sn_indeces = [i for i in range(len(ids)) if ids[i] in cepheid_sne]
        n_cepheid_paired_sne = len(cepheid_paired_sn_indeces)
        if n_cepheid_paired_sne == 0:
            return 0.0
        else:
            mu_offsets = [ mu_offsets_by_survey_dict[ceph_surveys[i]] + (extra_cepheid_offsets_by_survey[ceph_surveys[i]] if ceph_surveys[i] in extra_cepheid_offsets_by_survey.keys() else 0.0) for i in range(len(ceph_surveys))]
            #print ('luminosity_ratios for each of the Cepheid paired supernovae are: ' + str(luminosity_ratios))
            #mu_offset_uncertainties = [ mu_errs[i]  for i in cepheid_paired_sn_indeces ]
            #mean_luminosity_offset = np.average(mu_offsets, weights = np.array(mu_offset_uncertainties) ** -2.0 )
            mean_luminosity_offset = np.sum(np.dot(ceph_inv_cov_matrix, mu_offsets)) / np.sum(ceph_inv_cov_matrix)
            return mean_luminosity_offset

    def determineLuminosityRatioFromCepheids(self, mu_errs, ids, surveys, mu_offsets_by_survey_dict, cepheid_sne):
        """
        Get the fractional increase in inferred SNe absolute luminosity, if we offset the
            mus by some small amount.  This is actually not quite what we want, since
            we are interested in this ratio after we put it into a logarithm anyway.
        """
        cepheid_paired_sn_indeces = [i for i in range(len(ids)) if ids[i] in cepheid_sne]
        n_cepheid_paired_sne = len(cepheid_paired_sn_indeces)
        if n_cepheid_paired_sne == 0:
            return 1.0
        else:
            luminosity_ratios = [ 10.0 ** (- mu_offsets_by_survey_dict[surveys[i]] / 2.5 ) for i in cepheid_paired_sn_indeces ]
            #print ('luminosity_ratios for each of the Cepheid paired supernovae are: ' + str(luminosity_ratios))
            luminosity_ratio_uncertainties = [10.0 ** (-mu_offsets_by_survey_dict[surveys[i]] / 2.5) * np.log(10.0) / 2.5 * mu_errs[i] for i in cepheid_paired_sn_indeces ]
            mean_luminosity_ratio = np.average(luminosity_ratios, weights = np.array(luminosity_ratio_uncertainties) ** -2.0 )
            return mean_luminosity_ratio

    #def getPredictedMus(self, zs, ceph_muErrs, ceph_surveys, params, mu_offsets_by_survey_dict, cepheid_sne_indeces, WOfFunction, verbose = 0, extra_lum_offset = 0.0, extra_cepheid_offsets_by_survey = {}, forced_luminosity_offset = None ):
    def getPredictedMus(self, zs, ceph_inv_cov_matrix, ceph_surveys, params, mu_offsets_by_survey_dict, cepheid_sne_indeces, WOfFunction, verbose = 0, extra_lum_offset = 0.0, extra_cepheid_offsets_by_survey = {}, forced_luminosity_offset = None ):
        calc_mus = getMusForCosmology(zs, params, WOfFunction, astro_archive = self.astro_arch, cosmo_archive = self.cosmo_arch)
        #luminosity_ratio = self.determineLuminosityRatioFromCepheids(muErrs, ids, surveys, mu_offsets_by_survey_dict, cepheid_sne_list)
        """
        Note the "luminosity offset" is the derived mean magnitude of shifts in the Cepheid cepheid paired
            SNe.  So note that its sign shifts when we include it in the calculated mus.
        """
        luminosity_mu_offset = self.determineLuminosityOffsetFromCepheids(ceph_inv_cov_matrix, ceph_surveys, mu_offsets_by_survey_dict, cepheid_sne_indeces, extra_cepheid_offsets_by_survey = extra_cepheid_offsets_by_survey  )
        luminosity_mu_offset = luminosity_mu_offset + extra_lum_offset
        if forced_luminosity_offset != None:
            luminosity_mu_offset = forced_luminosity_offset
        if verbose: print ('luminosity_mu_offset = ' + str(luminosity_mu_offset))
        #calc_mus = (np.array(calc_mus) - 5.0 * np.log10(np.sqrt(luminosity_ratio))).tolist()
        #calc_mus = (np.array(calc_mus) + luminosity_mu_offset).tolist()
        calc_mus = (np.array(calc_mus) - luminosity_mu_offset).tolist()
        return calc_mus

    def readInSurveyMuOffsetCorrellationMatrix(self, correllation_file, delimiter = ','):
        """
        Reads in the correllation matrix in a reference file into
           the class memory.
        """
        cols = can.readInColumnsToList(correllation_file, delimiter = delimiter, n_ignore = 0)
        out_surveys = [survey.strip() for survey in cols[0][1:]]
        correllation_matrix = {out_surveys[i]:{col[0].strip():float(col[i+1]) for col in cols[1:]} for i in range(len(out_surveys))}
        return correllation_matrix

    def computeMuOffsetsFromCorrellationMatrix(self, uncorrellated_mu_offsets_by_survey_dict, correllation_matrix = None):
        """
        Multiply vector of mu offsets (numbers actually varied in
            MCMC) through the mu offsets correllation matrix by
            survey.  Also normalizes the offsets so that the
            total correllation in each matrix row is unity.
        No correllation means multiplying through by the unit
            matrix.
        """
        if correllation_matrix == None:
            correllation_matrix = self.survey_mu_correllation_matrix
        surveys_to_correllate = list(uncorrellated_mu_offsets_by_survey_dict.keys())
        correllated_mu_offsets_by_survey_dict = {}
        for survey in surveys_to_correllate:
            if survey in correllation_matrix:
                correllation_matrix_line = correllation_matrix[survey]
                normalization = np.sum([correllation_matrix_line[survey] if survey in correllation_matrix_line.keys() else 0.0 for survey in surveys_to_correllate])
                new_offset = np.sum([correllation_matrix_line[survey] * uncorrellated_mu_offsets_by_survey_dict[survey] if survey in correllation_matrix_line.keys() else 0.0 for survey in surveys_to_correllate])
                new_offset = new_offset / normalization
            else:
                new_offset = uncorrellated_mu_offsets_by_survey_dict[survey]
            correllated_mu_offsets_by_survey_dict[survey] = new_offset
        #print ('uncorrellated_mu_offsets_by_survey_dict = ' + str(uncorrellated_mu_offsets_by_survey_dict))
        #print ('correllated_mu_offsets_by_survey_dict = ' + str(correllated_mu_offsets_by_survey_dict))
        return correllated_mu_offsets_by_survey_dict

    def getMuOffsetsBySurvey(self, muOffsets, fix_mean_offset_survey = None):
        mu_offsets_by_survey_dict = { self.surveys[i]:muOffsets[i] for i in range(len(self.surveys)) }
        mu_offsets_by_survey_dict = self.computeMuOffsetsFromCorrellationMatrix(mu_offsets_by_survey_dict)
        if fix_mean_offset_survey != None:
            mu_offsets_by_mu = [mu_offsets_by_survey_dict[self.nonCeph_surveys[i]] if mu_offsets_by_survey_dict[self.nonCeph_surveys[i]] in mu_offsets_by_survey_dict[self.nonCeph_surveys[i]].keys() else 0.0 for i in range(len(self.nonCeph_mus))]
            weighted_mean_applied_offset = np.average(mu_offsets_by_mu, np.array(self.nonCeph_muErrs) ** -2.0 )
            corrective_offset = - np.sum([self.nonCeph_muErrs[i] ** -2.0 for i in range(len(self.nonCeph_muErrs)) if self.nonCeph_surveys[i] == fix_mean_offset_survey]) * np.sum(np.array(mu_offsets_by_mu) * np.array(self.nonCeph_muErrs) ** -2.0)
            mu_offsets_by_survey_dict[fix_mean_offset_survey] = corrective_offset
        return mu_offsets_by_survey_dict

    def FunctionToMinimizeAllParams(self, param_vals, params_to_vary, verbose = 0, return_prob = 0, subtract_off_mean = 0,
                                    fix_mean_offset_survey = None, extra_lum_offset = 0.0, fixed_survey_offsets = {}, forced_luminosity_offset = None  ):
         """
         Returns the log likelihood; useful for minimizing over in MCMCs (like emcee).
         """
         scalar_params, wOfFunctionParams, muOffsets = self.getAllFitParametersFromListOfParamsToVary(param_vals, params_to_vary)
         H0, OmegaM, OmegaLambda, OmegaR, Omega0 = scalar_params
         loadedWOfFunction = lambda zs: self.w_of_funct(zs, wOfFunctionParams[:])
         mu_offsets_by_survey_dict = self.getMuOffsetsBySurvey(muOffsets, fix_mean_offset_survey = fix_mean_offset_survey)
         if verbose: print ('mu_offsets_by_survey_dict = ' + str(mu_offsets_by_survey_dict))
         extra_offsets_by_survey = {survey:fixed_survey_offsets[survey][0] for survey in fixed_survey_offsets.keys() if fixed_survey_offsets[survey][1] in [0,1]}
         mu_offsets = np.array([mu_offsets_by_survey_dict[self.nonCeph_surveys[i]] + (extra_offsets_by_survey[self.nonCeph_surveys[i]] if self.nonCeph_surveys[i] in extra_offsets_by_survey.keys() else 0.0) for i in range(len(self.nonCeph_mus))])
         fixed_survey_offsets = fixed_survey_offsets
         extra_cepheid_offsets_by_survey = {survey:fixed_survey_offsets[survey][0] for survey in fixed_survey_offsets.keys() if fixed_survey_offsets[survey][1] in [0,2]}
         if verbose:
             print('fixed_survey_offsets = ' + str(fixed_survey_offsets))
             print ('extra_cepheid_offsets_by_survey = ' + str(extra_cepheid_offsets_by_survey))
         #calc_mus = self.getPredictedMus(self.nonCeph_zs, [self.sorted_muErrs[index] for index in self.cepheid_indeces], [self.sorted_surveys[index] for index in self.cepheid_indeces], scalar_params, mu_offsets_by_survey_dict, self.cepheid_indeces, loadedWOfFunction, verbose = verbose, extra_lum_offset = extra_lum_offset, extra_cepheid_offsets_by_survey = extra_cepheid_offsets_by_survey, forced_luminosity_offset = forced_luminosity_offset  )
         ceph_inv_cov_matrix = [ [self.sorted_inv_cov_matrix[col_index][row_index] for row_index in self.cepheid_indeces] for col_index in self.cepheid_indeces ]
         calc_mus = self.getPredictedMus(self.nonCeph_zs, ceph_inv_cov_matrix, [self.sorted_surveys[index] for index in self.cepheid_indeces], scalar_params, mu_offsets_by_survey_dict, self.cepheid_indeces, loadedWOfFunction, verbose = verbose, extra_lum_offset = extra_lum_offset, extra_cepheid_offsets_by_survey = extra_cepheid_offsets_by_survey, forced_luminosity_offset = forced_luminosity_offset  )
         if subtract_off_mean:
             mu_diffs = np.array(calc_mus) - np.array(self.nonCeph_mus)
             overall_wmean = can.weighted_mean(mu_diffs, self.nonCeph_muErrs)
             if verbose: print ('overall_wmean = ' + str(overall_wmean))
             overall_zerod_mus = np.array(self.nonCeph_mus) + overall_wmean
             #sorted_offset_mus = overall_zerod_mus - mu_offsets
             sorted_offset_mus = overall_zerod_mus - mu_offsets
         else:
             #sorted_offset_mus = np.array(self.sorted_mus) - mu_offsets
             sorted_offset_mus = np.array(self.nonCeph_mus) - mu_offsets
             #sorted_offset_mus = [ self.sorted_mus[i] - mu_offsets_by_survey_dict[self.sorted_surveys[i]] for i in range(len(self.sorted_mus)) ]
         #chi_square = np.sum( (np.array(calc_mus) - np.array(sorted_offset_mus)) ** 2.0 * np.array(self.nonCeph_mu_weights) )
         diffs = (np.array(sorted_offset_mus) - np.array(calc_mus))
         icov = self.nonCeph_inv_cov_matrix
         chi_square = np.dot(diffs, np.dot(icov, diffs))

         #log_likelihood = -0.5 * (chi_square + np.sum(np.log(2.0 * np.pi * np.array(self.nonCeph_mu_weights) )))
         log_likelihood = -0.5 * chi_square
         if np.isnan(log_likelihood): log_likelihood = -np.inf
         log_prior = self.mcmc_log_prior(param_vals, params_to_vary)
         dof = len(self.nonCeph_zs) - len(params_to_vary) - 1 #the -1 is the subtracted overall mean
         r_chi_square = chi_square / dof
         #chi_square = np.sum( (np.array(mu_diffs)) ** 2.0 * np.array(mu_weights) )

         #return the negative of the log likelihood, since we run minimization functions
         if return_prob:
             res_val = log_likelihood + log_prior
         else:
             res_val = r_chi_square
         if verbose:
             print ('[param_vals, params_to_vary] = ' + str([param_vals, params_to_vary]))
             print ('For: ' + str(self.scalar_param_strs) + ' | w(z) params | mu offsets = ')
             print (str(scalar_params) + ' | ' + str(wOfFunctionParams) + ' | ' + str( muOffsets))
             print ('[log_likelihood, log_prior] = ' + str([log_likelihood, log_prior]))


         #
         #

         #print ('param_vals = ' + str(param_vals.tolist()))
         #print ('[log_likelihood, log_prior, r_chi_square, res_val] = ' + str([log_likelihood, log_prior, r_chi_square, res_val] ))
         #f, axarr = plt.subplots(2)
         #axarr[0].scatter(self.nonCeph_zs, sorted_offset_mus , c = 'g', marker = '.')
         #axarr[0].scatter(self.nonCeph_zs, calc_mus, c = 'r', marker = '.')
         #axarr[1].scatter(self.nonCeph_zs, np.array(sorted_offset_mus ) - np.array(calc_mus), c = 'k', marker = '.')
         #plt.show()

         return res_val

    def doubleCheckParamsAreValid(self, param_vals, param_strs):
        valid_param_vals = []
        valid_param_strs = []
        for i in range(len(param_strs)):
            param_val = param_vals[i]
            param_str = param_strs[i]
            if param_str in self.scalar_param_strs:
                valid_param_strs = valid_param_strs + [param_str]
                valid_param_vals = valid_param_vals + [param_val]
            elif 'mu' in param_str:
                if param_str[2:] in self.surveys:
                    valid_param_strs = valid_param_strs + [param_str]
                    valid_param_vals = valid_param_vals + [param_val]
            elif param_str in self.w_id_strs:
                valid_param_strs = valid_param_strs + [param_str]
                valid_param_vals = valid_param_vals + [param_val]


        return [valid_param_vals, valid_param_strs]

    def calcQFromCosmology(self, redshift, scalar_cosmic_params, wOfFunctionParams):
        """
        Calculate the deceleration parameter, q(z), assuming that the dark energy EoS parameter
            is constant (though it doesn't needed to be -1).
        """
        H0, OmegaM, OmegaLambda, OmegaR, Omega0 = scalar_cosmic_params
        onePZ = 1.0 + redshift
        HOfZ = H0 * self.getScaledHOfZ(redshift, OmegaM, OmegaR, OmegaLambda, wOfFunctionParams)
        qDenominator = (2.0 * HOfZ ** 2.0)
        w = self.wOfFunction(redshift, wOfFunctionParams)
        qNumerator = onePZ * (3.0 * OmegaM * onePZ ** 2.0 + 4.0 * OmegaR * OnePZ ** 3.0 + 3 * (1 + w) * OmegaL * onePZ ** (3 * (1 + w) - 1))
        q = qNumerator / qDenominator - 1
        return q


    def calculateIntercept(self, redshift, varied_params, varied_param_vals, j0 = 1):
        scalar_params, wOfFunctionParams, muOffsets = self.getAllFitParametersFromListOfParamsToVary(varied_param_vals, varied_params)
        H0, OmegaM, OmegaLambda, OmegaR, Omega0 = scalar_params
        loadedWOfFunction = lambda zs: self.w_of_funct(zs, wOfFunctionParams[:])
        mu_offsets_by_survey_dict = self.getMuOffsetsBySurvey(muOffsets, fix_mean_offset_survey = fix_mean_offset_survey)
        q0 = self.calcQ0FromCosmology(self, 0.0, scalar_params, wOfFunctionParams)

        intercept = np.log10(1 + (1-q0) * redshift / 2.0 - (1 - q0 - 3 * q0 ** 2.0 + j0) * redshift ** 2.0 / 6.0) - 0.2 * mJ


    # ['w', 'w', 'w' 'H0', 'OmM', 'mu' ,'mu','mu', 'mu' ,'mu','mu', 'mu' ,'mu','mu', 'mu' ,'mu','mu', 'mu']
    def doFullCosmicFit( self, init_guess, params_to_vary,
                         do_init_minim = 0, n_mcmc_steps = 5000, n_mcmc_chains = 32, do_all_mu_offsets = 0, overwrite_mcmc = 1,
                         show_mcmc = 0, save_mcmc = 1, save_full_mcmc = 0,
                         fix_mean_offset_survey = None, additional_save_prefix = '', additional_save_suffix = '', verbose = 0,
                         overplot_basic_plot = 0, z_lims = None,
                         extra_luminosity_offset = 0.0, fixed_survey_offsets = {}, forced_luminosity_offset = None):
        """
        Params to vary can be:
        'w' => parameters for w of z function
        'H0' => H0
        'OmM' => OmegaM
        'OmR' => OmegaR
        'Om0' => Omega0
        'muOff' => muOffsets
        OmegaLambda is calculated as: Om0 - OmR - OmM
        'mu' + survey name => Survey offsets

        The fixed_luminosity_offset term is an amount by which we shift all SN that are cepheid paired in one direction.
            Useful as a diagnostic, as that should shift H0 in a predictable way, all else being equal.
        The fixed_survey_offsets specified fixed mu offsets, survey by survey.  The values of the dictionary should be a
            list, the first value of which is the offset and the second value is whether you should apply it to the
            Cepheid SNe of that survey only (0), the Hubble flow SNe of that survey (1), or all SNe from that survey (2).
            For example:
             fixed_survey_offsets = {'CFA3S':[0.05, 0]} applies a \Delta \mu of 0.05 to all SNe in CFA3S that have Cepheid pairs
             fixed_survey_offsets = {'CFA3S':[0.05, 1]} applies a \Delta \mu of 0.05 to all SNe in CFA3S do not have Cepheid pairs
             fixed_survey_offsets = {'CFA3S':[0.05, 2]} applies a \Delta \mu of 0.05 to all SNe in CFA3S
        """
        if z_lims != None and z_lims != self.z_lims:
            self.updateUsedSN(z_lims = z_lims)

        if do_all_mu_offsets:
            n_sn_by_survey = {survey:0 for survey in self.surveys}
            for survey in self.nonCeph_surveys:
                n_sn_by_survey[survey] = n_sn_by_survey[survey] + 1
            survey_with_max_visits = self.surveys[np.argmax([n_sn_by_survey[survey] for survey in self.surveys])]
            #params_to_vary = params_to_vary + ['mu' + survey for survey in self.surveys if survey != survey_with_max_visits]
            #init_guess = init_guess +  [self.default_mu_offset for survey in self.surveys if survey != survey_with_max_visits]
            params_to_vary = params_to_vary + ['mu' + survey for survey in self.surveys if survey != fix_mean_offset_survey]
            init_guess = init_guess +  [self.default_mu_offset for survey in self.surveys if survey != fix_mean_offset_survey]
        if len(params_to_vary) == 0:
            print ('Asked to run a minimization over no free parameters.  Returning empty list...')
            return []
        if do_init_minim:
            loaded_fit_funct = lambda params: self.FunctionToMinimizeAllParams(params, params_to_vary, verbose = 0, return_prob = 0, fix_mean_offset_survey = fix_mean_offset_survey,
                                                                               fixed_survey_offsets = fixed_survey_offsets, forced_luminosity_offset = forced_luminosity_offset)
            fit_bounds = []
            for param_to_vary in params_to_vary:
                if param_to_vary in self.scalar_param_strs:
                    new_param_bounds = self.cosmic_param_bounds[self.scalar_param_strs.index(param_to_vary)]
                elif param_to_vary in self.w_id_strs:
                    new_param_bounds = self.w_param_bounds[self.w_id_strs.index(param_to_vary)]
                elif param_to_vary[0:2] == 'mu':
                    new_param_bounds = self.mu_offset_bounds
                else:
                    print ('Param ' + str(param_to_vary) + ' not a known parameter/ ')
                    new_param_bounds = [-np.inf, np.inf]
                fit_bounds = fit_bounds + [new_param_bounds]
            #fit_bounds = ( [self.cosmic_param_bounds[self.scalar_param_strs.index(param_to_vary)] for param_to_vary in params_to_vary if param_to_vary in self.scalar_param_strs]
            #               + [self.w_param_bounds[self.w_id_strs.index(param_to_vary)] for param_to_vary in params_to_vary if param_to_vary in self.w_id_strs]
            #               # + (self.w_param_bounds if 'w' in params_to_vary else [])
            #               + [self.mu_offset_bounds for survey in self.surveys if 'mu' + survey in params_to_vary ] )
            #print ('[self.cosmic_param_bounds[self.scalar_param_strs.index(param_to_vary)] for param_to_vary in params_to_vary if param_to_vary in self.scalar_param_strs] = ' + str([self.cosmic_param_bounds[self.scalar_param_strs.index(param_to_vary)] for param_to_vary in params_to_vary if param_to_vary in self.scalar_param_strs]))
            #print ('self.cosmic_param_bounds = ' + str(self.cosmic_param_bounds))
            #print ('self.scalar_param_strs = ' + str(self.scalar_param_strs))
            if verbose: print ('[init_guess, fit_bounds, params_to_vary] = ' + str([init_guess, fit_bounds, params_to_vary] ))
            min_res = scipy.optimize.minimize(loaded_fit_funct, x0 = init_guess, bounds = fit_bounds)
            if verbose: print ('min_res = ' + str(min_res))
            init_guess = min_res.x
            print ('init_guess = ' + str(init_guess))
        init_guess, params_to_vary = self.doubleCheckParamsAreValid(init_guess, params_to_vary)
        loaded_fit_funct = lambda params: self.FunctionToMinimizeAllParams(params, params_to_vary, verbose = verbose, return_prob = 1, extra_lum_offset = extra_luminosity_offset, fixed_survey_offsets = fixed_survey_offsets, forced_luminosity_offset = forced_luminosity_offset )
        pos = init_guess + 1e-4 * np.random.randn(n_mcmc_chains, len(params_to_vary))
        nwalkers, ndim = pos.shape
        mcmc_sampler = emcee.EnsembleSampler(nwalkers, ndim, loaded_fit_funct, args=())
        #print ('self.nonCeph_surveys = ' + str(self.nonCeph_surveys ))
        mcmc_sampler.run_mcmc(pos, n_mcmc_steps, progress= True)
        #print ("min_res['x'].tolist() = " + str(min_res['x'].tolist() ))
        #fitted_scalar_params, fitted_wOfFunctionParams, fitted_muOffsets = self.getAllFitParametersFromListOfParamsToVary(min_res['x'], params_to_vary)
        #fitted_loadedWOfFunction = lambda zs: self.w_of_funct(zs, fitted_wOfFunctionParams[:])
        #calc_mus = self.getMusForCosmology(self.sorted_zs, fitted_scalar_params, wOfFunction = fitted_loadedWOfFunction)

        #print ('And for reference, the difference calculated weighted means for the these fit parameters are: ')
        #weighted_mean_diff_of_fitted_cosmology = self.calculateWeightedMeansBySurvey(self.sorted_mus, calc_mus, self.sorted_muErrs, self.sorted_surveys)
        #print (weighted_mean_diff_of_fitted_cosmology)
        full_fig_name = additional_save_prefix + 'MCMC_' + self.sn_data_type + ('_RAND' if self.randomize_sn else '_REAL' ) + '_'.join(params_to_vary) + '_NS' + str(n_mcmc_steps) + '_NC' + str(n_mcmc_chains) + additional_save_suffix + '.pdf'
        ref_fig_name = additional_save_prefix + 'MCMC_' + self.sn_data_type + ('_RAND' if self.randomize_sn else '_REAL' ) + '_OverPlot' + '_'.join(self.params_to_overplot) + '_NS' + str(n_mcmc_steps) + '_NC' + str(n_mcmc_chains) + additional_save_suffix + '.pdf'
        single_chain_title = additional_save_prefix + 'MCMC_' + self.sn_data_type + ('_RAND' if self.randomize_sn else '_REAL' ) + '_SingleChain' + '_'.join(params_to_vary) + '_NS' + str(n_mcmc_steps) + '_NC' + str(n_mcmc_chains) + additional_save_suffix + '.pdf'
        save_mcmc_file_name = additional_save_prefix + 'PosteriorStats_' + self.sn_data_type + ('_RAND' if self.randomize_sn else '_REAL' ) + '_'.join(params_to_vary) + '_NS' + str(n_mcmc_steps) + '_NC' + str(n_mcmc_chains) + additional_save_suffix + '.txt'
        save_full_mcmc_file_name = additional_save_prefix + 'MCMC_' + self.sn_data_type + ('_RAND' if self.randomize_sn else '_REAL' ) + '_'.join(params_to_vary) + '_NS' + str(n_mcmc_steps) + '_NC' + str(n_mcmc_chains) + additional_save_suffix + '.txt'
        if overwrite_mcmc:
            print ('overwriting mcmc... ')
            self.mcmc_sampler = mcmc_sampler
        if show_mcmc or save_mcmc:
            print ('I SHOULD BE SAVING!!!!')
            self.showMCMCResults(params_to_vary, mcmc_sampler = mcmc_sampler, save_fig = save_mcmc, show_fig = 0, full_fig_name = full_fig_name, single_chain_title = single_chain_title, ref_fig_name = ref_fig_name, overplot_basic_plot = overplot_basic_plot, fix_mean_offset_survey = fix_mean_offset_survey,  forced_luminosity_offset = forced_luminosity_offset, verbose = 0 )
        if save_mcmc:
            self.saveMCMCPosteriorMatrix(params_to_vary, save_mcmc_file_name, mcmc_sampler = mcmc_sampler, save_overplot_fit = overplot_basic_plot,)
        if save_full_mcmc:
            self.saveFullMCMC(params_to_vary, save_full_mcmc_file_name, mcmc_sampler = mcmc_sampler)
        return mcmc_sampler

    def saveFullMCMC(self, params_to_vary, save_file,
                     target_dir = None, mcmc_sampler = None, sep = ', ',
                     burn_in_steps = 1000, thinning = 5, save_rounding = 5,
                     rand_additions_by_col = {0:100, 1:1, 2:1}):
        """
        Save the full, flattened MCMC chains with the trimming and burn in applied
        """
        if target_dir == None:
            target_dir = self.root_dir + self.mcmc_results_dir
        if mcmc_sampler == None:
            flat_samples = self.mcmc_sampler.get_chain(discard=burn_in_steps, thin=thinning, flat=True)
        else:
            flat_samples = mcmc_sampler.get_chain(discard=burn_in_steps, thin=thinning, flat=True)
        header = sep.join(params_to_vary)
        flat_samples = np.transpose(flat_samples)
        additions_by_col = [(np.random.random() - 0.5) * rand_additions_by_col[i] if i in rand_additions_by_col.keys() else 0.0 for i in range(len(params_to_vary))]
        print('additions_by_col = ' + str(additions_by_col))
        flat_samples = [[can.round_to_n(elem * (1000 if params_to_vary[i][0:2] == 'mu' else 1) + additions_by_col[i], save_rounding) for elem in flat_samples[i]] + [additions_by_col[i]] for i in range(len(params_to_vary))]
        can.saveListsToColumns(flat_samples, save_file, target_dir, header = header, sep = sep)
        return 1

    def saveMCMCPosteriorMatrix(self, params_to_vary, save_mcmc_file_name,
                             target_dir = None, mcmc_sampler = None, burn_in_steps = 1000, thinning = 5,
                              sep = ', ', n_sig_figs = 5, save_overplot_fit = 0):
        """
        Save the covariance matrix of the given MCMC sampler
            for each of the varied parameters to a text file.
        """
        if target_dir == None:
            target_dir = self.root_dir + self.mcmc_results_dir
        if mcmc_sampler == None:
            flat_samples = self.mcmc_sampler.get_chain(discard=burn_in_steps, thin=thinning, flat=True)
        else:
            flat_samples = mcmc_sampler.get_chain(discard=burn_in_steps, thin=thinning, flat=True)
        flat_samples_mean = [np.mean(col) for col in flat_samples]
        header = sep.join(params_to_vary)
        ndim = len(params_to_vary)
        #best_fit_params = [ np.percentile(flat_samples[:, i], percentiles_to_print) for i in range(ndim) ]
        #best_fit_params = [ [can.round_to_n(param, n_sig_figs) for param in param_set] for param_set in best_fit_params ]
        flat_samples_t = np.transpose(flat_samples)
        matrix = np.cov(flat_samples_t)
        #for i in range(ndim):
        #    for j in range(ndim):
        #        col_i = [flat_sample[i] for flat_sample in flat_samples]
        #        col_j = [flat_sample[j] for flat_sample in flat_samples]
        #        print ('col_i = ' + str(col_i))
        #        print ('col_j = ' + str(col_j))
        #        rows_to_correlate = [[col_i[k], col_j[k]] for k in range(len(col_i))]
        #        print ('rows_to_correlate = ' + str(rows_to_correlate))
        #        corref_ij = np.corrcoef( rows_to_correlate )
        #        print ('[i,j,corref_ij] = ' + str([i,j,corref_ij]))

        #        matrix[i][j] = corref_ij
        print ('matrix = ' + str(matrix))
        save_rows = [ [can.round_to_n(row_elem, n_sig_figs) for row_elem in row] for row in matrix]
        print ('Saving covariance matrix to: ' + str(target_dir + save_mcmc_file_name))
        can.saveListToFile(save_rows, save_mcmc_file_name, save_dir = target_dir, sep = sep, append = False, header = header)

        return 1

    def showDistanceLadderFromFit(self, varied_params, param_vals, fix_mean_offset_survey = None, forced_luminosity_offset = None ):
        f, axarr = plt.subplots(2,2)
        scalar_params, w_params, mu_offsets = self.getAllFitParametersFromListOfParamsToVary(best_fit_params, varied_params)
        no_offset_mu_offsets_by_survey_dict = self.getMuOffsetsBySurvey([self.default_mu_offset for offset in mu_offsets], fix_mean_offset_survey = fix_mean_offset_survey)
        mu_offsets_by_survey_dict = self.getMuOffsetsBySurvey(mu_offsets, fix_mean_offset_survey = fix_mean_offset_survey)
        ceph_inv_cov_matrix = [ [self.sorted_inv_cov_matrix[col_index][row_index] for row_index in self.cepheid_indeces] for col_index in self.cepheid_indeces ]
        #calc_mus = self.getPredictedMus(self.nonCeph_zs, [self.sorted_muErrs[index] for index in self.cepheid_indeces], [self.sorted_surveys[index] for index in self.cepheid_indeces], best_scalar_params, no_offset_mu_offsets_by_survey_dict, self.cepheid_indeces, lambda zs: self.w_of_funct(zs, w_params), verbose = 1, forced_luminosity_offset = forced_luminosity_offset )
        calc_mus = self.getPredictedMus(self.nonCeph_zs, ceph_inv_cov_matrix, [self.sorted_surveys[index] for index in self.cepheid_indeces], best_scalar_params, no_offset_mu_offsets_by_survey_dict, self.cepheid_indeces, lambda zs: self.w_of_funct(zs, w_params), verbose = 1, forced_luminosity_offset = forced_luminosity_offset )
        shifted_meas_mus = np.array(self.sorted_mus) + np.array([mu_offsets_by_survey_dict[survey] for survey in self.sorted_surveys])
        axarr[0,1].scatter(shifted_meas_mus, calc_mus, color = 'k')
        axarr[0,1].set_xlabel(r'SN measured $\mu + \Delta \mu$')
        axarr[0,1].set_ylabel(r'SN inferred $\mu(z, H_0, \Omega_M)$')
        plt.show()

    def computeBestFitResultsFromMCMC(self, varied_params,
                                      percentiles_to_print = [(scipy.special.erf(n_sig / np.sqrt(2.0)) + 1 ) / 2 * 100 for n_sig in [-1,0,1]],
                                      n_sig_figs = 5,  mcmc_sampler = None,  fix_mean_offset_survey = None, burn_in_steps = 1000, thinning = 5,
                                      extra_lum_offset = 0.0, extra_cepheid_offsets_by_survey = {}, forced_luminosity_offset = None, verbose_weighted_means = 1
                                      ):
        if mcmc_sampler == None:
            flat_samples = self.mcmc_sampler.get_chain(discard=burn_in_steps, thin=thinning, flat=True)
            samples = self.mcmc_sampler.get_chain(discard=burn_in_steps, thin=thinning, flat=False)
        else:
            flat_samples = mcmc_sampler.get_chain(discard=burn_in_steps, thin=thinning, flat=True)
            samples = mcmc_sampler.get_chain(discard=burn_in_steps, thin=thinning, flat=False)
        ndim = len(varied_params)
        best_fit_params_and_errs = [np.percentile(flat_samples[:, i], [100.0 * np.e ** -2.0, 50.0, 100.0 - 100.0 * np.e ** -2.0]) for i in range(ndim) ]
        #print ('best_fit_params_and_errs = ' + str(best_fit_params_and_errs))
        best_fit_params = [param_set[1] for param_set in best_fit_params_and_errs]
        best_fit_param_errs = [[param_set[1] - param_set[0], param_set[2] - param_set[1]] for param_set in best_fit_params_and_errs]
        #print ('best_fit_params_and_errs = ' + str(best_fit_params_and_errs))
        best_scalar_params, best_w_params, best_offsets = self.getAllFitParametersFromListOfParamsToVary(best_fit_params, varied_params)
        #print ('[best_scalar_params, best_w_params, best_offsets] = ' + str([best_scalar_params, best_w_params, best_offsets]))
        #best_calc_mus = self.getMusForCosmology(self.sorted_zs, self.sorted_ids, self.sorted_mus, best_scalar_params, wOfFunction = lambda zs: self.w_of_funct(zs, best_w_params))
        no_shift_mu_offsets_by_survey_dict = self.getMuOffsetsBySurvey([self.default_mu_offset for offset in best_offsets], fix_mean_offset_survey = fix_mean_offset_survey)
        ceph_inv_cov_matrix = [ [self.sorted_inv_cov_matrix[col_index][row_index] for row_index in self.cepheid_indeces] for col_index in self.cepheid_indeces ]
        #best_calc_mus = self.getPredictedMus(self.nonCeph_zs, [self.sorted_muErrs[index] for index in self.cepheid_indeces], [self.sorted_surveys[index] for index in self.cepheid_indeces], best_scalar_params, no_shift_mu_offsets_by_survey_dict, self.cepheid_indeces, lambda zs: self.w_of_funct(zs, best_w_params), verbose = 0, extra_lum_offset = extra_lum_offset, extra_cepheid_offsets_by_survey = extra_cepheid_offsets_by_survey, forced_luminosity_offset = forced_luminosity_offset )
        best_calc_mus = self.getPredictedMus(self.nonCeph_zs, ceph_inv_cov_matrix, [self.sorted_surveys[index] for index in self.cepheid_indeces], best_scalar_params, no_shift_mu_offsets_by_survey_dict, self.cepheid_indeces, lambda zs: self.w_of_funct(zs, best_w_params), verbose = 0, extra_lum_offset = extra_lum_offset, extra_cepheid_offsets_by_survey = extra_cepheid_offsets_by_survey, forced_luminosity_offset = forced_luminosity_offset )
        weighted_mean_offsets = []
        weighted_mean_params = []
        truth_vals = []

        #These should be the best fit distance modulus offsets by survey - something that can be calculated analytically.
        #   These are used to spot check the MCMC posteriors.
        varied_surveys = [varied_param[2:] for varied_param in varied_params if varied_param[0:2] == 'mu']
        #weighted_mean_offsets_by_survey = self.calculateWeightedMeansBySurvey(self.sorted_mus, best_calc_mus, self.sorted_muErrs, varied_surveys)
        weighted_mean_offsets_by_survey = self.calculateWeightedMeansBySurvey(self.nonCeph_mus, best_calc_mus, self.nonCeph_muErrs, varied_surveys, [self.sorted_muErrs[i] for i in self.cepheid_indeces], [self.sorted_surveys[i] for i in self.cepheid_indeces], verbose = verbose_weighted_means )

        for i in range(ndim):
            varied_param = varied_params[i]
            mcmc = np.percentile( flat_samples[:, i],  percentiles_to_print )
            q = np.diff(mcmc)
            #print ('For param ' + varied_param + ': ' + str(mcmc[1]) + ' [ -' + str(q[0]) + ', +' + str(q[1]) + ']')
            if varied_param[0:2] == 'mu':
                survey = varied_param[2:]
                #mu_indeces_for_survey = [ j for j in range(len(self.sorted_surveys)) if self.sorted_surveys[j] == survey ]
                #calc_mus_for_survey, meas_mus_for_survey, mu_errs_for_survey = [np.array([best_calc_mus[index] for index in mu_indeces_for_survey]), np.array([overall_zerod_mus[index] for index in mu_indeces_for_survey]), np.array([self.sorted_muErrs[index] for index in mu_indeces_for_survey]) ]
                #calc_mus_for_survey, meas_mus_for_survey, mu_errs_for_survey = [np.array([best_calc_mus[index] for index in mu_indeces_for_survey]), np.array([self.sorted_mus[index] for index in mu_indeces_for_survey]), np.array([self.sorted_muErrs[index] for index in mu_indeces_for_survey]) ]
                #weighted_mean_mu_diff = -can.weighted_mean(calc_mus_for_survey - meas_mus_for_survey, mu_errs_for_survey )
                weighted_mean_mu_diff = weighted_mean_offsets_by_survey[survey]

                #print ('For reference, weighted mean offset for survey ' + varied_param[2:] + ' is: ' + str(weighted_mean_mu_diff ))
                #print ('The calculated mu offset for survey ' + varied_param[2:] + ' is: ' + str(mcmc[1]) )
                weighted_mean_offsets = weighted_mean_offsets + [weighted_mean_mu_diff]
                weighted_mean_params = weighted_mean_params + [ weighted_mean_mu_diff ]
                truth_vals = truth_vals + [weighted_mean_mu_diff]
            else:
                weighted_mean_params = weighted_mean_params + [best_fit_params[i] ]
                truth_vals = truth_vals + [self.truth_vals[varied_param]]
        best_fit_prob = self.FunctionToMinimizeAllParams(best_fit_params, varied_params, verbose = 1, return_prob = 1, forced_luminosity_offset = forced_luminosity_offset)
        weighted_mean_diff_prob = self.FunctionToMinimizeAllParams(weighted_mean_params, varied_params, verbose = 1, return_prob = 1, forced_luminosity_offset = forced_luminosity_offset)

        return [best_fit_params, best_fit_param_errs, best_fit_prob, weighted_mean_params, weighted_mean_diff_prob]

    def showMCMCResults(self, varied_params,
                        mcmc_sampler = None, target_dir = None, overplot_basic_plot = 0, bins = 50, n_single_chain_cols = 4, n_sig_figs = 5,
                        burn_in_steps = 1000, thinning = 5,
                         percentiles_to_print = [(scipy.special.erf(n_sig / np.sqrt(2.0)) + 1 ) / 2 * 100 for n_sig in [-1,0,1]],
                         save_fig = 0, save_single_chain = 0, show_fig = 1, contourf_cmap = 'plasma', contour_c = 'k', labelsize = 24, titlesize = 22, cols_to_center = [0,1,2],
                        full_fig_name = 'MCMC_fitter_results.pdf', ref_fig_name = 'MCMC_full_fitter_results_vs_ref.pdf', single_chain_title = 'MCMC_single_chain_results.pdf',
                        fitter_levels = can.niceReverse([0.0] + (1.0 - np.exp(-(np.arange(1.0, 3.1, 1.0) ** 2.0 )/ 2.0)).tolist()), background_fitter_levels = can.niceReverse([0.0] + (1.0 - np.exp(-(np.arange(1.0, 3.1, 1.0) ** 2.0 )/ 2.0)).tolist()),
                        fix_mean_offset_survey = None, ref_plot_alpha = 0.5, forced_plot_lims = None, overplot_elem_size = 5, forced_luminosity_offset = None, verbose = 1 ):
        percentiles_to_print = [can.round_to_n(percentile, n_sig_figs) for percentile in percentiles_to_print]
        if target_dir == None:
            target_dir = self.root_dir + self.plot_dir
        if forced_plot_lims == None:
            forced_plot_lims = self.forced_plot_lims

        ndim = len(varied_params)

        best_fit_params, best_fit_param_errs, best_fit_prob, weighted_mean_params, weighted_mean_diff_prob = self.computeBestFitResultsFromMCMC(varied_params, mcmc_sampler = mcmc_sampler,burn_in_steps = burn_in_steps, thinning = thinning,
                                                                                                                           percentiles_to_print = percentiles_to_print, n_sig_figs = n_sig_figs, forced_luminosity_offset = forced_luminosity_offset, verbose_weighted_means = verbose)
        #truth_vals = [weighted_mean_params[i] if varied_params[i][0:2] == 'mu' else self.truth_vals[varied_params[i]] for i in range(len(varied_params))]
        truth_vals = [0.0 if varied_params[i][0:2] == 'mu' else None for i in range(len(varied_params))]
        if mcmc_sampler == None:
            flat_samples = self.mcmc_sampler.get_chain(discard=burn_in_steps, thin=thinning, flat=True)
            samples = self.mcmc_sampler.get_chain(discard=burn_in_steps, thin=thinning, flat=False)
        else:
            flat_samples = mcmc_sampler.get_chain(discard=burn_in_steps, thin=thinning, flat=True)
            samples = mcmc_sampler.get_chain(discard=burn_in_steps, thin=thinning, flat=False)


        flat_samples = np.transpose(flat_samples)
        flat_sample_means = [np.mean(flat_sample) for flat_sample in flat_samples]
        flat_samples = np.array([np.array(flat_samples[i]) - flat_sample_means[i] if i in cols_to_center else np.array(flat_samples[i]) for i in range(len(flat_samples))])
        flat_samples = np.transpose(flat_samples)

        labels = [self.varied_param_plot_strs[param] for param in varied_params]
        if not(save_fig):
            fig_chains, axes = plt.subplots(ndim, figsize=(10, 7), sharex=True, squeeze = False)
            for i in range(ndim):
                ax = axes[i, 0]
                ax.plot(samples[:, :, i], "k", alpha=0.3)
                ax.set_xlim(0, len(samples))
                ax.set_ylabel(labels[i])
                ax.yaxis.set_label_coords(-0.1, 0.5)
            axes[-1, 0].set_xlabel("step number")

        fig_corner = corner.corner( flat_samples, bins = bins, labels=labels, truths = truth_vals, levels = fitter_levels,  label_kwargs = {'fontsize':labelsize}, title_kwargs  = {'fontsize':titlesize}, show_titles=True, titles = ['' for label in labels])
        axes_corner = np.array(fig_corner.axes).reshape((ndim, ndim))
        for i in range(len(varied_params)):
            param = varied_params[i]
            if param in list(forced_plot_lims.keys()):
                #Now need to loop through all cells in which this parameter is displayed
                lims = forced_plot_lims[param]
                x_plot_indeces = [[j, i] for j in range(i , len(varied_params))]
                y_plot_indeces = [ [i, j] for j in range(0, i) ]
                #print ('x_plot_indeces = ' + str(x_plot_indeces))
                #print ('y_plot_indeces = ' + str(y_plot_indeces))
                """
                for index_set in x_plot_indeces:
                    axes_corner[index_set[0], index_set[1]].set_xlim(lims)
                    #axes_corner[index_set[0], index_set[1]].text(0.5,0.5,'XLIMS: ' + str(index_set),fontsize = 20, horizontalalignment='center', verticalalignment='center', transform = axes_corner[index_set[0], index_set[1]].transAxes)
                for index_set in y_plot_indeces:
                    axes_corner[index_set[0], index_set[1]].set_ylim(lims)
                    #axes_corner[index_set[0], index_set[1]].text(0.5,0.25,'YLIMS: ' + str(index_set),fontsize = 20, horizontalalignment='center', verticalalignment='center', transform = axes_corner[index_set[0], index_set[1]].transAxes)
                """
        #fig_corner = corner.corner( flat_samples, labels=labels, truths = [self.truth_vals[var] for var in varied_params] )
        if overplot_basic_plot and len(self.params_to_overplot) > 1:
            """
            We need to make a new subplot in which we'll plot the contours.
            """
            n_overplot_dim = len(self.params_to_overplot)
            f_ref, axarr_ref = plt.subplots(n_overplot_dim - 1, n_overplot_dim - 1, squeeze = False, figsize = [overplot_elem_size * (len(self.params_to_overplot) - 1), overplot_elem_size * (len(self.params_to_overplot) - 1)])
            param_indeces_to_overplot = [varied_params.index(param) for param in self.params_to_overplot ]
            param_pairs_to_overplot = [[[self.params_to_overplot[j], self.params_to_overplot[i]] for i in range(j+1, len(param_indeces_to_overplot))] for j in range(len(param_indeces_to_overplot))]
            param_pairs_to_overplot = can.flattenListOfLists(param_pairs_to_overplot )
            print ('param_pairs_to_overplot = ' + str(param_pairs_to_overplot))
            axarr_indeces_to_overplot = [[[i - 1, j] for i in range(j+1, len(param_indeces_to_overplot))] for j in range(len(param_indeces_to_overplot))]
            axarr_indeces_to_overplot = can.flattenListOfLists(axarr_indeces_to_overplot )
            ref_flat_mcmc_samples = [np.transpose(self.basic_mcmc_sampler.get_chain(discard=burn_in_steps, thin=thinning, flat=True))[index] for index in param_indeces_to_overplot]
            #full_flat_mcmc_samples = np.transpose(self.mcmc_sampler.get_chain(discard=burn_in_steps, thin=thinning, flat=True))
            full_flat_mcmc_samples = np.transpose(flat_samples)
            for k in range(len(axarr_indeces_to_overplot)):
                axarr_indeces = axarr_indeces_to_overplot[k]
                param_pair = param_pairs_to_overplot[k]
                print ('axarr_indeces = ' + str(axarr_indeces))
                print ('param_pair = ' + str(param_pair))
                ref_mcmc_samples_to_plot = [ ref_flat_mcmc_samples[self.params_to_overplot.index(param_pair[0])], ref_flat_mcmc_samples[self.params_to_overplot.index(param_pair[1])] ]
                full_mcmc_samples_to_plot = [ full_flat_mcmc_samples[varied_params.index(param_pair[0])], full_flat_mcmc_samples[varied_params.index(param_pair[1])] ]
                #print ('ref_mcmc_samples_to_plot = ' + str(ref_mcmc_samples_to_plot))
                xlims =  axes_corner[axarr_indeces[0]+1, axarr_indeces[1]].get_xlim()
                ylims =  axes_corner[axarr_indeces[0]+1, axarr_indeces[1]].get_ylim()
                xticks =  axes_corner[axarr_indeces[0]+1, axarr_indeces[1]].get_xticks()
                yticks =  axes_corner[axarr_indeces[0]+1, axarr_indeces[1]].get_yticks()
                x_mesh, y_mesh, ref_mcmc_mesh = getContourMeshFromMCMCChainComponents(ref_mcmc_samples_to_plot, xlims, ylims, bins )
                x_mesh, y_mesh, full_mcmc_mesh = getContourMeshFromMCMCChainComponents(full_mcmc_samples_to_plot, xlims, ylims, bins )
                sorted_ref_mcmc_samples = np.flip(np.sort(ref_mcmc_mesh.flatten()))
                sorted_full_mcmc_samples = np.flip(np.sort(full_mcmc_mesh.flatten()))
                ref_levels = determineContourLevelsFromMCMCPosterior(sorted_ref_mcmc_samples, background_fitter_levels, samples_already_sorted = 1)
                full_levels = determineContourLevelsFromMCMCPosterior(sorted_full_mcmc_samples, fitter_levels, samples_already_sorted = 1)
                #sorted_ref_mcmc_samples = np.sorted(ref_mcmc_samples_to_plot)
                #sorted_full_mcmc_samples = np.sorted(full_mcmc_samples_to_plot)
                #axes_corner[axarr_indeces[0], axarr_indeces[1]].text(0.5,0.5,'[' + str(axarr_indeces[0]) + ', ' + str(axarr_indeces[1]) + ']',horizontalalignment='center', verticalalignment='center', transform = axes_corner[axarr_indeces[0], axarr_indeces[1]].transAxes)
                #axes_corner[axarr_indeces[0], axarr_indeces[1]].text(0.3, -1.2, '!!!!!!!!!!', color = 'r')
                print ('[ref_levels, full_levels] = ' + str([ref_levels, full_levels]))
                axarr_ref[axarr_indeces[0], axarr_indeces[1]].contourf(x_mesh, y_mesh, ref_mcmc_mesh, levels = ref_levels, alpha = ref_plot_alpha, cmap = contourf_cmap)
                #axarr_ref[axarr_indeces[0], axarr_indeces[1]].contourf(x_mesh, y_mesh, ref_mcmc_mesh, alpha = ref_plot_alpha, cmap = contourf_cmap)
                axarr_ref[axarr_indeces[0], axarr_indeces[1]].contour(x_mesh, y_mesh, full_mcmc_mesh, levels = full_levels, alpha = ref_plot_alpha, colors = contour_c)
                #axarr_ref[axarr_indeces[0], axarr_indeces[1]].contour(x_mesh, y_mesh, full_mcmc_mesh, alpha = ref_plot_alpha, colors = contour_c)
                axarr_ref[axarr_indeces[0], axarr_indeces[1]].axvline(truth_vals[varied_params.index(param_pair[0])], linestyle = '--', color = 'r', alpha = 0.25)
                axarr_ref[axarr_indeces[0], axarr_indeces[1]].axhline(truth_vals[varied_params.index(param_pair[1])], linestyle = '--', color = 'r', alpha = 0.25)
                axarr_ref[axarr_indeces[0], axarr_indeces[1]].set_xlabel(param_pair[0])
                axarr_ref[axarr_indeces[0], axarr_indeces[1]].set_ylabel(param_pair[1])
                axarr_ref[axarr_indeces[0], axarr_indeces[1]].set_xticks(xticks)
                axarr_ref[axarr_indeces[0], axarr_indeces[1]].set_yticks(yticks)
                axarr_ref[axarr_indeces[0], axarr_indeces[1]].set_xlim(xlims)
                axarr_ref[axarr_indeces[0], axarr_indeces[1]].set_ylim(ylims)
            #axarr_ref.sup_title('Full contours - ref fit; Lines - full fit')
            plt.tight_layout()
                #axarr_ref[axarr_indeces[0], axarr_indeces[1]].set_xlim(xlims)
                #axarr_ref[axarr_indeces[0], axarr_indeces[1]].set_ylim(ylims)

        if save_fig:
            fig_corner.savefig(target_dir + full_fig_name)
            if overplot_basic_plot and len(self.params_to_overplot) > 1:
                f_ref.savefig(target_dir + ref_fig_name)
        if show_fig:
            plt.show()
        else:
            plt.close('all')
        if save_single_chain:
            print ('Here 1 !!! single_chain_title: ' + single_chain_title  )
            if save_fig:
                n_single_chain_rows = int(np.ceil(len(varied_params) / n_single_chain_cols))
                f_s_chain, axarr_s_chain = plt.subplots(n_single_chain_rows, n_single_chain_cols, figsize = (5 * n_single_chain_cols, 2 * n_single_chain_rows ), squeeze = False)
                full_samples = self.mcmc_sampler.get_chain()
                for i in range(len(varied_params)):
                    axarr_s_chain[i // n_single_chain_cols, i % n_single_chain_cols].plot(full_samples[:, :, i], 'k', alpha = 0.3)
                    axarr_s_chain[i // n_single_chain_cols, i % n_single_chain_cols].set_xlim(0, len(full_samples))
                    axarr_s_chain[i // n_single_chain_cols, i % n_single_chain_cols].set_ylabel(varied_params[i])
                print ('target_dir + single_chain_title = ' + str(target_dir + single_chain_title) )
                plt.savefig(target_dir + single_chain_title )

        return 1

    # (self.sorted_mus, self.canon_mus, self.sorted_muErrs, self.sorted_surveys)
    def calculateWeightedMeansBySurvey(self, meas_mus, calc_mus, mu_errs, varied_surveys, ceph_muErrs, ceph_surveys, n_sig_figs = 5, verbose = 1):
        """
        To calculate what the weighted mean offsets SHOULD BE, we need to minimize the chi square with respect to
           each offset.  This requires solving a system of linear equations since the derivative of chi_square
           with respect to each survey offset does depend on the survey offsets of the other surveys.
        """
        full_weighted_mean_diff = np.sum((np.array([calc_mus]) - np.array(meas_mus)) * np.array(mu_errs) ** -2.0)
        n_cepheid_indeces = len(self.cepheid_indeces)
        n_sn = len(meas_mus)
        n_varied = len(varied_surveys)
        weighted_mu_mean_diffs = {survey:self.default_mu_offset for survey in self.surveys}
        #Otherwise, we need to simultaneously solve an entangled set of equations.  Does python support matrix algebra? It does
        linear_eqns_to_solve = [ [ 0.0 for j in range(n_varied) ] for i in range(n_varied) ]
        linear_offsets = [0.0 for j in range(n_varied)]
        if verbose:
            print ('varied_surveys = ' + str(varied_surveys))
        for j in range(n_varied):
            survey = varied_surveys[j]
            cepheid_survey_indeces = [ i for i in range(len(ceph_surveys)) if ceph_surveys[i] == survey ]
            survey_indeces = [ k for k in range(n_sn) if self.nonCeph_surveys[k] == survey ]
            #survey_and_cepheid_indeces = can.intersection(self.cepheid_indeces, survey_indeces)
            #survey_not_cepheid_indeces = [index for index in survey_indeces if not(index in self.cepheid_indeces)]

            meas_mus_in_survey = np.array( [meas_mus[index] for index in survey_indeces ])
            calc_mus_in_survey = np.array( [calc_mus[index] for index in survey_indeces ])
            mu_errs_in_survey = np.array( [mu_errs[index] for index in survey_indeces ] )
            n_in_survey_in_hf = len(meas_mus_in_survey)

            if n_cepheid_indeces > 0:
                dLumOffsetDSurveyOffset = np.sum([ ceph_muErrs[index] ** -2.0 for index in cepheid_survey_indeces ]) / np.sum(np.array(ceph_muErrs) ** -2.0 )
            else:
                dLumOffsetDSurveyOffset = 0

            weighted_mean_diff_of_survey = np.sum( [ (calc_mus_in_survey[i] - meas_mus_in_survey[i]) * mu_errs_in_survey[i] ** -2.0 for i in range(n_in_survey_in_hf) ] )

            additive_term = weighted_mean_diff_of_survey + full_weighted_mean_diff * dLumOffsetDSurveyOffset

            #all_survey_offsets_terms = [ np.sum([dLumOffsetDSurveyOffset * (calc_mus[k] - meas_mus[k]) * mu_errs[k] ** -2.0 for k in range(len(meas_mus)) if self.sorted_surveys[k] == each_survey]) for each_survey in varied_surveys ]
            all_survey_offsets_terms = [ dLumOffsetDSurveyOffset *  np.sum([mu_errs[k] ** -2.0 for k in range(n_sn) if self.nonCeph_surveys[k] == each_survey]) for each_survey in varied_surveys ]
            lum_offset_coefficient = np.sum([mu_errs[i] ** -2.0 for i in survey_indeces]) + dLumOffsetDSurveyOffset * np.sum([mu_errs[i] ** -2.0 for i in range(n_sn)])
            if n_cepheid_indeces > 0:
                lum_offsets_by_survey = [np.sum([ ceph_muErrs[cepheid_index] ** -2.0 for cepheid_index in range(len(ceph_surveys)) if ceph_surveys[cepheid_index] == each_survey ]) / np.sum(np.array(ceph_muErrs) ** -2.0 ) for each_survey in varied_surveys ]
                #lum_offsets_by_survey = [np.sum([ mu_errs[cepheid_index] ** -2.0 for cepheid_index in self.cepheid_indeces if self.sorted_surveys[cepheid_index] == each_survey ]) / np.sum([ mu_errs[cepheid_index] ** -2.0 for cepheid_index in self.cepheid_indeces ]) for each_survey in varied_surveys ]
            else:
                lum_offsets_by_survey = [0.0 for each_survey in varied_surveys ]
            lum_offsets_terms = [ lum_offset_coefficient * lum_offset for lum_offset in lum_offsets_by_survey ]
            this_survey_terms = np.array([np.sum([mu_errs[k] ** -2.0 for k in survey_indeces]) if each_survey == survey else 0.0 for each_survey in varied_surveys])
            offset_coefficients = this_survey_terms + np.array(all_survey_offsets_terms) + np.array(lum_offsets_terms)
            if verbose:
                print ('For survey: ' + str(survey))
                print ('[can.round_to_n(coef, n_sig_figs) for coef in this_survey_terms] = ' + str([can.round_to_n(coef, n_sig_figs) for coef in this_survey_terms]))
                print ('[can.round_to_n(coef, n_sig_figs) for coef in all_survey_offsets_terms] = ' + str([can.round_to_n(coef, n_sig_figs) for coef in all_survey_offsets_terms]))
                print ('[can.round_to_n(coef, n_sig_figs) for coef in lum_offsets_terms] = ' + str([can.round_to_n(coef, n_sig_figs) for coef in lum_offsets_terms]))
                print ('[can.round_to_n(coef, n_sig_figs) for coef in offset_coefficients] = ' + str([can.round_to_n(coef, n_sig_figs) for coef in offset_coefficients]))
                print ('additive_term = ' + str(additive_term))
            linear_eqns_to_solve[j] = offset_coefficients
            linear_offsets[j] = -additive_term

        linear_eqns_to_solve = np.array(linear_eqns_to_solve)
        if verbose:
            print ('linear_eqns_to_solve (rounded to ' + str(n_sig_figs) + ' sig figs):')
            for j in range(n_varied):
                print ('Survey ' + str(varied_surveys[j]) + ': ' + ' + '.join([str(can.round_to_n(coef, n_sig_figs)) for coef in linear_eqns_to_solve[j]]) + ' + ' + str(can.round_to_n(linear_offsets[j], n_sig_figs)) + ' = 0')
        #print ('linear_offsets = ' + str(linear_offsets))
        if n_varied > 1:
            try:
                if verbose: print ('[linear_eqns_to_solve, linear_offsets] = ' + str([linear_eqns_to_solve, linear_offsets]) )
                solved_survey_offsets = np.linalg.solve(linear_eqns_to_solve, linear_offsets)
            except np.linalg.LinAlgError as err:
                if 'Singular matrix' in str(err):
                    solved_survey_offsets = [ 0.0 for survey in self.surveys ]
                else:
                    raise
            #print ('solved_survey_offsets = ' + str(solved_survey_offsets))
            for j in range(n_varied):
                survey = varied_surveys[j]
                weighted_mu_mean_diffs[survey] = solved_survey_offsets[j]
        elif n_varied > 0:
            offset = linear_offsets[0] / linear_eqns_to_solve[0][0]
            weighted_mu_mean_diffs[varied_surveys[0]] = offset
        if verbose: print ('weighted_mu_mean_diffs = ' + str(weighted_mu_mean_diffs))
        return weighted_mu_mean_diffs

    def __init__(self, w_of_funct_str = 'w0', default_mu_offset = 0.0,
                 root_dir = 'stubbs/variableMuFits/', plot_dir = 'plots/', fit_res_dir = 'mcmcResults/', dir_base = '',
                 cepheid_file_name = 'calibratorset.txt', randomize_sn = 0, force_rand_fit_to_background_cosmology = 0,
                 covariance_file = 'stubbs/SNIsotropyProject/OriginalSNDataFiles/Pantheon+SH0ES_122221_1.cov', #covmat_NOSYS.txt',
                 sn_toy_data_file = 'stubbs/variableMuFits/ArtificialSurveys/ArtificialSNe_C.csv',
                 survey_mu_correllation_file = 'stubbs/SNIsotropyProject/OriginalSNDataFiles/PantheonPlusZeropointCorrellationsBySurvey.txt',
                 canonical_cosmic_vals = {'H0':70.0, 'L':1.0, 'OmM':0.3, 'OmL':0.7, 'OmR': 0.0001, 'Om0': 1.0, 'w0':-1,'wa':0}, #These values of H0 and OmM are the best fit ones if we only allow them to vary in a cosmic fit. So they are good starting values
                 H0_bounds = [30.0, 110.0], OmegaM_bounds = [0.0, 1.0], OmegaLambda_bounds = [0.0, 1.0], Omega0_bounds = [0.5, 1.5], OmegaR_bounds = [0.0, 0.1],
                 sn_data_type = 'pantheon_plus', mu_offset_bounds = [-1.0, 1.0], forced_plot_lims = {'OmM':[0.05, 0.55], 'w0':[-1.45, -0.45], 'H0':[67.5, 72.5] },
                 init_guess_params_to_overplot = None, params_to_overplot = ['OmM','w0'], overplot_mcmc_steps = 5000, overplot_mcmc_chains = 10,
                 mu_prior_type = 'gauss', muOffsetPriors  = {'ASASSN':[0.0, 0.05], 'CANDELS':[0.0, 0.05], 'CFA1':[0.0, 0.05],  'CFA2':[0.0, 0.05], 'CFA3K':[0.0, 0.05],      'CFA3S':[0.0, 0.05], 'CFA4p1':[0.0, 0.05],
                                    'CFA4p2':[0.0, 0.05],  'CSP':[0.0, 0.05],   'DES':[0.0, 0.05],  'FOUNDATION':[0.0, 0.05], 'HST':[0.0, 0.05],   'KAIT':[0.0, 0.05], 'KAITM':[0.0, 0.05],
                                    'LOWZ':[0.0, 0.05],    'PS1MD':[0.0, 0.05], 'SDSS':[0.0, 0.05],  'SNAP':[0.0, 0.05],      'SNLS':[0.0, 0.05],  'SWIFT':[0.0, 0.05], 'SWIFTNEW':[0.0, 0.05]},
                 z_lims = [-np.inf, np.inf] ):
        H0, OmegaM, OmegaLambda, OmegaR, Omega0, w0, wa = [canonical_cosmic_vals[param_str] for param_str in ['H0', 'OmM', 'OmL', 'OmR', 'Om0', 'w0','wa']]
        print ('[H0, OmegaM, OmegaLambda, OmegaR, Omega0, w0, wa] = ' + str([H0, OmegaM, OmegaLambda, OmegaR, Omega0, w0, wa]))
        if init_guess_params_to_overplot == None and params_to_overplot != None:
            init_guess_params_to_overplot = [ canonical_cosmic_vals[param] for param in params_to_overplot ]
        self.dir_base = dir_base
        self.astro_arch = apa.AstronomicalParameterArchive()
        self.cosmo_arch = cpa.CosmologicalParameterArchive(H0 = H0, OmegaM = OmegaM, OmegaLambda = OmegaLambda, Omega0 = Omega0, OmegaR = OmegaR, params_source = 'pantheon')
        self.scalar_params = [self.cosmo_arch.getH0()[0], self.cosmo_arch.getOmegaM()[0], self.cosmo_arch.getOmegaLambda()[0], self.cosmo_arch.getOmegaR()[0], self.cosmo_arch.getOmega0()[0]]
        self.age_of_universe = self.cosmo_arch.getAgeOfUniverse( units = 'yr' )[0]
        self.forced_plot_lims = forced_plot_lims
        self.z_lims = z_lims
        self.root_dir = dir_base + root_dir
        self.plot_dir = plot_dir
        self.mcmc_results_dir = fit_res_dir
        self.default_w_of_funct = lambda zs: -1.0 if type(zs) in [int, float] else -1.0 + np.zeros(np.shape(zs))

        self.survey_mu_correllation_matrix = self.readInSurveyMuOffsetCorrellationMatrix(dir_base  + survey_mu_correllation_file)

        self.canonical_cosmic_vals = canonical_cosmic_vals
        self.cosmic_param_bounds = [H0_bounds, OmegaM_bounds, OmegaLambda_bounds, OmegaR_bounds, Omega0_bounds]

        self.H0_in_inv_Myr = self.cosmo_arch.getH0(units = 'year')[0] * 10.0 ** 6.0
        self.km_s_Mpc_in_Myr = self.astro_arch.getKmPerSToPcPerYr()
        self.cepheid_file_name = cepheid_file_name
        self.H0, self.OmegaM, self.OmegaLambda, self.Omega0, self.OmegaR = [self.cosmo_arch.getH0()[0], self.cosmo_arch.getOmegaM()[0], self.cosmo_arch.getOmegaLambda()[0], self.cosmo_arch.getOmega0()[0], self.cosmo_arch.getOmegaR()[0]]
        self.scalar_param_strs = ['H0','OmM', 'OmL', 'OmR', 'Om0']
        self.mu_offset_bounds = mu_offset_bounds
        self.sn_data_type = sn_data_type
        self.muOffsetPriors = muOffsetPriors
        self.mu_prior_type = mu_prior_type
        self.sn_toy_data_file = dir_base + sn_toy_data_file
        self.covariance_file = dir_base + covariance_file
        print ('self.sn_toy_data_file = ' + str(self.sn_toy_data_file))
        #if self.sn_data_type == 'art_pantheon':
        #    self.offsets_by_survey_order = ['CANDELS', 'CFA1', 'CFA2', 'CFA3K', 'CFA3S', 'CFA4p1', 'CFA4p2', 'CSP', 'DES', 'FOUNDATION', 'PS1MD', 'SDSS', 'SNLS']
        #elif self.sn_data_type in ['realplus', 'real_plus', 'realp', 'real_p', 'pantheon_plus']:
        #    self.offsets_by_survey_order = ['CANDELS', 'CFA1', 'CFA2', 'CFA3K', 'CFA3S', 'CFA4p1', 'CFA4p2', 'CSP', 'HST', 'PS1MD', 'SDSS', 'SNAP', 'SNLS', 'FOUNDATION', 'KAIT', 'SWIFT', 'SWIFTNEW', 'LOWZ', 'DES']
        #else:
        #    self.offsets_by_survey_order = ['CANDELS', 'CFA1', 'CFA2', 'CFA3K', 'CFA3S', 'CFA4p1', 'CFA4p2', 'CSP', 'HST', 'PS1MD', 'SDSS', 'SNAP', 'SNLS']

        self.default_mu_offset = default_mu_offset
        self.color_to_survey_dict = {}
        self.lookback_time_type = 'tau'

        self.truth_vals = canonical_cosmic_vals
        self.w_of_funct, self.default_w_params, self.w_param_bounds, self.w_id_strs  = self.getVariableWCosmicFunction(w_of_funct_str = w_of_funct_str)
        #self.sorted_sn, self.sorted_zs, self.sorted_lookBacks, self.sorted_mus, self.sorted_muErrs, self.sorted_resids, self.sorted_resids_noBySurveyZeroing, self.sorted_surveys, self.sorted_ids, self.sorted_plot_colors, self.cepheid_sne = self.readInSNData([0.0 for survey in self.offsets_by_survey_order ], randomize_sn = randomize_sn)
        #self.sorted_sn, self.sorted_zs, self.sorted_mus, self.sorted_muErrs, self.sorted_surveys, self.sorted_ids, self.sorted_plot_colors, self.cepheid_sne = self.readInSNData([0.0 for survey in self.offsets_by_survey_order ], randomize_sn = randomize_sn)
        self.full_sorted_zs, self.full_sorted_sn, self.full_sorted_muErrs, self.full_sorted_mus, self.full_sorted_surveys, self.full_sorted_ids, self.full_sorted_plot_colors, self.full_sorted_cov_matrix, self.cepheid_sne, self.full_cepheid_indeces = self.readInSNData(randomize_sn = randomize_sn, force_to_background_cosmology = force_rand_fit_to_background_cosmology)
        self.sorted_zs, self.sorted_sn, self.sorted_mus, self.sorted_muErrs, self.sorted_plotColors, self.sorted_surveys, self.sorted_ids, self.sorted_cov_matrix, self.cepheid_indeces = [self.full_sorted_zs[:], self.full_sorted_sn[:], self.full_sorted_mus[:], self.full_sorted_muErrs[:], self.full_sorted_plot_colors[:], self.full_sorted_surveys[:], self.full_sorted_ids[:], self.full_sorted_cov_matrix[:], self.full_cepheid_indeces[:]]
        self.sorted_inv_cov_matrix = np.linalg.inv(self.sorted_cov_matrix)
        self.varied_param_plot_strs = {'H0':r'$\Delta H_0$', 'OmM':r'$\Delta \Omega_M$', 'w0':r'$\Delta w$', 'wa':r'$\Delta w_a$'}
        for survey in self.sorted_surveys:
            self.varied_param_plot_strs['mu' + survey] = r'$\Delta \mu$' + ''.join([r'$_{}$'.format(char) for char in survey] )
        self.non_cepheid_indeces = [index for index in range(len(self.sorted_sn)) if not(index in self.cepheid_indeces)]
        self.sorted_mu_weights = 1.0 / np.array(self.sorted_muErrs) ** 2.0
        self.nonCeph_zs, self.nonCeph_sn, self.nonCeph_mus, self.nonCeph_muErrs, self.nonCeph_mu_weights, self.nonCeph_plotColors, self.nonCeph_surveys, self.nonCeph_ids = [[arr[index] for index in self.non_cepheid_indeces] for arr in [self.sorted_zs, self.sorted_sn, self.sorted_mus, self.sorted_muErrs, self.sorted_mu_weights, self.sorted_plotColors, self.sorted_surveys, self.sorted_ids]]
        self.nonCeph_cov_matrix = [[self.sorted_cov_matrix[col_index][row_index] for row_index in self.non_cepheid_indeces] for col_index in self.non_cepheid_indeces]
        self.nonCeph_inv_cov_matrix = np.linalg.inv(self.nonCeph_cov_matrix)
        self.randomize_sn = randomize_sn
        #self.cepheid_indeces = [index for index in self.cepheid_indeces if self.sorted_surveys[index] == 'SWIFT'] #To debug, we want to make no correction for the cepheid indeces
        #self.cepheid_sne = [self.sorted_ids[index] for index in self.cepheid_indeces]
        self.cepheid_surveys = np.unique([self.sorted_surveys[index] for index in self.cepheid_indeces  ])
        self.full_cepheid_surveys = self.cepheid_surveys[:]
        unique_surveys = np.unique(self.sorted_surveys)
        mean_redshifts_by_survey = [np.mean([self.sorted_zs[i] for i in range(len(self.sorted_zs)) if self.sorted_surveys[i] == survey]) for survey in unique_surveys]
        self.surveys = can.safeSortOneListByAnother(mean_redshifts_by_survey, [unique_surveys])[0]
        self.all_surveys = np.unique(self.full_sorted_surveys)
        self.n_surveys = len(self.surveys)
        for survey in self.surveys:
            self.truth_vals['mu' + survey] =  default_mu_offset

        print ('[self.H0, self.OmegaM, self.OmegaLambda, self.OmegaR , self.Omega0, self.default_w_params] = ' + str([self.H0, self.OmegaM, self.OmegaLambda, self.OmegaR , self.Omega0, self.default_w_params]))
        self.canon_mus = getMusForCosmology(self.nonCeph_zs, [self.H0, self.OmegaM, self.OmegaLambda, self.OmegaR , self.Omega0], lambda zs: self.w_of_funct(zs, self.default_w_params), astro_archive = self.astro_arch, cosmo_archive = self.cosmo_arch) #calcMuForw.ResidualMuCalculatorForArbitraryWofT(wOfFunction = self.w_of_funct, initial_zs = self.sorted_zs).getMus()
        self.null_chi_square = np.sum((np.array(self.nonCeph_mus) - np.array(self.canon_mus)) ** 2.0 / np.array(self.nonCeph_muErrs) ** 2.0  )

        self.init_guess_params_to_overplot = init_guess_params_to_overplot
        self.params_to_overplot = params_to_overplot
        self.overplot_mcmc_steps = overplot_mcmc_steps
        self.overplot_mcmc_chains = overplot_mcmc_chains
        if params_to_overplot != None and len(params_to_overplot) > 0:
            print ('Running fitting for basic cosmological model...')
            self.params_to_overplot = params_to_overplot
            print ('[init_guess_params_to_overplot, params_to_fit_and_overplot] = ' + str([init_guess_params_to_overplot, params_to_overplot] ))
            self.basic_mcmc_sampler = self.doFullCosmicFit( init_guess_params_to_overplot, params_to_overplot, do_init_minim = 0, n_mcmc_steps = overplot_mcmc_steps, n_mcmc_chains = overplot_mcmc_chains, do_all_mu_offsets = 0, overwrite_mcmc = 0, show_mcmc = 0, save_mcmc = 0, overplot_basic_plot = 0, verbose = 0)
        else:
            self.params_to_overplot = []
            self.basic_mcmc_sampler = None

if __name__ == "__main__" :
    """
    Make plot of Hubble constant uncertainties based on supernovae offset residuals by survey.
    """
    cl_args = sys.argv[1:]
    save_files_suffix = cl_args[0]
    print ('Starting to fit Hubble constant from the command line... ')
    dir_base = '/Users/sashabrownsberger/Documents/Harvard/physics/'
    results_dir = dir_base  + 'stubbs/variableMuFits/mcmcResults/'
    plot_dir = dir_base  + 'stubbs/variableMuFits/plots/PosteriorWidths/'
    #mu_offsets = [0.001, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6]
    #mu_offsets = [0.35, 0.45]
    #mu_offsets = [0.01, 0.25, 0.3, 0.4, 0.5, 0.6]
    mu_offsets = [0.001, 0.014, 0.018, 0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5,0.55, 0.6]
    #mu_offsets = [0.002, 0.005, 0.008, 0.015, 0.025, 0.03, 0.035]
    #mu_offsets = [0.001, 0.018, 0.05]
    ref_survey = 'PS1MD'
    #Mostly pulled from Table 3 of https://arxiv.org/pdf/2112.03864.pdf
    mu_offsets_scalings_dict = {'SWIFT':np.mean([0.012, 0.011]),
                                'ASASSN':np.mean([0.022, 0.021, 0.02, 0.021, 0.020, 0.022, 0.021, 0.021]),
                                'CFA1':np.mean([0.012, 0.01, 0.011, 0.011]), #Set sams as CFA3S
                                'CFA2':np.mean([0.012, 0.01, 0.011, 0.011]),
                                'LOWZ':np.mean([0.05]), #From https://iopscience.iop.org/article/10.1086/512054/pdf, "estimate the intrinsic dispersion about the relation is 0.05 mag"
                                'KAITM':np.mean([0.012, 0.011, 0.010, 0.011, 0.012, 0.010, 0.011, 0.010, 0.013, 0.011, 0.010, 0.011, 0.012, 0.010, 0.010, 0.011]),
                                'CFA4p2':np.mean([0.011, 0.011, 0.010, 0.011]),
                                'KAIT':np.mean([0.012, 0.011, 0.010, 0.011, 0.012, 0.010, 0.011, 0.01, 0.013, 0.011, 0.010, 0.011, 0.012, 0.010, 0.010, 0.011]),
                                'CFA3S':np.mean([0.012, 0.01, 0.011, 0.011]),
                                'CSP':np.mean([0.011, 0.010, 0.010, 0.011, 0.011, 0.011, 0.011, 0.011]),
                                'CFA3K':np.mean([0.012, 0.010, 0.012, 0.010]),
                                'CFA4p1':np.mean([0.012, 0.010, 0.011, 0.012]),
                                'PS1MD':np.mean([0.006, 0.006, 0.006, 0.006]),
                                'SDSS':np.mean([0.006, 0.005, 0.006, 0.006]),
                                'DES':np.mean([0.006, 0.006, 0.006, 0.006]),
                                'SNLS':np.mean([0.005, 0.005, 0.005, 0.006]),
                                'HST':np.mean([0.006, 0.006, 0.006, 0.006]),
                                'SNAP':np.mean([0.006, 0.006, 0.006, 0.006]),
                                'CANDELS':np.mean([0.006, 0.006, 0.006, 0.006])}
    mu_offset_normalization = mu_offsets_scalings_dict[ref_survey]
    mu_offsets_scalings_dict = {key:mu_offsets_scalings_dict[key] / mu_offset_normalization for key in mu_offsets_scalings_dict.keys()}
    print ('mu_offset_scalings_dict = ' + str(mu_offsets_scalings_dict))
    #mu_offsets = [0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6]
    #mu_offsets = [0.001, 0.01, 0.04, 0.08, 0.1, 0.14, 0.18, 0.2, 0.3, 0.4, 0.5, 0.6]
    #mu_offsets = [0.1]
    #mu_offsets = [0.001, 0.1]
    z_lims = [0.0, np.inf]
    z_lims_str = 'zLim_Full'
    #SH0ES redshift cuts; if you use these, you should probably set OmM = 0.3, w = -1
    #z_lims = [0.023, 0.15]
    #z_lims_str = 'zLim_0p023_0p15'
    H0_plot_file_name = 'H0UnertaintyVsMuOffsetPrior_' + z_lims_str + '.pdf'
    OmM_plot_file_name = 'OmMUnertaintyVsMuOffsetPrior_' + z_lims_str + '.pdf'
    w0_plot_file_name = 'w0UnertaintyVsMuOffsetPrior_' + z_lims_str + '.pdf'
    save_file = 'PosteriorFitsDifferentOffsets_' + z_lims_str + '_' + save_files_suffix + '.txt'
    plot_file_name = 'CosmicParamsUncertaintyVsMuOffsetPrior_' + z_lims_str + '_' + save_files_suffix + '.pdf'
    #mu_offsets = [0.0002, 0.2]
    posterior_fits = [0.0 for i in range(len(mu_offsets))]
    overplot_steps = 10000
    overplot_steps = 10000
    overplot_n_chains = 6
    full_steps = 10000
    #full_steps = 1500
    full_n_chains = 50
    #full_n_chains = 10
    #surveys = ['CANDELS', 'CFA1',  'CFA2', 'CFA3K', 'CFA3S', 'CFA4p1', 'CFA4p2', 'CSP', 'DES', 'FOUNDATION', 'HST', 'KAIT', 'LOWZ', 'PS1MD', 'SDSS', 'SNAP', 'SNLS', 'SWIFT', 'SWIFTNEW']
    cosmic_params_to_vary = ['H0', 'OmM', 'w0']
    cosmic_params_to_vary_init_guess = [70.0, 0.3, -1.0]
    #cosmic_params_to_vary = ['H0','OmM']
    #cosmic_params_to_vary_init_guess = [70.0, 0.3]
    #surveys = ['CANDELS', 'CFA1', 'CFA2', 'CFA3K', 'CFA3S', 'CFA4p1', 'CFA4p2', 'CSP', 'DES', 'HST', 'KAIT', 'LOWZ', 'PS1MD', 'SDSS', 'SNAP', 'SNLS', 'SWIFT']
    #{'CANDELS': 1.5791014285714287, 'CFA1': 0.023386, 'CFA2': 0.019155937499999998, 'CFA3K': 0.03073, 'CFA3S': 0.025453863636363636, 'CFA4p1': 0.03031708333333333, 'CFA4p2': 0.021545384615384616, 'CSP': 0.02634633663366336, 'DES': 0.3892227979274612, 'FOUNDATION': 0.03817131428571429, 'HST': 1.0970235294117647, 'KAIT': 0.022085254237288136, 'LOWZ': 0.0233895, 'PS1MD': 0.2865038078291815, 'SDSS': 0.20657667525773193, 'SNAP': 1.2256033333333334, 'SNLS': 0.6388950854700854, 'SWIFT': 0.005737, 'SWIFTNEW': 0.0155971875}
    surveys = ['SWIFT', 'ASASSN', 'CFA2', 'CFA1', 'KAITM', 'LOWZ', 'CFA4p2', 'KAIT', 'CFA3S', 'CSP', 'CFA3K', 'CFA4p1', 'PS1MD', 'SDSS', 'DES', 'SNLS', 'HST', 'SNAP', 'CANDELS']
    #surveys = ['CANDELS', 'CFA1' ]
    mu_prior_type = 'gauss'
    #surveys = ['SDSS']
    rand_fitter = CosmicFitter(w_of_funct_str = 'w0', randomize_sn = 0, sn_data_type = 'pantheon_plus', params_to_overplot = cosmic_params_to_vary, overplot_mcmc_steps = overplot_steps, overplot_mcmc_chains = overplot_n_chains, mu_prior_type = mu_prior_type, dir_base = dir_base)
    #rand_fitter = CosmicFitter(w_of_funct_str = 'w0', randomize_sn = 1, sn_data_type = 'pantheon_plus', params_to_overplot = None, overplot_mcmc_steps = overplot_steps, overplot_mcmc_chains = overplot_n_chains, mu_prior_type = mu_prior_type)
    rand_fitter.updateUsedSN(z_lims = z_lims, surveys_to_include = surveys )
    surveys = list(rand_fitter.surveys)
    for i in range(len(mu_offsets)):
        mu_offset = mu_offsets[i]
        for survey in surveys:
            mu_offset_scaling = mu_offsets_scalings_dict[survey]
            rand_fitter.muOffsetPriors[survey] = [0.0, mu_offset * mu_offset_scaling]
        rand_fitter.doFullCosmicFit(cosmic_params_to_vary_init_guess + [0.0 for survey in surveys], cosmic_params_to_vary + ['mu' + survey for survey in surveys], n_mcmc_steps = full_steps, n_mcmc_chains = full_n_chains , verbose = 0, additional_save_prefix = 'GaussPriorWidth' + str(int(1000 * mu_offset)) + 'mMags_' + z_lims_str + '_', additional_save_suffix = '_' + save_files_suffix, save_mcmc = 1, show_mcmc = 0, save_full_mcmc = 1 )
        posterior_fits[i] = rand_fitter.computeBestFitResultsFromMCMC(cosmic_params_to_vary + ['mu' + survey for survey in surveys], verbose_weighted_means = 0)

    cols_to_save = [mu_offsets]
    print ('posterior_fits = ' + str(posterior_fits))
    for i in range(len(cosmic_params_to_vary + surveys)):
        print ('i = ' + str(i))
        cols_to_save = cols_to_save + [[can.round_to_n(posterior_fit[0][i], 6) for posterior_fit in posterior_fits] , [can.round_to_n(posterior_fit[1][i][0],6) for posterior_fit in posterior_fits], [can.round_to_n(posterior_fit[1][i][1], 6) for posterior_fit in posterior_fits]]
    print ('cols_to_save = ' + str(cols_to_save))
    header = 'DeltaMu, ' + ', '.join([param + 'mu, -' + param + 'sig, +' + param + 'sig' for param in ['H0', 'OmM'] + ['mu' + survey for survey in surveys]])
    can.saveListsToColumns(cols_to_save, save_file, results_dir, sep = ', ', append = False, header = header, type_casts = None)

    """
    f, axarr = plt.subplots(3, 1, figsize = (5,6))
    H0_ax = axarr[0]
    Om_ax = axarr[1]
    w0_ax = axarr[2]
    makePlotOfPosteriorWidths([save_file], 1, ax = H0_ax, xlabel = '$\Delta \mu_S$ Gaussian prior width (mags)', ylabel = r'$H_0$ MCMC $\sigma$ (km/s/Mpc)',)
    makePlotOfPosteriorWidths([save_file], 2, ax = Om_ax, xlabel = '$\Delta \mu_S$ Gaussian prior width (mags)', ylabel = r'$\Omega_M$ MCMC $\sigma$',)
    makePlotOfPosteriorWidths([save_file], 3, ax = w0_ax, xlabel = '$\Delta \mu_S$ Gaussian prior width (mags)', ylabel = r'$w_{\Lambda}$ MCMC $\sigma$',)
    plt.subplots_adjust (wspace=0.5, hspace=0.5)

    plt.savefig(plot_dir +  plot_file_name)
    """

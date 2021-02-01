import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np
import math
from loadSN import loadSN
import cantrips as c
import matplotlib.pyplot as plt
from CosmologicalParameterArchive import CosmologicalParameterArchive
import scipy.integrate as integrate
import scipy.optimize as optimize
import scipy.interpolate as interpolate
import binData as bd
import time
import randomSortData as rsd
import os
import sys
import matplotlib.gridspec as gridspec
import matplotlib
import cantrips as c


def getDirs():
    plt_save_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/SNTwoPointCorrelationsProject/plots/'
    res_save_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/SNTwoPointCorrelationsProject/fits/'
    return [plt_save_dir, res_save_dir]

def plotRChiSqr(res_file, res_save_dir = None):
    if res_save_dir is None:
        plt_save_dir, res_save_dir = getDirs()
    raw_fit_lines = c.readInFileLineByLine(res_file, file_dir = res_save_dir, n_ignore = 0)
    n_lines = len(raw_fit_lines)
    real_rchisqrs = 0
    rand_rchisqrs = [[] for line in raw_fit_lines]
    real_fits = []
    rand_fits = [[] for line in raw_fit_lines]
    for i in range(n_lines):
        line = raw_fit_lines[i]
        info, fit_str = line.split('=')
        info = info.strip()
        fit_str = fit_str.strip()

        print ('[info, fit_str] = ' + str([info, fit_str]))
        fit = []
        #for fit_elem in fit_str.split(','):
        #    if fit_
        #real_or_rand = info.split(':')[0]
        #if real_or_rand.lower() in ['real']:
        #    if real_rchisqrs > 0:
        #        print ('Multiple measurements on real data found on this fit!  That likely means multiple (possibly different) fits on non-randomized data wrote to this same file. We keep only the most recent.')
        #    real
    return 1

def functToMinimize(fitting_funct, autocorrs, autocorr_errs, comoving_perps, comoving_parallels, dof, param_scalings, param_functs, verbose = 0, ref_rchisqr = 1.0 ):
    #f, axarr = plt.subplots(2,1)
    params = [param_functs[i](param_scalings[i]) for i in range(len(param_functs))]
    funct_output = fitting_funct(np.array(comoving_perps), np.array(comoving_parallels), params )
    #axarr[0].scatter(comoving_perps, funct_output)
    #axarr[1].scatter(comoving_parallels, funct_output)
    #plt.show()
    shift = c.weighted_mean(funct_output - np.array(autocorrs), autocorr_errs)
    #print ('[funct_output, shift] = ' + str([funct_output, shift]))
    #print ('shift = ' + str(shift))
    rchisqr = np.sum(((funct_output - shift) - np.array(autocorrs)) ** 2.0 / (np.array(autocorr_errs) ** 2.0)) / dof

    #if verbose: print ('[param_scalings, params, shift, dof, ref_rchisqr, rchisqr / ref_rchisqr] = ' + str([param_scalings, params, shift, dof, ref_rchisqr, rchisqr / ref_rchisqr]))

    return [rchisqr / ref_rchisqr, np.append(np.array(param_scalings), shift) ]


def runMinimization(fig_name = None, fit_results_file_name = None, fitting_funct_str = 'polynomial_decay',
                    verbose = 1, save_fig = 0, show_fig = 1, save_fit = 0, extinction_correction_functs = [lambda ext: 0.0],
                    n_z_bins = 2, autocorr_to_plot = 'dl', surveys_to_load = ['all'], surveys_to_excise = [], randomize = 0,
                    plt_save_dir = None, res_save_dir = None, plot_perp_parallel = 1, zHD = 1, only_good_ext = 0, z_bin_ranges = [[0.0, 0.2], [0.2, np.inf]] ):
    if verbose: print ('Starting minimization of correlations...')
    start = time.time()
    if len(extinction_correction_functs) == 1: extinction_correction_functs = [extinction_correction_functs[0] for z_bin in range(n_z_bins)]
    # autocorr_to_plot can be 'dl','ext','mu','dl-ext', 'dl,ext' 'dl,dl-dl,ext'


    fig_name_prefixes = {'dl':'DeltaDLCorrelations_',
                         'ext':'ExtinctionCorrelations_',
                         'mu':'DeltaMuCorrelations_',
                         'dl-ext':'DeltaDLSubExtinctionCorrelations',
                         'dl,ext':'DeltaDLExtinctionsCrossCorrelations',
                         'dl,dl-dl,ext':'DeltaDLCorrSubDeltaDLExtinctionsCrossCorrelations'}
    extintion_cut_str = ('goodExtinctions' if only_good_ext else 'allSN')
    zHD_str = ('zHD' if zHD else 'zCMB')
    if fig_name is None:
        fig_name = fig_name_prefixes[autocorr_to_plot] + '_gaussSmoothing_' + extintion_cut_str + 'PantheonDataReal_TwoZBins_Even_fullFrame_' + fitting_funct_str + zHD_str + '.png'

    if fit_results_file_name is None:

        fit_results_file_name = fig_name_prefixes[autocorr_to_plot] + '_gaussSmoothing_' + extintion_cut_str + 'PantheonDataReal_TwoZBins_' + fitting_funct_str + zHD_str +'.txt'


    if plt_save_dir is None:
        plt_save_dir, dummy = getDirs()
    if res_save_dir is None:
        dummy, res_save_dir = getDirs()

    smoothing = 'gauss' # 'rect_bin', 'gauss'

    colorbar_labels = {'dl':r'Mean $(\delta D_{L,i} \times \delta D_{L,j})$ in bin',
                       'mu':r'Mean $(\delta \mu_{i} \times \delta \mu_{j})$ in bin',
                       'ext':r'Mean $(ext_i \times ext_j)$ in bin',
                       'dl-ext': r'Mean $((D_{L_i} \times D_{L_j}) - (ext_i \times ext_j))$ in bin' ,
                       'dl,ext':r'Mean $((D_{L_i} \times ext_{j}))$ in bin',
                       'dl,dl-dl,ext':r'Mean $((D_{L_i}\times D_{L_j} - D_{L_i} \times ext_{j}))$ in bin'}
    colorbar_label = colorbar_labels[autocorr_to_plot]

    all_sn = loadSN(1, surveys_to_load, pull_extinctions = 0, zHD = zHD)
    if only_good_ext:
        all_sn = [sn for sn in all_sn if not(np.isnan(sn['extinction'])) ]
    all_sn = [sn for sn in all_sn if not(sn['survey'] in surveys_to_excise)]
    if verbose: print ('len(all_sn) = ' + str(len(all_sn)))
    all_zHDs = [sn['zHD'] for sn in all_sn]
    all_zCMBs = [sn['zCMB'] for sn in all_sn]
    all_zs = [sn['z'] for sn in all_sn]
    print ('len(all_zs) = ' + str(len(all_zs)))
    #all_zs = [sn['zHD'] if sn['zHD'] <= 0.08 else sn['zCMB'] for sn in all_sn]
    print ('all_zs[0] = ' + str(all_zs[0]))
    sorted_zs, sorted_sn = c.safeSortOneListByAnother(all_zs, [all_zs, all_sn])
    mid_z = sorted_zs[len(sorted_zs) // 2]
    min_z = np.min(sorted_zs)
    max_z = np.max(sorted_zs)

    cosmo_arch = CosmologicalParameterArchive()
    OmM = cosmo_arch.getOmegaM()[0]
    OmL = cosmo_arch.getOmegaLambda()[0]
    OmR = cosmo_arch.getOmegaR()[0]
    H0 = cosmo_arch.getH0()[0]
    h = H0 / 100
    speedoflight = cosmo_arch.getc()


    n_stds_for_levels = {'dl':[0.3 for z_bin_range in z_bin_ranges],
                         'ext':[0.8 for z_bin_range in z_bin_ranges],
                         'mu':[0.8 for z_bin_range in z_bin_ranges],
                         'dl-ext':[0.3 for z_bin_range in z_bin_ranges],
                         'dl,ext':[0.5 for z_bin_range in z_bin_ranges],
                         'dl,dl-dl,ext':[0.5 for z_bin_range in z_bin_ranges]}
    n_std_for_levels = n_stds_for_levels[autocorr_to_plot]
    x_lims_in_zranges = [[0.0, 800 / h], [0.0, 800 / h]]
    x_lims_in_zranges = [[0.0, 'max'], [0.0, 'max']]
    y_lims_in_zranges = [[0.0, 200 / h], [0.0, 200 / h]]
    y_lims_in_zranges = [[0.0, 'max'], [0.0, 'max']]
    n_x_pixels = 50
    #n_comoving_bins = [50, int(50 * max_comoving_parallels / max_comoving_perps)]
    #comoving_bin_size = [(max_comoving_perps - 0.0 ) / n_comoving_bins[0], (max_comoving_parallels - 0.0 ) / n_comoving_bins[1]]

    if save_fig or show_fig:
        figsize = (3.5 * len(z_bin_ranges), 6)
        f = plt.figure(figsize=figsize)
        gs_z_bins = gridspec.GridSpec(1, len(z_bin_ranges), figure=f, bottom=9.0, top = 10.0)

    min_rchisqrs = [0.0 for z_bin_num in range(len(z_bin_ranges))]
    min_param_scalings_sets = [[] for z_bin_num in range(len(z_bin_ranges))]
    min_params_sets = [[] for z_bin_num in range(len(z_bin_ranges))]
    min_shifts = [0.0 for z_bin_num in range(len(z_bin_ranges))]
    for z_bin_num in range(len(z_bin_ranges)):

        z_bin_range = z_bin_ranges[z_bin_num]
        extinction_correction_funct = extinction_correction_functs[z_bin_num]

        if verbose: print ('z_bin_range = ' + str(z_bin_range))
        #binned_sn = [sn for sn in all_sn if sn['z'] >= z_bin_range[0] and sn['z'] < z_bin_range[1]]
        binned_sn = [ sn for sn in all_sn if sn['z'] >= z_bin_range[0] and sn['z'] < z_bin_range[1] ]
        if verbose: print ('len(binned_sn) = ' + str(len(binned_sn)))
        if len(binned_sn) <= 1: break
        binned_zs = [sn['z'] for sn in binned_sn]
        binned_exts = [sn['extinction'] for sn in binned_sn]
        binned_muDiffs = [sn['muDiff'] for sn in binned_sn]
        binned_muDiffs = [binned_muDiffs[i] - extinction_correction_funct(binned_exts[i]) for i in range(len(binned_muDiffs))]
        binned_muErrs = [sn['muErr'] for sn in binned_sn]
        #if randomize:
        #    rand_binned_muDiffs, rand_binned_muErrs = rsd.randomShuffleListOfLists([binned_muDiffs, binned_muErrs])
        #    binned_muDiffs, binned_muErrs = [rand_binned_muDiffs, rand_binned_muErrs]
        binned_deltaDLs = [10.0 ** (binned_muDiff / 5.0) - 1 for binned_muDiff in binned_muDiffs]
        binned_deltaDLErrs = [abs(binned_muErrs[i] * 10.0 ** (binned_muDiffs[i] / 5.0) * np.log(10) / 5.0) for i in range(len(binned_deltaDLs))]
        binned_RAs = [sn['RA'] for sn in binned_sn]
        binned_Decs = [sn['Dec'] for sn in binned_sn]
        binned_surveys = [sn['survey'] for sn in binned_sn]
        sorted_binned_exts = np.array(sorted(binned_exts))
        #nan_indeces = [index for index in range(len(binned_exts)) if np.isnan(binned_exts[index])]

        print ('[speedoflight, H0] = ' + str([speedoflight, H0] ))
        comoving_dist = [speedoflight / H0 * integrate.quad(lambda z_int: 1.0 / math.sqrt( OmM * (1.0 + z_int) ** 3.0 + OmL + OmR * (1.0 + z_int) ** 4.0 ) , 0 , z)[0] for z in binned_zs ]
        comoving_centers = c.flattenListOfLists([[ 0.5 * (comoving_dist[i] ** 2.0
                                                          + 2.0 * comoving_dist[i] * comoving_dist[j] * np.cos(c.measureAngularSeparationOnSky([binned_RAs[i], binned_Decs[i]], [binned_RAs[j], binned_Decs[j]], return_radian = 1) )
                                                          + comoving_dist[j] ** 2.0) ** 0.5
                                                   for j in range(i+1, len(comoving_dist)) ] for i in range(len(comoving_dist))],  )

        delta_DL_products = c.flattenListOfLists([[binned_deltaDLs[i] * binned_deltaDLs[j] for j in range(i+1, len(binned_deltaDLs)) ] for i in range(len(binned_deltaDLs))] )
        delta_DL_product_errs = c.flattenListOfLists([[ math.sqrt((binned_deltaDLs[i] * binned_deltaDLErrs[j]) ** 2.0 + (binned_deltaDLs[j] * binned_deltaDLErrs[i]) ** 2.0)  for j in range(i+1, len(binned_muDiffs)) ] for i in range(len(binned_muDiffs))] )
        delta_mu_product_norms = c.flattenListOfLists([[(binned_muDiffs[i] * binned_muDiffs[j]) / (binned_muErrs[i] * binned_muErrs[j]) for j in range(i+1, len(binned_muDiffs)) ] for i in range(len(binned_muDiffs))] )
        ext_product_norms = c.flattenListOfLists([[(binned_exts[i] * binned_exts[j]) for j in range(i+1, len(binned_exts)) ] for i in range(len(binned_exts))] )
        ext_delta_DL_cross_products = c.flattenListOfLists([[(binned_exts[i] * binned_deltaDLs[j]) for j in range(i+1, len(binned_exts)) ] for i in range(len(binned_exts))] )
        ext_delta_mu_cross_products = c.flattenListOfLists([[(binned_exts[i] * binned_muDiffs[j]) for j in range(i+1, len(binned_exts)) ] for i in range(len(binned_exts))] )
        ang_seps = c.flattenListOfLists([[c.measureAngularSeparationOnSky([binned_RAs[i], binned_Decs[i]], [binned_RAs[j], binned_Decs[j]], return_radian = 1) for j in range(i+1, len(binned_RAs))] for i in range(len(binned_RAs))])

        if autocorr_to_plot in ['dl']: # 'dl','ext','mu','dl-ext'
            autocorrs = delta_DL_products
            autocorr_errs = np.zeros(np.shape(autocorrs)) + 1.0
        elif autocorr_to_plot in ['ext']: # 'dl','ext','mu','dl-ext'
            autocorrs = ext_product_norms
            autocorr_errs = np.zeros(np.shape(autocorrs)) + 1.0
        elif autocorr_to_plot in ['mu']: # 'dl','ext','mu','dl-ext'
            autocorrs = delta_mu_product_norms # could also be delta_mu_product_norms, delta_DL_products, ext_product_norms
            #autocorr_errs = delta_mu_product_norms
            autocorr_errs = np.zeros(np.shape(autocorrs)) + 1.0
        elif autocorr_to_plot in ['dl-ext']: # 'dl','ext','mu','dl-ext'
            dl_ext_product_scaling = np.sum(np.array(ext_product_norms) * np.array(ext_product_norms) / (np.zeros(np.shape(delta_DL_products)) + 1.0)) / np.sum(np.array(ext_product_norms) ** 2.0 / (np.zeros(np.shape(delta_DL_products)) + 1.0)) #To minimze sum_i ((a_i - A * b_i) / sig_i) ** 2.0, set A = (sum_i a_i * b_i / sig_i ** 2.0) / (sum_i b_i * b_i / sig_i ** 2.0)
            print ('dl_ext_product_scaling = ' + str(dl_ext_product_scaling))
            autocorrs = np.array(delta_DL_products) - np.array(ext_product_norms) * dl_ext_product_scaling
            autocorr_errs = np.zeros(np.shape(autocorrs)) + 1.0
        elif autocorr_to_plot in ['dl,ext']:
            autocorrs = ext_delta_DL_cross_products
            autocorr_errs = np.zeros(np.shape(autocorrs)) + 1.0
        elif autocorr_to_plot in ['dl,dl-dl,ext']:
            #First, compute the scaling of extinction that minimizes correlation function of magnitudes (since we assume magnitudes are linearly dependent on extinction)
            min_scaling = np.sum(np.array(delta_mu_product_norms) * np.array(ext_delta_mu_cross_products)) / np.sum(np.array(ext_delta_mu_cross_products) ** 2.0)
            print ('min_scaling = ' + str(min_scaling))
            #Then subtract the extinctions from the distance modulus residuals with this minimum scaling
            ext_corrected_muDiffs = np.array(binned_muDiffs) - np.array(binned_exts) * min_scaling
            ext_corrected_deltaDLs = [10.0 ** (muDiff / 5.0) - 1 for muDiff in ext_corrected_muDiffs]
            autocorrs = c.flattenListOfLists([[ext_corrected_deltaDLs[i] * ext_corrected_deltaDLs[j] for j in range(i+1, len(ext_corrected_deltaDLs)) ] for i in range(len(ext_corrected_deltaDLs))] )
            autocorr_errs = np.zeros(np.shape(autocorrs)) + 1.0
        elif autocorr_to_plot in ['random']:
            autocorrs = np.random.normal(0.0, 1.0, np.shape(ang_seps))
            autocorr_errs = np.zeros(np.shape(autocorrs)) + 1.0

        plt.hist(autocorrs, bins = 101)
        plt.show()

        comoving_to_redshift_interp = interpolate.interp1d([ speedoflight / H0 * integrate.quad(lambda z_int: 1.0 / math.sqrt( OmM * (1.0 + z_int) ** 3.0 + OmL + OmR * (1.0 + z) ** 4.0 ) , 0 , z)[0] for z in np.arange(0.0, 3.0, 0.001) ], np.arange(0.0, 3.0, 0.001) )
        redshift_centers = [comoving_to_redshift_interp(center) for center in comoving_centers]
        comoving_parallels = [0 for center in comoving_centers]
        pair_num = 0
        for i in range(len(comoving_dist)):
            for j in range(i+1, len(comoving_dist)):
                comoving_parallels[pair_num] = 0.5 * abs(comoving_dist[i] ** 2.0 - comoving_dist[j] ** 2.0) / comoving_centers[pair_num]
                pair_num = pair_num + 1

        comoving_perps = [0 for center in comoving_centers]
        pair_num = 0
        for i in range(len(comoving_dist)):
            for j in range(i+1, len(comoving_dist)):
                comoving_perps[pair_num] = np.abs(np.sin(ang_seps[pair_num])) * comoving_dist[i] * comoving_dist[j] / comoving_centers[pair_num]
                pair_num = pair_num + 1

        ps = ( np.array(comoving_parallels) ** 2.0 + np.array(comoving_perps) ** 2.0 ) ** 0.5

        if plot_perp_parallel:
            xs = comoving_perps
            ys = comoving_parallels
        else:
            xs = comoving_seps
            ys = redshift_centers

        max_xs = max(xs)
        min_xs = min(xs)
        max_ys = max(ys)
        min_ys = min(ys)
        min_autocorr = min(autocorrs)
        max_autocorr = max(autocorrs)

        x_lims = x_lims_in_zranges[z_bin_num]
        if x_lims[1] in ['max', 'Max','MAX']:
            x_lims[1] = max_xs
        y_lims = y_lims_in_zranges[z_bin_num]
        if y_lims[1] in ['max', 'Max','MAX']:
            y_lims[1] = max_ys

        if plot_perp_parallel: n_xy_bins = [ n_x_pixels, int(n_x_pixels * (y_lims[1] - y_lims[0]) / (x_lims[1] - x_lims[0])) ]
        else: n_xy_bins = [n_x_pixels, n_x_pixels]

        trimmed_xs = [xs[i] for i in range(len(autocorrs)) if (xs[i] > x_lims[0] and xs[i] < x_lims[1] and ys[i] > y_lims[0] and ys[i] < y_lims[1]) ]
        trimmed_ys = [ys[i] for i in range(len(autocorrs)) if (xs[i] > x_lims[0] and xs[i] < x_lims[1] and ys[i] > y_lims[0] and ys[i] < y_lims[1]) ]
        trimmed_autocorrs = [autocorrs[i] for i in range(len(autocorrs)) if (xs[i] > x_lims[0] and xs[i] < x_lims[1] and ys[i] > y_lims[0] and ys[i] < y_lims[1]) ]
        trimmed_autocorr_errs = [autocorr_errs[i] for i in range(len(autocorrs)) if (xs[i] > x_lims[0] and xs[i] < x_lims[1] and ys[i] > y_lims[0] and ys[i] < y_lims[1]) ]
        xs, ys, autocorrs, autocorr_errs = [trimmed_xs, trimmed_ys, trimmed_autocorrs, trimmed_autocorr_errs]
        if randomize:
            rand_autocorrs, rand_autocorr_errs = rsd.randomShuffleListOfLists([autocorrs, autocorr_errs])
            autocorrs = rand_autocorrs
            autocorr_errs = rand_autocorr_errs

        ref_dof = len(autocorrs) - 1
        null_shift = c.weighted_mean(np.array(autocorrs), autocorr_errs)
        ref_rchisqr = np.sum(((np.array(autocorrs) - null_shift) / np.array(autocorr_errs)) ** 2.0 ) / ref_dof
        print ('ref_rchisqr = ' + str(ref_rchisqr))


        #poly_for_DL_of_ang_sep = np.polyfit(ang_seps, np.array(delta_DL_products) * 10.0 ** 6.0, 2, w = 1.0 / np.array(delta_DL_product_errs) * 10.0 ** 6.0)
        fit_funct = lambda thetas, A, theta0: A / (1.0 + thetas / theta0)
        theta0s = np.arange(0.5, 10.0, 0.5)
        As = np.linspace(-0.5 * 10.0 ** -5.0, 0.5 * 10.0 ** -5.0, 3)
        As = []
        fit_params = [[0.0 for theta0 in theta0s] for A in As]
        fit_chisqrs = [[0.0 for theta0 in theta0s] for A in As]
        for i in range(len(As)):
            A = As[i]
            if verbose: print ('Working on A = ' + str(A))
            for j in range(len(theta0s)):
                dof = len(delta_mu_product_norms) - 3
                theta0 = theta0s[j]
                curve = fit_funct(np.array(ang_seps), A, theta0)
                offset = c.weighted_mean(curve - np.array(autocorrs), autocorr_errs)
                r_chi_sqr = np.sum(((curve - offset) - np.array(autocorrs)) ** 2.0 / (1 ** 2.0)) / dof
                fit_params[i][j] = [A, theta0, offset]
                fit_chisqrs[i][j] = r_chi_sqr
        bin_sizes = [(x_lims[1] - x_lims[0] ) / n_xy_bins[0], (y_lims[1] - y_lims[0] ) / n_xy_bins[1]]
        print ('[n_xy_bins, bin_sizes] = ' +str([n_xy_bins, bin_sizes]))
        #if verbose: print ('comoving_bin_size = ' + str(comoving_bin_size))
        x_bin_edges = (np.arange(x_lims[0], x_lims[1], bin_sizes[0])).tolist()
        y_bin_edges = (np.arange(y_lims[0], y_lims[1], bin_sizes[1])).tolist()
        if y_bin_edges[-1] < max_ys: y_bin_edges = y_bin_edges + [max_ys]
        if x_bin_edges[-1] < max_xs: x_bin_edges = x_bin_edges + [max_xs]
        x_mesh, y_mesh = np.meshgrid([(x_bin_edges[i] + x_bin_edges[i+1])/2.0 for i in range(len(x_bin_edges)-1)],
                                               [(y_bin_edges[i] + y_bin_edges[i+1])/2.0 for i in range(len(y_bin_edges)-1)])
        #print ('[perp_mesh, parallel_mesh] = ' + str([perp_mesh, parallel_mesh]))


        if fitting_funct_str in ['polynomial_decay']:
            fitting_funct = lambda x, y, params: params[0] / (1.0 + (x / params[1] ) ** params[3] + (y / params[2] ) ** params[4] ) #params = [A, r0_perp, r0_parallel, alpha_perp, alpha_parallel]
            param_bounds = [[-1.0 * np.nanstd(autocorrs), 1.0 * np.nanstd(autocorrs)], [(x_bin_edges[1] + x_bin_edges[0]), (x_bin_edges[-1] + x_bin_edges[-2])], [(y_bin_edges[1] + y_bin_edges[0]), (y_bin_edges[-1] + y_bin_edges[-2])], [0.0, 10.0], [0.0, 10.0] ]
            param_functs = [lambda param_index, bound = bound: bound[0] + (bound[1] - bound[0]) * param_index for bound in param_bounds]
            n_fit_params = len(param_functs)
            As = [0.5]
            x0s = np.linspace(0.0, 1.0, 6)
            y0s = np.linspace(0.0, 1.0, 6)
            alpha_perps = [0.1]
            alpha_parallels = [0.1]
            autocorr_product_fits = [ [ [ [ [[] for m in range(len(alpha_parallels)) ] for l in range(len(alpha_perps)) ] for k in range(len(r0_parallels)) ] for j in range(len(x0s)) ] for i in range(len(As)) ]
            autocorr_rchisqrs = [ [ [ [ [0.0 for m in range(len(alpha_parallels)) ] for l in range(len(alpha_perps)) ] for k in range(len(r0_parallels)) ] for j in range(len(y0s)) ] for i in range(len(As)) ]
            dof = len(autocorrs) - (n_fit_params + 1)
            if verbose: print ('Taking coarse fits...')
            for i in range(len(As)):
                A = As[i]
                #print ('A = ' + str(A) )
                for j in range(len(x0s)):
                    x0 = x0s[j]
                    #print ('muperp = ' + str(muperp) )
                    for k in range(len(y0s)):
                        y0 = y0s[k]
                        #print ('muparallel= ' + str(muparallel) )
                        for l in range(len(alpha_perps)):
                            alpha_perp = alpha_perps[l]
                            for m in range(len(alpha_parallels)):
                                alpha_parallel = alpha_parallels[m]
                                param_scalings = [A, x0, y0, alpha_perp, alpha_parallel]
                                rchisqr, best_fit_param_scalings = functToMinimize(fitting_funct, autocorrs, autocorr_errs, xs, ys, dof, param_scalings, param_functs, verbose = verbose, ref_rchisqr = ref_rchisqr)
                                rchisqr = rchisqr * ref_rchisqr

                                autocorr_product_fits[i][j][k][l][m] = best_fit_param_scalings
                                autocorr_rchisqrs[i][j][k][l][m] = rchisqr

                    if verbose: print ('Taking coarse fits: ' + str(c.round_to_n((i * len(x0s) + j+1) / (len(As) * len(x0s)) * 100, 5)) + '% done.')
            flattened_rchisqr = c.flattenListOfLists(c.flattenListOfLists(c.flattenListOfLists(c.flattenListOfLists(autocorr_rchisqrs))))
            flattened_fits = c.flattenListOfLists(c.flattenListOfLists(c.flattenListOfLists(c.flattenListOfLists(autocorr_product_fits))))
            min_rchisqr_argmin = np.argmin(flattened_rchisqr)
            coarse_min_rchisqr = flattened_rchisqr[min_rchisqr_argmin]
            #coarse_min_param_scalings = flattened_fits[min_rchisqr_argmin][0:-1]
            coarse_min_param_scalings = flattened_fits[min_rchisqr_argmin][0:-1]
            if verbose: print ('coarse_min_param_scalings = ' + str(coarse_min_param_scalings ))
            #Use the above as a coarse guide.  Now do python numerical minimization. That requires that I scale these params so that we're minimizing around unity
            if verbose: print ('Starting minimization...')
            min_fit = optimize.minimize(lambda test_param_scalings: functToMinimize(fitting_funct, autocorrs, autocorr_errs, xs, ys, dof, test_param_scalings, param_functs, verbose = verbose, ref_rchisqr = ref_rchisqr)[0],  coarse_min_param_scalings, bounds = [[0.0, 1.0] for scaling in coarse_min_param_scalings])
            if verbose: print ('Finished minimization...')
            min_param_scalings = min_fit['x']
            min_rchisqr, min_param_scalings = functToMinimize(fitting_funct, autocorrs, autocorr_errs, xs, ys, dof, min_param_scalings, param_functs, verbose = verbose, ref_rchisqr = ref_rchisqr)
            min_param_scalings = min_param_scalings.tolist()
            min_rchisqr = min_rchisqr * ref_rchisqr
            #min_shift = min_param_scalings[-1]
            min_shift = min_param_scalings[-1]
            #min_params = [param_functs[j](min_param_scalings[j]) for j in range(len(param_functs)) ]
            min_params = [param_functs[j](min_param_scalings[j]) for j in range(len(param_functs)) ]
        elif fitting_funct_str in ['gauss','Gaussian','gaussian','gauss']:
            fitting_funct = lambda x, y, params: param[0] * np.exp(- ((x - params[1]) ** 2.0 / (2.0 * params[3] ** 2.0)) - ((y - params[2]) ** 2.0 / (2.0 * (params[3] * params[4]) ** 2.0))) # A, mu_perp, mu_parallel, sig_perp, el
            param_bounds = [[-0.2 * np.nanstd(autocorrs), 0.2 * np.nanstd(autocorrs)], [(x_bin_edges[1] + x_bin_edges[0]) / 2.0, (x_bin_edges[-1] + x_bin_edges[-2]) / 2.0 / 2]
                           [(y_bin_edges[1] + y_bin_edges[0]) / 2.0, (y_bin_edges[-1] + y_bin_edges[-2]) / 2.0 / 2], [20.0, 1000.0], [-8.0, 10.0]] # A, muperp, muparallel, sig, el
            param_functs = [lambda param_index, bound = bound: bound[0] + (bound[1] - bound[0]) * param_index for bound in param_bounds]
            n_fit_params = len(param_functs)
            As = [0.25, 0.5, 0.75] #np.linspace(0.3, 0.0.0, 1.0, 6)
            if verbose: print ('As = ' + str(As))
            mu_xs = np.linspace(0.0, 1.0, 21)
            mu_ys = np.linspace(0.0, 1.0, 21)
            sigs = [0.5] #np.linspace(0.0, 1.0, 6)
            el_proxies = [0.5] #np.linspace(0.0, 1.0, 3)
            autocorr_product_fits = [ [ [ [ [[] for m in range(len(el_proxies)) ] for l in range(len(sigs)) ] for k in range(len(mu_xs)) ] for j in range(len(mu_ys)) ] for i in range(len(As)) ]
            autocorr_rchisqrs = [ [ [ [ [0.0 for m in range(len(el_proxies)) ] for l in range(len(sigs)) ] for k in range(len(mu_xs)) ] for j in range(len(mu_ys)) ] for i in range(len(As)) ]
            dof = len(autocorrs) - (n_fit_params + 1)
            if verbose: print ('Taking coarse fits...')
            for i in range(len(As)):
                A = As[i]
                for j in range(len(mu_xs)):
                    mu_x = mu_xs[j]
                    for k in range(len(mu_ys)):
                        mu_y = mu_ys[k]
                        for l in range(len(alpha_perps)):
                            sigma_perp =sigma_perps[l]
                            for m in range(len(el_proxies)):
                                el_proxy = el_proxies[m]
                                param_scalings = [A, mu_x, mu_y, sig_perp, el_proxy]
                                rchisqr, best_fit_param_scalings = functToMinimize(fitting_funct, autocorrs, autocorr_errs, xs, ys, dof, param_scalings, param_functs, verbose = verbose, ref_rchisqr = ref_rchisqr)
                                rchisqr = rchisqr * ref_rchisqr
                                autocorr_product_fits[i][j][k][l][m] = best_fit_param_scalings
                                autocorr_rchisqrs[i][j][k][l][m] = rchisqr
            flattened_rchisqr = c.flattenListOfLists(c.flattenListOfLists(c.flattenListOfLists(c.flattenListOfLists(autocorr_rchisqrs))))
            flattened_fits = c.flattenListOfLists(c.flattenListOfLists(c.flattenListOfLists(c.flattenListOfLists(autocorr_product_fits))))
            min_rchisqr_argmin = np.argmin(flattened_rchisqr)
            coarse_min_rchisqr = flattened_rchisqr[min_rchisqr_argmin]
            #coarse_min_param_scalings = flattened_fits[min_rchisqr_argmin][0:-1]
            coarse_min_param_scalings = flattened_fits[min_rchisqr_argmin][0:-1]
            if verbose: print ('coarse_min_param_scalings = ' + str(coarse_min_param_scalings ))
            #Use the above as a coarse guide.  Now do python numerical minimization. That requires that I scale these params so that we're minimizing around unity
            if verbose: print ('Starting minimization...')
            min_fit = optimize.minimize(lambda test_param_scalings: functToMinimize(fitting_funct, autocorrs, autocorr_errs, xs, ys, dof, test_param_scalings, param_functs, verbose = verbose, ref_rchisqr = ref_rchisqr)[0],  coarse_min_param_scalings, bounds = [[0.0, 1.0] for scaling in coarse_min_param_scalings])
            if verbose: print ('Ending minimization...')
            min_param_scalings = min_fit['x']
            min_rchisqr, min_param_scalings = functToMinimize(fitting_funct, autocorrs, autocorr_errs, xs, ys, dof, min_param_scalings, param_functs, verbose = verbose, ref_rchisqr = ref_rchisqr)
            min_param_scalings = min_param_scalings.tolist()
            min_rchisqr = min_rchisqr * ref_rchisqr
            #min_shift = min_param_scalings[-1]
            min_shift = min_param_scalings[-1]
            #min_params = [param_functs[j](min_param_scalings[j]) for j in range(len(param_functs)) ]
            min_params = [param_functs[j](min_param_scalings[j]) for j in range(len(param_functs)) ]
        elif fitting_funct_str[0:-1] in ['inv_poly', 'inverse_polynomial']:
            order = int(fitting_funct_str[-1])
            n_coefs = sum([order + 1 for order in range(order+1)]) - 1
            dof = len(autocorrs) - (n_coefs + 1)
            term_powers = [[i % (order+1),  i // (order+1)] for i in range((order + 1) ** 2 )  if (i % (order+1)+ i // (order+1) <=order and i % (order+1)+ i // (order+1) > 0)]
            inv_autocorrs = 1.0 / np.array(autocorrs)
            inv_autocorr_errs = autocorr_errs / (np.array(autocorrs) ** 2.0)
            #fit_res = c.polyfit2d(comoving_perps, comoving_parallels, inv_autocorrs / inv_autocorr_errs, kx = order, ky = order, order = order)[0]
            #fitting_funct = lambda r_perp, r_parallel, poly_terms: ((1.0 / np.sum([poly_terms[i] * r_perp ** (term_powers[i][0]) * r_parallel ** (term_powers[i][1]) for i in range(n_coefs) ]) - autocorrs) / (autocorr_errs)) ** 2.0
            middle_x = (max_xs + min_xs) / 2.0
            middle_y = (max_ys + min_ys) / 2.0
            fitting_funct = lambda x, y, poly_terms, term_powers = term_powers, middle_x = middle_x, middle_y = middle_y: poly_terms[-1] / ( 1.0 + np.sum([poly_terms[i] * (x / middle_x) ** (term_powers[i][0]) * (y / middle_y) ** (term_powers[i][1]) for i in range(n_coefs) ], axis = 0))
            #param_bounds = [[-1.0 * np.nanstd(autocorrs), 1.0 * np.nanstd(autocorrs)] for coef in range(n_coefs)] + [[-1.0 * np.nanstd(autocorrs), 1.0 * np.nanstd(autocorrs)]]
            param_bounds = [[-5.0, 5.0] for coef in range(n_coefs)] + [[-1.0 * np.nanstd(autocorrs), 1.0 * np.nanstd(autocorrs)]]
            param_functs = [lambda param_index, bound = bound: bound[0] + (bound[1] - bound[0]) * param_index for bound in param_bounds]

            if verbose: print ('Starting minimization...')
            fit_res = optimize.minimize(lambda test_param_scalings: functToMinimize(fitting_funct, autocorrs, autocorr_errs, xs, ys, dof, test_param_scalings, param_functs, verbose = verbose, ref_rchisqr = ref_rchisqr)[0],  [0.5 for coef in range(n_coefs)] + [0.51] , bounds = [[0.0, 1.0] for scaling in param_bounds])
            if verbose: print ('Ending minimization...')

            if verbose: print ('fit_res = ' + str(fit_res))
            #rchisqr = np.sum(((fitting_funct(comoving_perps, comoving_parallels, fit_res) - autocorrs) / (autocorr_errs ** 1.0)) ** 2.0 )['fun']
            min_param_scalings = fit_res['x'].tolist()
            min_rchisqr, min_param_scalings = functToMinimize(fitting_funct, autocorrs, autocorr_errs, xs, ys, dof, min_param_scalings, param_functs, verbose = verbose, ref_rchisqr = ref_rchisqr)
            min_param_scalings = min_param_scalings.tolist()
            #print ('[min_rchisqr, min_param_scalings] = ' +str([min_rchisqr, min_param_scalings]))
            min_rchisqr = min_rchisqr * ref_rchisqr
            #min_shift = min_param_scalings[-1]
            min_shift = min_param_scalings[-1]
            #min_params = [param_functs[j](min_param_scalings[j]) for j in range(len(param_functs)) ]
            min_params = [param_functs[j](min_param_scalings[j]) for j in range(len(param_functs)) ]
            print ('[min_rchisqr, min_params, min_shift] = ' + str([min_rchisqr, min_params, min_shift]))

        else : #Just compute the weighted shift from a constant
            if verbose: print ('Just fitting a constant offset to correlations.')
            fitting_funct = lambda x, y, A: x * 0.0
            param_bounds = [ [0.0, 1.0] ]
            param_scalings = [0.5]
            autocorr_product_fits =  [[] for m in range(len(param_scalings)) ]
            autocorr_rchisqrs = [0.0 for m in range(len(param_scalings)) ]
            param_functs = [lambda param_index, bound = bound: bound[0] + (bound[1] - bound[0]) * param_index for bound in param_bounds]
            dof = len(autocorrs) - 1
            rchisqr, best_fit_param_scalings = functToMinimize(fitting_funct, autocorrs, autocorr_errs, xs, ys, dof, param_scalings, param_functs, verbose = verbose, ref_rchisqr = ref_rchisqr)
            rchisqr = rchisqr * ref_rchisqr
            autocorr_product_fits[0] = best_fit_param_scalings
            autocorr_rchisqrs[0] = rchisqr
            min_rchisqr = autocorr_rchisqrs[0]
            min_fit = autocorr_product_fits[0]
            min_shift = min_fit[-1]


        #min_params = [0.00002, 1000.0, 100.0, 500, 1.0, 0.0]
        if verbose: print ('[min_rchisqr, min_param_scalings, min_params, min_shift] = ' + str([min_rchisqr, min_param_scalings, min_params, min_shift]))

        fitted_autocorrs = fitting_funct(x_mesh, y_mesh, (min_params)) - min_shift
        #fitted_ext_product = fitting_funct(perp_mesh, parallel_mesh, *(min_ext_params)) - min_ext_shift

        if (show_fig or save_fig) :
            if verbose: print ('I am binning/smoothing the measured covariances...')
            if verbose: print ('smoothing = ' + str(smoothing))
                #Bin in such a way that points not in any bin have an index of -1
            xy_binning = [[sum([j if (xs[k] >= x_bin_edges[j-1] and xs[k] < x_bin_edges[j]) else 0 for j in range(1, len(x_bin_edges)) ]) - 1,
                                  sum([i if (ys[k] >= y_bin_edges[i-1] and ys[k] < y_bin_edges[i]) else 0 for i in range(1, len(y_bin_edges)) ]) - 1]
                              for k in range(len(xs))]
            #print ('xy_binning = ' + str(xy_binning))
            autocorrs_inBins = [[[] for j in range(len(x_bin_edges)-1) ] for i in range(len(y_bin_edges)-1)]
            autocorr_errs_inBins = [[[] for j in range(len(x_bin_edges)-1) ] for i in range(len(y_bin_edges)-1)]

            for k in range(len(xy_binning)):
                if (k % 10000 == 0 and verbose): print (str(c.round_to_n(k / len(xy_binning) * 100, 3)) + '% done.')
                x_bin, y_bin = xy_binning[k]
                if x_bin >= 0 and y_bin >= 0:
                    autocorr = autocorrs[k]
                    autocorr_err = autocorr_errs[k]
                    autocorrs_inBins[y_bin][x_bin] = autocorrs_inBins[y_bin][x_bin] + [autocorr]
                    autocorr_errs_inBins[y_bin][x_bin] = autocorr_errs_inBins[y_bin][x_bin] + [autocorr_err]

            binned_autocorrs = [[c.weighted_mean(autocorrs_inBins[i][j], autocorr_errs_inBins[i][j]) for j in range(len(x_bin_edges)-1) ] for i in range(len(y_bin_edges)-1)]
            binnedAutocorr_errs = [[math.sqrt(len(autocorrs_inBins[i][j]) / sum([err ** -2.0 for err in autocorrs_inBins[i][j]])) /np.sqrt(len(autocorrs_inBins[i][j]))  if len(autocorrs_inBins[i][j]) >=1 else 0.0 for j in range(len(x_bin_edges)-1) ] for i in range(len(y_bin_edges)-1)] #standard error calculation
            if smoothing in ['gauss', 'gaussian', 'Gauss', 'Gaussian', 'GAUSS', 'GAUSSIAN']:
                smooth_funct = lambda x, y, mu_x, mu_y, sig_x, sig_y: 1.0 / np.sqrt(np.pi * sig_x ** 2.0 * sig_y ** 2.0) * np.exp(-((x - mu_x) ** 2.0 / (2.0 * sig_x ** 2.0) + (y - mu_y) ** 2.0 / (2.0 * sig_y ** 2.0)  ) )
                smoothing_widths = [bin_sizes[0], bin_sizes[1] ]
                smoothed_autocorrs = [ [0.0 for j in range(len(x_bin_edges)-1) ] for i in range(len(y_bin_edges)-1) ]
                smoothedAutocorr_errs = [ [0.0 for j in range(len(x_bin_edges)-1) ] for i in range(len(y_bin_edges)-1) ]
                for i in range(len(y_bin_edges)-1):
                    if verbose: print ('Gaussian smoothing is ' + str(c.round_to_n(i / (len(y_bin_edges)-1), 3) * 100) + '% done.' )
                    for j in range(len(x_bin_edges)-1):
                        #print ('Binning [i,j] = ' + str([i,j]))
                        x = (x_bin_edges[j+1] + x_bin_edges[j]) / 2.0
                        y = (y_bin_edges[i+1] + y_bin_edges[i]) / 2.0
                        correlation_weights = smooth_funct(np.array(xs), np.array(ys), x, y, *smoothing_widths)
                        smoothed_autocorrs[i][j] = np.sum(np.array(correlation_weights) * np.array(autocorrs)) / np.sum(correlation_weights)
                        smoothedAutocorr_errs[i][j] = np.sqrt(np.sum((np.array(correlation_weights) * np.array(autocorr_errs)) ** 2.0))
            else:
                smoothed_autocorrs = binned_autocorrs
                smoothedAutocor_errs = binnedautocorr_errs
            if verbose: print ('I am done binning/smoothing the measured covariances.')
            smoothedAutocorrDeviations = np.array(smoothed_autocorrs) / np.array(smoothedAutocorr_errs)

            n_levels = 21
            bin_display_bounds = [np.nanmean(binned_autocorrs) - n_std_for_levels[z_bin_num] * np.nanstd(binned_autocorrs), np.nanmean(binned_autocorrs) + n_std_for_levels[z_bin_num] * np.nanstd(binned_autocorrs)]
            bin_display_bounds = [-0.0005, 0.0025] #Force it when you need to compare two plots with the same scaling

            gs_this_z_bin = gridspec.GridSpecFromSubplotSpec(3, 5, subplot_spec=gs_z_bins[z_bin_num])

            ticksize = 6
            matplotlib.rc('xtick', labelsize=ticksize)
            matplotlib.rc('ytick', labelsize=ticksize)
            ax1 = f.add_subplot(gs_this_z_bin[0, 0:4])
            ax1.set_title(str(z_bin_range[0]) + r'$< z <$' + str(z_bin_range[1]))
            ax2 = f.add_subplot(gs_this_z_bin[1, 0:4])
            ax3 = f.add_subplot(gs_this_z_bin[2, 0:4])
            ax_contour_cbar = f.add_subplot(gs_this_z_bin[0:2, 4])
            ax_fit_cbar = f.add_subplot(gs_this_z_bin[2, 4])
            im1 = ax1.imshow(binned_autocorrs, vmin = bin_display_bounds[0], vmax = bin_display_bounds[1] )
            #im1 = axarr[0,0].imshow(smoothedAutocorrDeviations, vmin = np.nanmean(smoothedAutocorr) - 0.5 * np.nanstd(smoothedAutocorr), vmax = np.nanmean(smoothedAutocorr) + 0.5 * np.nanstd(smoothedAutocorr) )
            #ax1.set_xlabel(r'$r_{\perp}$ [Mpc]')
            if plot_perp_parallel: ax1.set_ylabel(r'$r_{\parallel}$ [Mpc/h]')
            else: ax1.set_ylabel(r'z')
            ax1.set_xticks(np.linspace(x_lims[0]/bin_sizes[0], x_lims[1]/bin_sizes[0], 9) )
            ax1.set_yticks(np.linspace(y_lims[0]/bin_sizes[1], y_lims[1]/bin_sizes[1], 9) )
            ax1.set_xticklabels([c.round_to_n(tick* bin_sizes[0] * h, 3)  for tick in ax1.get_xticks()])
            print ('[ax1.get_yticks, bin_sizes[1], h] = ' + str([ax1.get_yticks, bin_sizes[1], h]))
            ax1.set_yticklabels( [c.round_to_n(tick * bin_sizes[1] * h, 3)  for tick in ax1.get_yticks()] )
            n_contours = 11
            ax2.imshow(smoothed_autocorrs, vmin = bin_display_bounds[0], vmax = bin_display_bounds[1])
            #ax2.set_xlabel(r'$r_{\perp}$ [Mpc]')
            if plot_perp_parallel: ax2.set_ylabel(r'$r_{\parallel}$ [Mpc/h]')
            else: ax2.set_ylabel(r'z')
            ax2.set_xticks(np.linspace(x_lims[0]/bin_sizes[0], x_lims[1]/bin_sizes[0], 6) )
            ax2.set_yticks(np.linspace(y_lims[0]/bin_sizes[1], y_lims[1]/bin_sizes[1], 6) )
            ax2.set_xticklabels([c.round_to_n(tick * bin_sizes[0] * h, 3)  for tick in ax2.get_xticks()])
            ax2.set_yticklabels( [c.round_to_n(tick * bin_sizes[1] * h, 3)  for tick in ax2.get_yticks()] )
            fit_levels = sorted(list(set(np.linspace(np.min(fitted_autocorrs), np.max(fitted_autocorrs), 11))))
            if len(fit_levels) == 1: fit_levels = [fit_levels[0], fit_levels[0] + 1]
            #print ('fit_levels = ' + str(fit_levels))
            #im2 = ax3.contourf(perp_mesh, parallel_mesh, fitted_autocorrs, levels = fit_levels) #levels = np.linspace(np.nanmean(smoothedAutocorr) - 0.5 * np.nanstd(smoothedAutocorr), np.nanmean(smoothedAutocorr) + 0.5 * np.nanstd(smoothedAutocorr), 21))
            im2 = ax3.imshow(fitted_autocorrs, vmin = fit_levels[0], vmax = fit_levels[-1])
            if plot_perp_parallel:
                ax3.set_xlabel(r'$r_{\perp}$ [Mpc/h]')
                ax3.set_ylabel(r'$r_{\parallel}$ [Mpc/h]')
            else:
                ax3.set_xlabel(r'$r_{sep}$ [Mpc/h]')
                ax3.set_ylabel(r'$z$')
            ax3.set_xticks(np.linspace(x_lims[0]/bin_sizes[0], x_lims[1]/bin_sizes[0], 6) )
            ax3.set_yticks(np.linspace(y_lims[0]/bin_sizes[1], y_lims[1]/bin_sizes[1], 6) )
            ax3.set_xticklabels([c.round_to_n(tick * bin_sizes[0] * h, 3)  for tick in ax3.get_xticks()])
            ax3.set_yticklabels( [c.round_to_n(tick * bin_sizes[1] * h, 3)  for tick in ax3.get_yticks()] )
            cbar_contour = f.colorbar(im1, cax=ax_contour_cbar)
            cbar_fit = f.colorbar(im2, cax=ax_fit_cbar)
            cbar_contour.set_label(colorbar_label, rotation=270, labelpad=20)
            #gs_this_z_bin.tight_layout()
            gs_this_z_bin
            gs_z_bins.tight_layout(f)
        print ('min_params = ' + str(min_params))

        min_rchisqrs[z_bin_num] = min_rchisqr
        min_param_scalings_sets[z_bin_num] = min_param_scalings
        min_params_sets[z_bin_num] = min_params[:]
        min_shifts[z_bin_num] = min_shift

    if save_fig:
        print ('Saving figure to file ' + plt_save_dir + fig_name)
        plt.savefig(plt_save_dir + fig_name)
    if show_fig:
        plt.show()

    if save_fit:
        if not(os.path.exists(res_save_dir + fit_results_file_name)):
            fit_file = open(res_save_dir + fit_results_file_name, 'w')
            #fit_file.write('param_bounds = [[-0.2 * np.nanstd(delta_mu_product_norms), 0.2 * np.nanstd(delta_mu_product_norms)], [(perp_bin_edges[1] + perp_bin_edges[0]) / 2.0, (perp_bin_edges[-1] + perp_bin_edges[-2]) / 2.0 / 2], [(parallel_bin_edges[1] + parallel_bin_edges[0]) / 2.0, (parallel_bin_edges[-1] + parallel_bin_edges[-2]) / 2.0 / 2], [20.0, 1000.0], [-8.0, 10.0]]' + os.linesep)
            #fit_file.write('As = [0.15, 0.35, 0.5, 0.65, 0.85], muperps = np.linspace(0.0, 1.0, 21), muparallels = np.linspace(0.0, 1.0, 21), sigs = [0.5], el_proxies = [1.0] '+ os.linesep)
        else:
            fit_file = open(res_save_dir + fit_results_file_name, 'a')
        print ('[min_rchisqrs, min_param_scalings_sets, min_params_sets, min_shifts] = ' + str([min_rchisqrs, min_param_scalings_sets, min_params_sets, min_shifts]))
        fit_file.write(('random : ' if randomize else 'real: ') + '[min_rchisqrs, min_param_scalings, min_params, min_shift] = ' + str([min_rchisqrs, min_param_scalings_sets, min_params_sets, min_shifts])  + os.linesep)
        fit_file.close()

    end = time.time()
    print ('Full analysis took ' + str(end - start) + 's')

    print ('min_rchisqrs = ' + str(min_rchisqrs))
    return [min_rchisqrs, min_params_sets]

def superFuncToMinimize(params, z_bin, n_z_bins, fitting_funct, extinction_correction_fit_funct):

    rchisqr = runMinimization(fitting_funct_str = 'constant', show_fig = 0, save_fig = 0, verbose = 0, n_z_bins = n_z_bins,
                              extinction_correction_functs = [lambda ext: extinction_correction_fit_funct(ext, *params) for j in range(n_z_bins)] )[0][z_bin]
    print ('[params, rchisqr] = ' + str([params, rchisqr]))
    return rchisqr


def runSuperMinimizationOnExtinction():
    all_sn = loadSN(1, ['all'], pull_extinctions = 1)
    all_sn = [sn for sn in all_sn if not(np.isnan(sn['extinction'])) ]
    sn_ext = [sn['extinction'] for sn in all_sn]
    sn_ext_errs = [sn['extinctionErr'] for sn in all_sn]
    mean_ext = c.weighted_mean(sn_ext, sn_ext_errs)
    print ('weighted mean extinction = ' +str(mean_ext))
    n_z_bins = 2
    #benchmark rchisqr = [rchisqr for <dL, dL> autocorrelations after removing weighted mean] =

    fitting_funct = 'constant'
    #extinction_correction_fit_funct = lambda ext, E0, E1, E2:  ((ext - mean_ext) / mean_ext) ** 2.0 * E2 + (ext - mean_ext) / mean_ext * E1 + E0 * mean_ext
    extinction_correction_fit_funct = lambda ext, E0, E1:   (ext - mean_ext) / mean_ext * E1 + E0 * mean_ext

    minimizations = [[] for n_z_bin in range(n_z_bins)]
    for z_bin in range(n_z_bins):
        minimizations[z_bin] = optimize.minimize( lambda params: superFuncToMinimize(params, z_bin, n_z_bins, fitting_funct, extinction_correction_fit_funct), [0.0, 0.0])
        print ('!!!!!!!!Finished one minimization!!!!!!!!')
    print ('minimizations = ' + str(minimizations))

    autocorr_to_plot = 'dl'
    fig_name = 'DeltaDLCorrelations_LinearExtinctionCorrection' + '_gaussSmoothing_goodExtinctionsPantheonDataReal_TwoZBins_gaussfitA.png'
    runMinimization(fitting_funct_str = 'polynomial_decay', show_fig = 0, save_fig = 1, verbose = 1, n_z_bins = n_z_bins,
                    extinction_correction_functs = [lambda ext: extinction_correction_fit_funct(ext, *(minimization['x'])) for minimization in minimizations],
                    autocorr_to_plot = autocorr_to_plot, fig_name = fig_name)

    plt.close('all')



if __name__ == "__main__":
    sys_args = sys.argv[1:]
    if len(sys_args) > 0:
        fig_name = sys_args[0]
    else:
        fig_name = None
    print ('fig_name = ' + str(fig_name))
    surveys_to_load = ['all']
    surveys_to_excise = []
    zHD = 1
    n_z_bins = 1
    only_good_ext = 0
    all_sn = loadSN(1, surveys_to_load, pull_extinctions = 0, zHD = zHD)
    if only_good_ext:
        all_sn = [sn for sn in all_sn if not(np.isnan(sn['extinction'])) ]
    all_sn = [sn for sn in all_sn if not(sn['survey'] in surveys_to_excise)]
    all_zs = [sn['z'] for sn in all_sn]
    sorted_zs, sorted_sn = c.safeSortOneListByAnother(all_zs, [all_zs, all_sn])

    z_bin_ranges = [[sorted_zs[i * (len(sorted_zs) // n_z_bins)], sorted_zs[(i+1) * (len(sorted_zs) // n_z_bins) -1]]  for i in range(n_z_bins)]
    print ('z_bin_ranges = ' + str(z_bin_ranges))
    #z_bin_ranges = [[0.0, 0.2], [0.2, max(sorted_zs)]]
    #print ('z_bin_ranges = ' + str(z_bin_ranges))
    print (runMinimization(fig_name = fig_name, fitting_funct_str = 'inv_poly1', show_fig = 1, save_fig = 0, save_fit = 0, verbose = 1, autocorr_to_plot = 'dl', n_z_bins = 2, surveys_to_load = surveys_to_load, plot_perp_parallel = 1, zHD = zHD, z_bin_ranges = z_bin_ranges, surveys_to_excise = surveys_to_excise ))

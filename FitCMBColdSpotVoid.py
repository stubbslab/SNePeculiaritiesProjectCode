import cantrips as can
import matplotlib.pyplot as plt
import makePlotOfPS1MDFieldsClass as mpfc
import numpy as np
import time

if __name__ == "__main__":
    z_tests = [0.155] # [0.35] #[0.155]
    frac_comoving_dists_to_bin = 2.5 / 5
    central_coords = [[48.29990842, -20.43730437]] # [[35.0, -5.0]] # [[48, -19]] #RA and Dec, in degrees
    cutoff_radius_Mpc = np.inf
    ticksize = 8
    labelsize = 10

    do_randomization_by_field = 0
    do_randomization_by_survey = 0
    randomize_all_sn = 0
    gal_dens_weighting = 0.0 #can be 0.5
    save_plot = 1
    z_range = [-0.1, 3.0]
    overdensity_param = 200 #Delta, the fraction of background mass density that the average halo mass density must be
    n_good_fits_to_show = 5
    resid_profile_funct = 'exp_void'
    bounds = [[-20.0, 0.0], [0.5, 2.5]]
    n_grid_samples = [101, 81] #[121, 41]# [101, 61]
    sn_data_type = 'pantheon_plus' #'pantheon_plus' #real - for old Pantheon data
    zHD = 1
    init_guess = [0.0, 2]

    sdssdir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/SNIsotropyProject/SDSSGalaxies/'
    fulldata_fastRead = can.readInColumnsToList('SDSS_fullCoverage_SDSSGals_pzAll.csv', sdssdir, delimiter = ',', n_ignore = 2, all_np_readable = 1)
    field_plotter = mpfc.PanStarsFieldManager(1, full_sdss_gal_data_file = 'SDSS_fullCoverage_SDSSGals_pzAll.csv', preloaded_sdss_gals = fulldata_fastRead, gal_dens_weighting = gal_dens_weighting, z_range = z_range, sn_data_type = sn_data_type, zHD = zHD, cutoff_radius_Mpc = cutoff_radius_Mpc, NFW_overdensity_param = overdensity_param, resid_profile_funct = resid_profile_funct, randomize_all_sn = randomize_all_sn)

    comoving_bin = field_plotter.r_of_z_interp(min(z_tests)) * frac_comoving_dists_to_bin
    comoving_bin = 250
    print ('comoving_bin = ' + str(comoving_bin))
    min_n_sn = 1
    points_to_test = list(can.cartesian([z_tests, list(range(len(central_coords)))]))
    for i in range(len(points_to_test)) :
         print ('points_to_test[i] = ' + str(points_to_test[i]))
         points_to_test[i] = [points_to_test[i][0]] + central_coords[int(points_to_test[i][1])]
    print ('points_to_test = ' + str(points_to_test))

    fits = [ field_plotter.doSingleFitAroundCoord(point[0], point[1:], comoving_bin, min_n_sn, one_d_fit = 0, init_guess = init_guess, bounds = bounds, n_grid_samples = n_grid_samples, method = 'grid', print_params = 0, show_coarse_fit = 0) for point in points_to_test ]
    n_randomizations = 300
    rand_fit_ratios = [[] for point in points_to_test]
    rand_field_plotter = mpfc.PanStarsFieldManager(1, full_sdss_gal_data_file = 'SDSS_fullCoverage_SDSSGals_pzAll.csv', preloaded_sdss_gals = fulldata_fastRead, gal_dens_weighting = gal_dens_weighting, z_range = z_range, sn_data_type = sn_data_type, zHD = zHD, cutoff_radius_Mpc = cutoff_radius_Mpc, NFW_overdensity_param = overdensity_param, resid_profile_funct = resid_profile_funct, randomize_all_sn = 1)

    for i in range(len(points_to_test)):
        point = points_to_test[i]
        start_time = time.time()
        prev_time = start_time
        for j in range(n_randomizations):
            rand_field_plotter.initializeSN(rand_field_plotter.z_range, ['all'], [], randomize_all = 1)
            rand_field_plotter.all_muResids = [sn['mu'] - rand_field_plotter.mu_of_z(sn['z']) for sn in rand_field_plotter.all_sns]
            rand_fit = rand_field_plotter.doSingleFitAroundCoord(point[0], point[1:], comoving_bin, min_n_sn, one_d_fit = 0, init_guess = init_guess, bounds = bounds, n_grid_samples = n_grid_samples, method = 'grid', print_params = 0, show_coarse_fit = 0)
            rand_fit_ratios[i] = rand_fit_ratios[i] + [rand_fit[1] / rand_fit[2]]
            curr_time = time.time()
            print ('After [i, j] = ' + str([i,j]) + ', expecting ' + str(can.round_to_n((curr_time - prev_time) * (n_randomizations - j), 5)) + 's more for i = ' + str(i))
            prev_time = curr_time

    print ('fits = ' + str(fits))
    if n_randomizations > 0:
        real_fit_ratio = fits[0][1] / fits[0][2]
        print ('real fit ratio = ' + str(real_fit_ratio))
        print ('rand_fit_ratios[0] = ' + str(rand_fit_ratios[0]))
        n_better = len([ratio for ratio in rand_fit_ratios[0] if ratio < real_fit_ratio])
        hist = plt.hist(rand_fit_ratios, bins = max(10, n_randomizations // 5), edgecolor = 'k', color = 'white')
        hist_levels = hist[0]
        hist_centers= hist[1]
        hist_max = np.max(hist_levels)
        max_center = np.max(hist_centers)
        plt.axvline(real_fit_ratio, linestyle = '--', c = 'r')
        plt.text(real_fit_ratio, hist_max / 2, 'Real improvement over null', rotation = 90, c = 'r', horizontalalignment = 'right', verticalalignment = 'center')
        plt.text(max_center, hist_max * 9 / 10, 'Real data > ' + str(len(rand_fit_ratios[0]) - n_better) + '/' + str(len(rand_fit_ratios[0])) + ' of bootstrapped data', c = 'k', horizontalalignment = 'right', verticalalignment = 'center')
        plt.xlabel('Best fit vs null ' + r'$\chi^2_\nu$' + ' ratio')
        plt.ylabel('N randomizations in bin')
        plt.title('Permutation test Void fits to CMB coldspot')
        plt.savefig('/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/SNIsotropyProject/plots/' + 'CMBColdSpot_HistOfPermutationTests2.pdf')
        plt.close('all')
    save_prefix = 'ColdSpot_supervoid'

    f, axarr  = plt.subplots(2, len(z_tests), squeeze = False, figsize = (5 * len(z_tests), 8))
    for i in range(len(z_tests)):
        fitted_unitless_density_param, fitted_radius_power = fits[i][0]
        #fitted_unitless_density_param, fitted_radius_power = [-20.0, 2]
        central_coord = [z_tests[i]] + central_coords[i]
        sn_indeces_to_plot = field_plotter.getSNWithinComovingDistance(z_tests[i], central_coords[i], comoving_bin)
        #print ('sn_indeces_to_plot = ' + str(sn_indeces_to_plot))
        zs_to_plot = [field_plotter.all_zs[index] for index in sn_indeces_to_plot]
        center_line_zs = np.linspace(min(zs_to_plot), max(zs_to_plot), 101)
        muResids_to_plot = [field_plotter.all_muResids[index] for index in sn_indeces_to_plot]
        mu_errs_to_plot = [field_plotter.all_mu_errs[index] for index in sn_indeces_to_plot]
        RAs_to_plot = [field_plotter.all_RAs[index] for index in sn_indeces_to_plot]
        Decs_to_plot = [field_plotter.all_Decs[index] for index in sn_indeces_to_plot]
        surveys_to_plot = [field_plotter.all_surveys[index] for index in sn_indeces_to_plot]
        unique_surveys_to_plot = np.unique(surveys_to_plot).tolist()

        colors_to_plot = [field_plotter.survey_to_color_dict[survey] for survey in surveys_to_plot]
        RA_offset_to_center = (180.0 - central_coord[1])
        RAs_for_resids = field_plotter.centerRAs(RAs_to_plot, RA_offset_to_center)
        central_RA, central_Dec = [field_plotter.centerRAs([central_coord[1]],  RA_offset_to_center )[0],  central_coord[2]]
        fitted_resids = field_plotter.muDiff_of_z_funct(zs_to_plot, RAs_for_resids, Decs_to_plot, zs_to_plot, RAs_to_plot, Decs_to_plot, muResids_to_plot, mu_errs_to_plot, [central_RA, central_Dec], [central_coord[0], 0.0, 0.0, fitted_unitless_density_param, fitted_radius_power])
        #print ('[len(zs_to_plot), len(RAs_to_plot), len(Decs_to_plot)] = ' + str([len(zs_to_plot), len(RAs_to_plot), len(Decs_to_plot)]))
        fitted_center_line_resids = field_plotter.muDiff_of_z_funct(center_line_zs, field_plotter.centerRAs([central_coord[1] for z in center_line_zs], RA_offset_to_center), [central_coord[2] for z in center_line_zs], zs_to_plot, RAs_to_plot, Decs_to_plot, muResids_to_plot, mu_errs_to_plot, [central_RA, central_Dec], [central_coord[0], fitted_unitless_density_param, fitted_radius_power, 0.0, 0.0], field_plotter.fit_funct)
        null_center_line_resids = field_plotter.muDiff_of_z_funct(center_line_zs, field_plotter.centerRAs([central_coord[1] for z in center_line_zs], RA_offset_to_center), [central_coord[2] for z in center_line_zs], zs_to_plot, RAs_to_plot, Decs_to_plot, muResids_to_plot, mu_errs_to_plot, [central_RA, central_Dec], [central_coord[0], 0.0, fitted_radius_power, 0.0, 0.0], field_plotter.fit_funct)
        null_resids = field_plotter.muDiff_of_z_funct(zs_to_plot, RAs_for_resids, Decs_to_plot, zs_to_plot, RAs_to_plot, Decs_to_plot, muResids_to_plot, mu_errs_to_plot, [central_RA, central_Dec], [central_coord[0], 0.0, 0.0, 0.0, 2.0], field_plotter.fit_funct)
        ax0, ax1 = [axarr[0][i], axarr[1][i]]
        void_curve = ax0.plot(center_line_zs, fitted_center_line_resids, c = 'k', alpha = 0.5)[0]
        null_curve = ax0.plot(center_line_zs, null_center_line_resids, c = 'k', alpha = 0.5, linestyle = '--')[0]
        scats = []
        for survey in unique_surveys_to_plot:
            zs_in_survey = [zs_to_plot[i] for i in range(len(zs_to_plot)) if surveys_to_plot[i] == survey]
            muResids_in_survey = [muResids_to_plot[i] for i in range(len(zs_to_plot)) if surveys_to_plot[i] == survey]
            mu_errs_in_survey = [mu_errs_to_plot[i] for i in range(len(zs_to_plot)) if surveys_to_plot[i] == survey]
            scats = scats + [ax0.scatter(zs_in_survey, muResids_in_survey, c = field_plotter.survey_to_color_dict[survey], marker = 'o')]
            ax0.errorbar(zs_in_survey, muResids_in_survey, yerr = mu_errs_in_survey, colors = field_plotter.survey_to_color_dict[survey], fmt = 'none', ecolor = field_plotter.survey_to_color_dict[survey])
        #ax0.scatter(zs_to_plot, muResids_to_plot, c = colors_to_plot, marker = 'o')
        fitted_data = ax0.scatter(zs_to_plot, fitted_resids, c = 'k', marker = 'x', alpha = 0.75)
        #[ ax0.annotate(str(j + 1), (zs_to_plot[j], muResids_to_plot[j]), color = 'k', fontsize = ticksize, verticalalignment = 'center', horizontalalignment = 'center') for j in range(len(sn_indeces_to_plot)) ]
        void_center_line = ax0.axvline(central_coord[0], color = 'k', linestyle = '--')
        ax0.legend([null_curve, void_curve, fitted_data, void_center_line], ['Fit through center - no void', 'Fit through center - with void', 'SNe Residuals from Void', 'Void center'])
        ax1.scatter(RAs_to_plot, Decs_to_plot, c = colors_to_plot, marker = 'o')
        #[ ax1.annotate(str(j + 1), (RAs_to_plot[j], Decs_to_plot[j]), color = 'k', fontsize = ticksize, verticalalignment = 'center', horizontalalignment = 'center' ) for j in range(len(sn_indeces_to_plot)) ]
        void_RA_Dec = ax1.scatter(central_coord[1], central_coord[2], c = 'k', marker = 'x', s = 100 )
        ax0.set_xlabel(r'$z$', fontsize = labelsize)
        ax1.set_xlabel(r'RA (deg)', fontsize = labelsize)
        ax0.set_ylabel(r'$\Delta \mu$ (mag)', fontsize = labelsize)
        ax1.set_ylabel(r'Dec (deg)', fontsize = labelsize)
        #ax0.set_xscale('log')
        ax0.tick_params(axis='both', labelsize= ticksize)
        ax1.tick_params(axis='both', labelsize= ticksize)
        label_strs = field_plotter.param_info_dict['label_strs']
        unit_strs = field_plotter.param_info_dict['units']
        fitted_params = field_plotter.param_conversion_funct([fitted_unitless_density_param, fitted_radius_power])
        val_text = '\n'.join([label_strs[i] + r'$=$' + str(can.round_to_n(fitted_params[i], 3)) + ' ' + unit_strs[i]  for i in range(len(fitted_params))]) + '\n' + r'$\chi^2_\nu / \chi^2_{\nu, 0} = $' + str(can.round_to_n(fits[i][1] / fits[i][2], 5))
        ax0.text(0.05, 0.03, val_text, fontsize = labelsize, transform = ax0.transAxes, verticalalignment = 'bottom', horizontalalignment = 'left', color = 'k')
        scats_for_legend = scats + [void_RA_Dec]
        surveys_for_legend = unique_surveys_to_plot + ['Void center']
        ax1.legend(scats_for_legend, surveys_for_legend)
        #plt.tight_layout()
        plt.suptitle('CMB Cold Spot Void Fit, Comoving bin of '+ str(can.round_to_n(comoving_bin, 3)) + ' Mpc')
        plt.savefig('/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/SNIsotropyProject/plots/' + 'CMBColdSpotFit_ComovingBin_' + str(comoving_bin) + 'Mpc_delta0_' + str(can.round_to_n(fitted_unitless_density_param, 2))  + 'r0_' + str(can.round_to_n(10.0 ** fitted_radius_power, 2)) + 'MpcVoid.pdf')
        plt.show()
        #field_plotter.makePlotOfOneDFits(points_to_test, fits, fig_size_unit = 2.5, save_plot_prefix = save_prefix, comoving_bin = comoving_bin)


import CosmicFitterClass as cfc
import numpy as np
import matplotlib.pyplot as plt
import cantrips as can
import time
import sys
import os

if __name__ == "__main__" :
    """
    Run MCMCs on artificial Pantheon+ - like data sets of toy SNe Ia.

    $ python MeasureCosmicParamContoursFromToyModels.py LowZ 0 0.01 Revision0 Revision0
    $ python MeasureCosmicParamContoursFromToyModels.py MidZ 0 0.01 Revision0 Revision0
    $ python MeasureCosmicParamContoursFromToyModels.py HighZ 0 0.01 Revision0 Revision0
    $ python MeasureCosmicParamContoursFromToyModels.py LowZ 0 0.01 Revision0 Revision1 # <=== This run is to generate an result for no additional surveys; the 0 after LowZ tells the code to run with no additional surveys
    $ python MeasureCosmicParamContoursFromToyModels.py LowZ 1 0.01 Revision0 Revision1
    $ python MeasureCosmicParamContoursFromToyModels.py MidZ 1 0.01 Revision0 Revision1
    $ python MeasureCosmicParamContoursFromToyModels.py HighZ 1 0.01 Revision0 Revision1

    [REPEAT FOR ALL MU OFFSETS OF INTEREST - TYPICALLY 0.005, 0.01, 0.025, 0.05]
    """
    cl_args = sys.argv[1:]
    #top_survey = cl_args[0]
    #top_surveys = [top_survey]
    toy_sn_data_file_suffix = cl_args[0]
    extra_survey_to_include_in_sequence = int(cl_args[1])
    extra_surveys_to_include_in_sequence = [extra_survey_to_include_in_sequence ]
    mu_offset = float(cl_args[2])
    mu_offsets = [mu_offset]
    print ('mu_offsets = ' + str(mu_offsets))
    base_cosmology_suffix = cl_args[3]
    run_id = cl_args[4]
    run_ids = [run_id]

    fiducial_n_scalings = 1
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


    #print ('top_surveys = ' + str(top_surveys))
    print ('Starting to fit Hubble constant from the command line... ')
    dir_base = '/Users/sashabrownsberger/Documents/Harvard/physics/'
    results_dir = dir_base + 'stubbs/variableMuFits/mcmcResults/ToyModel/'
    plot_dir = dir_base + 'stubbs/variableMuFits/plots/ToyModel/'
    #mu_offsets = [0.001, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.3, 0.5]
    #mu_offsets = [0.001, 0.003, 0.005, 0.007, 0.009, 0.011, 0.013, 0.015, 0.017, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2]
    #mu_offsets = [0.001, 0.005, 0.01, 0.015, 0.02, 0.05]
    z_lims = [0.0, np.inf]
    z_lims_str = 'zLim_Full'
    n_sn_per_survey = 300
    #z_lims = [0.0, 0.15]
    #z_lims_str = 'zLim_0p0_0p15'
    cosmicParams_plot_file_name = 'H0OmMw0_Unertainties_VsMuOffsetPrior_' + z_lims_str + '_ToyData.pdf'
    #mu_offsets = [0.0002, 0.2]
    posterior_fits = [0.0 for i in range(len(mu_offsets))]
    overplot_steps = 1500
    #overplot_steps = 500
    overplot_n_chains = 6
    full_steps = 10000
    full_steps = 1500
    full_n_chains = 50
    #full_n_chains = 10
    #surveys = ['CANDELS', 'CFA1',  'CFA2', 'CFA3K', 'CFA3S', 'CFA4p1', 'CFA4p2', 'CSP', 'DES', 'FOUNDATION', 'HST', 'KAIT', 'LOWZ', 'PS1MD', 'SDSS', 'SNAP', 'SNLS', 'SWIFT', 'SWIFTNEW']
    #n_iterations = 10
    cosmic_params_to_vary = ['H0','OmM','w0']
    #cosmic_params_to_vary = ['H0']
    init_guess_params_to_vary = [70.0, 0.3, -1.0]
    #init_guess_params_to_vary = [70.0]

    #all_surveys = ['A','B','C','D','E','F','G','H','I','J', 'H', 'I','J','K','L','M','N','O','P', 'Q']
    #base_surveys = ['CANDELS', 'CFA1', 'CFA2', 'CFA3K', 'CFA3S', 'CFA4p1', 'CFA4p2', 'CSP', 'DES', 'HST', 'KAIT', 'LOWZ', 'PS1MD', 'SDSS', 'SNAP', 'SNLS', 'SWIFT']
    base_surveys = ['SWIFT', 'ASASSN', 'CFA2', 'CFA1', 'KAITM', 'LOWZ', 'CFA4p2', 'KAIT', 'CFA3S', 'CSP', 'CFA3K', 'CFA4p1', 'PS1MD', 'SDSS', 'DES', 'SNLS', 'HST', 'SNAP', 'CANDELS']
    #base_surveys = ['SWIFT',    'ASASSN',   'CFA1',   'CFA2',   'LOWZ',    'KAITM',   'CFA4p2',  'KAIT',  'CFA3S',  'CSP',  'CFA3K',  'CFA4p1',   'PS1MD',  'SDSS',  'DES',  'SNLS',  'HST',  'SNAP',  'CANDELS']
    #base_surveys = ['CFA1', 'CFA2', 'CFA3K', 'CFA3S', 'CFA4p1', 'CFA4p2', 'CSP', 'KAIT','LOWZ', 'PS1MD', 'SWIFT']
    #base_surveys = ['SWIFT', 'ASASSN', 'CFA1']
    base_params_to_vary = cosmic_params_to_vary + ['mu' + survey for survey in base_surveys]
    if toy_sn_data_file_suffix.lower() == 'lowz':
        all_extra_surveys = [toy_sn_data_file_suffix]
        z_ranges = [[0.01, 0.1]]
    elif toy_sn_data_file_suffix.lower() == 'midz':
        all_extra_surveys = [toy_sn_data_file_suffix]
        z_ranges = [[0.1, 0.4]]
    elif toy_sn_data_file_suffix.lower() == 'highz':
        all_extra_surveys = [toy_sn_data_file_suffix]
        z_ranges = [[0.4, 1.5]]
    else:
        all_extra_surveys = ['A','B','C','D','E','F','G','H','I','J']
        z_ranges = [[-np.inf, np.inf]]
    print ('[all_extra_surveys, z_ranges] = ' + str([all_extra_surveys, z_ranges] ))
    #top_surveys = ['A','B']
    survey_errs_file = dir_base + 'stubbs/variableMuFits/ArtificialSurveys/ArtificialSNe_' + 'BaseErrs' + '.csv'
    base_cosmology_file = 'ToySNe_GPW' + str(int(1000 * mu_offset)) + 'mMags_' + z_lims_str + '_' + 'MCMC_' + 'toy_surveys' + '_REAL' + '_'.join(base_params_to_vary) + '_NS' + str(full_steps) + '_NC' + str(full_n_chains) + '_' + toy_sn_data_file_suffix + base_cosmology_suffix + '.txt'
    mcmc_outputs_dir = dir_base + 'stubbs/variableMuFits/mcmcResults/ToyModel/'
    if extra_survey_to_include_in_sequence > 0 and not(os.path.isfile(mcmc_outputs_dir + base_cosmology_file)):
        print ('We do not have the basic mcmc result file for mu offset ' + str(mu_offset) + 'mags, when we run under these same conditions without any additional surveys. ')
        print ('Looking for file ' + str(mcmc_outputs_dir + base_cosmology_file))
        print ('You should rerun the analysis as: $ python MeasureCosmicParamContoursFromToyModels.py ' + str(toy_sn_data_file_suffix) +  ' 0 ' + str(mu_offset) + ' ' + base_cosmology_suffix + ' ' + base_cosmology_suffix)
        sys.exit()
    elif extra_survey_to_include_in_sequence == 0:
        #If we're just running the analysis with no additional surveys,
        # we can just use dummy cosmic parameters.
        H0, OmegaM, w = [70.0, 0.3, 1.0]
    else:
        #Otherwise, we need to read in the best fit parameters from the reference file
        reference_data = can.readInColumnsToList(mcmc_outputs_dir + base_cosmology_file, delimiter = ', ', n_ignore = 1)
        H0_col, OmM_col, w_col = [reference_data[0], reference_data[1], reference_data[2]]
        H0_col = [float(elem) - float(H0_col[-1]) for elem in H0_col[0:-1]]
        OmM_col = [float(elem) - float(OmM_col[-1]) for elem in OmM_col[0:-1]]
        w_col = [float(elem)  - float(w_col[-1]) for elem in w_col[0:-1]]
        H0 = np.median(H0_col)
        OmegaM = np.median(OmM_col)
        w = np.median(w_col)

    surveys_file = 'stubbs/variableMuFits/ArtificialSurveys/ArtificialSNe_muOffset' + str(int(1000 * mu_offset)) + '_' + toy_sn_data_file_suffix + '.csv'
    covariance_file = 'stubbs/variableMuFits/ArtificialSurveys/ArtificialSNe_covariances_extra' + str(len(all_extra_surveys) * n_sn_per_survey) + '.txt'
    #surveys_file = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/variableMuFits/ArtificialSurveys/ArtificialSNe' + '_' + toy_sn_data_file_suffix + '.csv'
    if not(os.path.isfile(surveys_file)):
        #Need to make a new cosmology file
        print ('Making my own, new data file of SNe...')
        cfc.generateToySNeData(randomize_orig_sn = 0,
                               n_sn_per_survey = n_sn_per_survey, surveys = all_extra_surveys,
                               H0 = H0, OmegaM = OmegaM, Omega0 = 1.0, OmegaR = 0.0001, w_params = [w],
                               sn_data_type = 'pantheon_plus', cepheid_file_name = 'calibratorset.txt',
                               HF_surveys_to_include = 'all', z_ranges = z_ranges, mu_err = 'median',
                               colors = ['r','b','g','orange','cyan','magenta','purple','darkblue','darkgreen', 'skyblue',],
                               art_data_file = surveys_file, dir_base = dir_base)
    if not(os.path.isfile(dir_base + covariance_file)):
        #Need to make a new cosmology file
        print ('Making my own, new data file of SNe...')
        cfc.generateToySNeCovarianceFile(len(all_extra_surveys) * n_sn_per_survey, covariance_file, dir_base = dir_base)

    print ('surveys_file = ' + str(surveys_file))
    toy_data = can.readInColumnsToList(dir_base + surveys_file, delimiter = ', ', n_ignore = 1)
    surveys_col = 1
    colors_col = 5
    colors_dict = {toy_data[surveys_col][i]:toy_data[colors_col][i] for i in range(len(toy_data[0]))}
    colors = [colors_dict[survey] for survey in all_extra_surveys]
    print ('colors = ' + str(colors))
    mu_prior_type = 'gauss'
    all_cols_to_save = []
    #surveys = ['SDSS']
    #extra_surveys_to_include_in_sequence = [0] + [all_extra_surveys.index(survey)+1 for survey in top_surveys]
    #print ('extra_surveys_to_include_in_sequence = ' + str(extra_surveys_to_include_in_sequence))
    f, axarr = plt.subplots(3,1)
    H0_scats = []
    OmM_scats = []
    w_scats = []
    for j in range(len(extra_surveys_to_include_in_sequence)):
        n_extra_surveys = extra_surveys_to_include_in_sequence[j]
        color = colors[n_extra_surveys-1]
        extra_surveys_to_include = all_extra_surveys[0:n_extra_surveys]
        surveys_to_include = base_surveys + extra_surveys_to_include
        params_to_vary = cosmic_params_to_vary + ['mu' + survey for survey in surveys_to_include]
        #params_to_vary = ['H0','OmM','w0']
        save_files_dict = {mu_offset:['ToySNe_GPW' + str(int(1000 * mu_offset)) + 'mMags_' + z_lims_str + '_' + 'MCMC_' + 'toy_surveys' + '_REAL' + '_'.join(params_to_vary) + '_NS' + str(full_steps) + '_NC' + str(full_n_chains) + '_' + toy_sn_data_file_suffix + str(run_ids[iter]) + '.txt' for iter in range(len(run_ids))]}
        for iter in range(len(run_ids)):
            run_id = run_ids[iter]
            start_time = time.time()
            print ('Working on surveys ' + str(j) + ' of ' + str(len(extra_surveys_to_include_in_sequence)))
            save_fig = 'PosteriorFitsCosmicParamsToyData_' + 'NSurveys' + str(n_extra_surveys) + '_' + z_lims_str + '.pdf'
            toy_fitter = cfc.CosmicFitter(w_of_funct_str = 'w0', randomize_sn = 0, sn_data_type = 'toy_surveys', sn_toy_data_file = surveys_file,
                                          params_to_overplot = cosmic_params_to_vary, overplot_mcmc_steps = overplot_steps, overplot_mcmc_chains = overplot_n_chains, mu_prior_type = mu_prior_type,
                                          muOffsetPriors = {survey:[0.0, 0.05] for survey in base_surveys + all_extra_surveys}, fit_res_dir = 'mcmcResults/ToyModel/',
                                          covariance_file = covariance_file, dir_base = dir_base)

            for survey in base_surveys:
                mu_offset_scaling = mu_offsets_scalings_dict[survey]
                toy_fitter.muOffsetPriors[survey] = [0.0, fiducial_n_scalings * mu_offset_normalization * mu_offset_scaling]
            print ('toy_fitter.surveys = ' + str(toy_fitter.surveys))
            print ('surveys_to_include = ' + str(surveys_to_include))
            toy_fitter.updateUsedSN(z_lims = z_lims, surveys_to_include = surveys_to_include)
            #surveys = list(toy_fitter.surveys)
            for i in range(len(mu_offsets)):
                mu_offset = mu_offsets[i]
                print ('Working on mu_offset = ' + str(mu_offset))
                for survey in surveys_to_include:
                    toy_fitter.muOffsetPriors[survey] = [0.0, mu_offset]
                toy_fitter.doFullCosmicFit(init_guess_params_to_vary + [0.0 for survey in surveys_to_include], cosmic_params_to_vary + ['mu' + survey for survey in surveys_to_include], n_mcmc_steps = full_steps, n_mcmc_chains = full_n_chains , verbose = 0, additional_save_prefix = 'ToySNe_GPW' + str(int(1000 * mu_offset)) + 'mMags_' + z_lims_str + '_', additional_save_suffix = '_' + toy_sn_data_file_suffix + str(run_id), show_mcmc = 0, save_mcmc = 1, save_full_mcmc = 1)
                #toy_fitter.doFullCosmicFit([70.0, 0.3, -1.0] , ['H0', 'OmM', 'w0'] , n_mcmc_steps = full_steps, n_mcmc_chains = full_n_chains , verbose = 0, additional_save_prefix = 'ToySNe_GaussPriorWidth' + str(int(1000 * mu_offset)) + 'mMags_' + z_lims_str + '_', additional_save_suffix = '_' + str(run_id), show_mcmc = 0, save_mcmc = 1, save_full_mcmc = 1)
                #posterior_fits[i] = toy_fitter.computeBestFitResultsFromMCMC(['H0','OmM', 'w0'] + ['mu' + survey for survey in surveys], verbose_weighted_means = 0)

            """
            cols_to_save = [mu_offsets]
            for i in range(len(['H0', 'OmM', 'w0'] + surveys)):
                print ('i = ' + str(i))
                cols_to_save = cols_to_save + [[can.round_to_n(posterior_fit[0][i], 6) for posterior_fit in posterior_fits] , [can.round_to_n(posterior_fit[1][i][0],6) for posterior_fit in posterior_fits], [can.round_to_n(posterior_fit[1][i][1], 6) for posterior_fit in posterior_fits]]
            header = 'DeltaMu, ' + ', '.join([param + 'mu, -' + param + 'sig, +' + param + 'sig' for param in ['H0', 'OmM'] + ['mu' + survey for survey in surveys]])
            can.saveListsToColumns(cols_to_save, save_file, results_dir, sep = ', ', append = False, header = header, type_casts = None)
            """
        print ('save_files_dict = ' + str(save_files_dict))
        H0_scat = cfc.makePlotOfPosteriorWidths(save_files_dict, 0, ax = axarr[0], xlabel = '$\Delta \mu_S$ Gaussian prior width (mags)', ylabel = r'$H_0$ MCMC $\sigma$ (km/s/Mpc)',
                                      data_load_dir = results_dir, plot_save_dir = plot_dir, save_fig = 0, save_fig_name = save_fig, plot_marker = 'x', plot_color = color, show_all_chain_posteriors = 1)
        OmM_scat = cfc.makePlotOfPosteriorWidths(save_files_dict, 1, ax = axarr[1], xlabel = '$\Delta \mu_S$ Gaussian prior width (mags)', ylabel = r'$\Omega_M$ MCMC $\sigma$',
                                      data_load_dir = results_dir, plot_save_dir = plot_dir, save_fig = 0, save_fig_name = save_fig, plot_marker = 'x', plot_color = color, show_all_chain_posteriors = 1)
        w_scat = cfc.makePlotOfPosteriorWidths(save_files_dict, 2, ax = axarr[2], xlabel = '$\Delta \mu_S$ Gaussian prior width (mags)', ylabel = r'$w$ MCMC $\sigma$',
                                      data_load_dir = results_dir, plot_save_dir = plot_dir, save_fig = 0, save_fig_name = save_fig, plot_marker = 'x', plot_color = color, show_all_chain_posteriors = 1)
        H0_scats = H0_scats + [H0_scat]
        OmM_scats = OmM_scats + [OmM_scat]
        w_scats = w_scats + [w_scat]
        #all_cols_to_save = all_cols_to_save + [cols_to_save]
        end_time = time.time()
        print ('Took ' + str(end_time - start_time) + 's')

    #If this is our first time doing this, we need ot delete the data file.
    if extra_survey_to_include_in_sequence == 0:
        os.remove(dir_base + surveys_file)

    #axarr[0].legend(H0_scats, [str(n_surveys) + ' surveys' for n_surveys in surveys_to_include_in_sequence])
    #axarr[1].legend(OmM_scats, [str(n_surveys) + ' surveys' for n_surveys in surveys_to_include_in_sequence])
    #axarr[2].legend(w_scats, [str(n_surveys) + ' surveys' for n_surveys in surveys_to_include_in_sequence])
    #plt.savefig(plot_dir + save_fig)
    #plt.show()
    #plt.close('all')

    """
    f, axarr = plt.subplots( len(all_surveys), 3, figsize = (12.0, 3.0 * len(all_surveys)) )
    for n_surveys in range(1, len(all_surveys) + 1):
        cols_to_save = all_cols_to_save[n_surveys - 1]
        H0_neg_errs = cols_to_save[2]
        H0_plus_errs = cols_to_save[3]
        OmM_neg_errs = cols_to_save[5]
        OmM_plus_errs = cols_to_save[6]
        w0_neg_errs = cols_to_save[8]
        w0_plus_errs = cols_to_save[9]
        axarr[n_surveys - 1, 0].scatter(mu_offsets, (np.array(H0_neg_errs) + np.array(H0_plus_errs)) / 2.0, c = 'k', marker = 'x')
        axarr[n_surveys - 1, 0].text(0.5, 0.85, r'$N_{S}$=' + str(n_surveys), transform=axarr[n_surveys - 1, 0].transAxes)
        axarr[n_surveys - 1, 0].set_xlabel(r'$\Delta \mu_S$ prior width (mags)', fontsize = 12)
        axarr[n_surveys - 1, 0].set_ylabel(r'$H_0$ MCMC $\sigma$ (km/s/Mpc)', fontsize  = 12)
        axarr[n_surveys - 1, 1].scatter(mu_offsets, (np.array(OmM_neg_errs) + np.array(OmM_plus_errs)) / 2.0, c = 'k', marker = 'x')
        axarr[n_surveys - 1, 1].set_xlabel(r'$\Delta \mu_S$ prior width (mags)', fontsize = 12)
        axarr[n_surveys - 1, 1].set_ylabel(r'$\Omega_M$ MCMC $\sigma$', fontsize  = 12)
        axarr[n_surveys - 1, 1].text(0.5, 0.85, r'$N_{S}$=' + str(n_surveys), transform=axarr[n_surveys - 1, 1].transAxes)
        axarr[n_surveys - 1, 2].scatter(mu_offsets, (np.array(w0_neg_errs) + np.array(w0_plus_errs)) / 2.0, c = 'k', marker = 'x')
        axarr[n_surveys - 1, 2].set_xlabel(r'$\Delta \mu_S$ prior width (mags)', fontsize = 12)
        axarr[n_surveys - 1, 2].set_ylabel(r'$w_{\Lambda}$ MCMC $\sigma$', fontsize  = 12)
        axarr[n_surveys - 1, 2].text(0.5, 0.85, r'$N_{S}$=' + str(n_surveys), transform=axarr[n_surveys - 1, 2].transAxes)
    """

    plt.subplots_adjust (wspace=0.5, hspace=0.5)
    plt.savefig(plot_dir +  cosmicParams_plot_file_name)

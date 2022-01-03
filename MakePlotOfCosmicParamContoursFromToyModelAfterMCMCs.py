
import CosmicFitterClass as cfc
import numpy as np
import matplotlib.pyplot as plt
import cantrips as can
import time
import sys

if __name__=="__main__":
    """
    $ python MakePlotOfCosmicParamContoursFromToyModelAfterMCMCs.py 0 A B C
    """
    cl_args = sys.argv[1:]
    surveys_to_plot = cl_args
    results_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/variableMuFits/mcmcResults/ToyModel/'
    plot_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/variableMuFits/plots/ToyModel/'
    mu_offsets = [0.001, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.3, 0.5]
    #mu_offsets = [0.001, 0.003, 0.005, 0.007, 0.009, 0.011, 0.013, 0.015, 0.017, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2]
    mu_offsets = [0.005, 0.01, 0.02, 0.06 ]
    z_lims_str = 'zLim_Full'
    #z_lims = [0.0, 0.15]
    #z_lims_str = 'zLim_0p0_0p15'
    cosmicParams_plot_file_name = 'H0OmMw0_Unertainties_VsMuOffsetPrior_' + z_lims_str + '_ToyData.pdf'
    #mu_offsets = [0.0002, 0.2]
    full_steps = 10000
    #full_steps = 2000
    full_n_chains = 50
    #full_n_chains = 10
    #surveys = ['CANDELS', 'CFA1',  'CFA2', 'CFA3K', 'CFA3S', 'CFA4p1', 'CFA4p2', 'CSP', 'DES', 'FOUNDATION', 'HST', 'KAIT', 'LOWZ', 'PS1MD', 'SDSS', 'SNAP', 'SNLS', 'SWIFT', 'SWIFTNEW']
    n_iterations = 1

    #all_surveys = ['A','B','C','D','E','F','G','H','I','J', 'H', 'I','J','K','L','M','N','O','P', 'Q']
    base_surveys = ['SWIFT', 'CFA2', 'KAIT', 'CFA1', 'LOWZ', 'CFA4p2', 'CFA3S', 'CSP' ,'CFA4p1', 'CFA3K', 'PS1MD', 'SDSS', 'DES', 'SNLS', 'HST', 'SNAP', 'CANDELS']
    all_surveys_set = [['LowZ'], ['MidZ'], ['HighZ']]

    surveys_files = ['/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/variableMuFits/ArtificialSurveys/ArtificialSNe_LowZ.csv',
                     '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/variableMuFits/ArtificialSurveys/ArtificialSNe_MidZ.csv',
                     '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/variableMuFits/ArtificialSurveys/ArtificialSNe_HighZ.csv' ]
    toy_data_sets = [can.readInColumnsToList(surveys_file, delimiter = ', ', n_ignore = 1) for surveys_file in surveys_files]
    surveys_col = 1
    colors_col = 5
    colors_dicts = [{toy_data[surveys_col][i]:toy_data[colors_col][i] for i in range(len(toy_data[0]))} for toy_data in toy_data_sets]
    colors = [[colors_dicts[j][survey] for survey in all_surveys_set[j]] for j in range(len(all_surveys_set))]
    print ('colors = ' + str(colors))
    mu_prior_type = 'gauss'
    all_cols_to_save = []
    #surveys = ['SDSS']
    extra_surveys_to_include_in_sequence_set = [ [] for surveys_to_plot in all_surveys_set ]
    """
    for i in range(len(all_surveys_set)):
        surveys_to_plot = all_surveys_set[i]
        for survey in surveys_to_plot:
            if survey in all_surveys :
                extra_surveys_to_include_in_sequence_set[i] = extra_surveys_to_include_in_sequence_set[i] + [all_surveys.index(survey) + 1]
            else:
                extra_surveys_to_include_in_sequence_set[i] = extra_surveys_to_include_in_sequence_set[i]  + [0]
    """ 
    #extra_surveys_to_include_in_sequence = [all_surveys.index(survey)+1 for survey in surveys_to_plot]
    print ('extra_surveys_to_include_in_sequence_set = ' + str(extra_surveys_to_include_in_sequence_set))
    f, axarr = plt.subplots(3,1, figsize = (5,6))
    H0_ylims, OmM_ylims, w_ylims = [ [0.21, 0.34], [0.28, 0.62], [0.07, 0.21] ]
    labelsize = 10
    H0_scats = []
    OmM_scats = []
    w_scats = []

    for j in range(len(extra_surveys_to_include_in_sequence_set)):
        n_extra_surveys = extra_surveys_to_include_in_sequence_set[j][0]
        if n_extra_surveys> 0:
            color = colors[j][n_extra_surveys-1]
        else:
            color = 'k'

        surveys_to_include = all_surveys[j][0:n_extra_surveys]
        params_to_vary = ['H0','OmM','w0'] + ['mu' + survey for survey in base_surveys] + ['mu' + survey for survey in surveys_to_include]
        save_files_dict = {mu_offset:['ToySNe_GaussPriorWidth' + str(int(1000 * mu_offset)) + 'mMags_' + z_lims_str + '_' + 'MCMC_' + 'toy_surveys' + '_REAL' + '_'.join(params_to_vary) + '_NS' + str(full_steps) + '_NC' + str(full_n_chains) + '_' + str(iter) + '.txt' for iter in range(1, n_iterations + 1)] for mu_offset in mu_offsets}
        print ('save_files_dict = ' + str(save_files_dict))

        axarr[0].set_ylim(H0_ylims)
        #axarr[1].set_ylim(OmM_ylims)
        axarr[2].set_ylim(w_ylims)
        #H0_scat = cfc.makePlotOfPosteriorWidths(save_files_dict, 0, ax = axarr[0], xlabel = '$\Delta \mu_S$ Gaussian prior width (mags)', ylabel = r'$H_0$ MCMC $\sigma$ (km/s/Mpc)',
        #                              ref_param_vals_for_other_ax = (None if j < len(extra_surveys_to_include_in_sequence) - 1 else [70.0, 67.4, 0.5, r"$H_0$ tension (# of $\sigma$s)"]),
        #                              data_load_dir = results_dir, plot_save_dir = plot_dir, save_fig = 0, save_fig_name = 0, plot_marker = 'x', plot_color = color, show_all_chain_posteriors = 0,
        #                              zero_uncertainties_at_zero_prior = 0, labelsize = labelsize, ylims = H0_ylims)
        H0_scat = cfc.makePlotOfPosteriorWidths(save_files_dict, 0, ax = axarr[0], xlabel = '$\Delta \mu_S$ Gaussian prior width (mags)', ylabel = r'$H_0$ MCMC $\sigma$ (km/s/Mpc)',
                                      data_load_dir = results_dir, plot_save_dir = plot_dir, save_fig = 0, save_fig_name = 0, plot_marker = 'x', plot_color = color, show_all_chain_posteriors = 0,
                                      zero_uncertainties_at_zero_prior = 0, labelsize = labelsize, ylims = H0_ylims)
        H0_scats = H0_scats + [H0_scat]

        OmM_scat = cfc.makePlotOfPosteriorWidths(save_files_dict, 1, ax = axarr[1], xlabel = '$\Delta \mu_S$ Gaussian prior width (mags)', ylabel = r'$\Omega_M$ MCMC $\sigma$',
                                      data_load_dir = results_dir, plot_save_dir = plot_dir, save_fig = 0, save_fig_name = 0, plot_marker = 'x', plot_color = color, show_all_chain_posteriors = 0,
                                       zero_uncertainties_at_zero_prior = 0, labelsize = labelsize, ylims = OmM_ylims)
        OmM_scats = OmM_scats + [OmM_scat]
        w_scat = cfc.makePlotOfPosteriorWidths(save_files_dict, 2, ax = axarr[2], xlabel = '$\Delta \mu_S$ Gaussian prior width (mags)', ylabel = r'$w$ MCMC $\sigma$',
                                      data_load_dir = results_dir, plot_save_dir = plot_dir, save_fig = 0, save_fig_name = 0, plot_marker = 'x', plot_color = color, show_all_chain_posteriors = 0,
                                      zero_uncertainties_at_zero_prior = 0, labelsize = labelsize, ylims = w_ylims)

        w_scats = w_scats + [w_scat]

        #all_cols_to_save = all_cols_to_save + [cols_to_save]
    axarr[0].legend(H0_scats, ['+' + str(n_surveys) + ' surveys' for n_surveys in extra_surveys_to_include_in_sequence], ncol = 2)
    axarr[1].legend(OmM_scats, ['+' + str(n_surveys) + ' surveys' for n_surveys in extra_surveys_to_include_in_sequence], ncol = 2)
    axarr[2].legend(w_scats, ['+' + str(n_surveys) + ' surveys' for n_surveys in extra_surveys_to_include_in_sequence], ncol = 2)
    plt.subplots_adjust (wspace=0.5, hspace=0.5)
    plt.savefig(plot_dir + cosmicParams_plot_file_name)
    plt.show()

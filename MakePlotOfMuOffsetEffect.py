import numpy as np
import matplotlib.pyplot as plt
import random
import loadSN as lsn
import CosmicFitterClass as cfc
import SNDataArchive as snda


if __name__ == "__main__":
    sn_data_type = 'pantheon_plus'
    results_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/variableMuFits/mcmcResults/'
    plot_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/variableMuFits/plots/'
    z_lims = [0.0, np.inf]
    z_lims_str = 'zLim_Full'
    #z_lims = [0.0, 0.15]
    #z_lims_str = 'zLim_0p0_0p15'
    plot_file_name = 'H0UnertaintyVsMuOffsetPrior_' + z_lims_str + '.pdf'
    save_file = 'PosteriorFitsDifferentOffsets_' + z_lims_str + '.txt'
    overplot_steps = 1500
    overplot_n_chains = 4
    full_steps = 5000
    full_n_chains = 40
    sn_data_archive = snda.SNDataArchive()
    colors_by_survey = sn_data_archive.survey_color_map
    #surveys = ['CANDELS', 'CFA1',  'CFA2', 'CFA3K', 'CFA3S', 'CFA4p1', 'CFA4p2', 'CSP', 'DES', 'FOUNDATION', 'HST', 'KAIT', 'LOWZ', 'PS1MD', 'SDSS', 'SNAP', 'SNLS', 'SWIFT', 'SWIFT']
    surveys = ['CANDELS', 'CFA1',  'CFA2', 'CFA3K', 'CFA3S', 'CFA4p1', 'CFA4p2', 'CSP', 'DES', 'HST', 'KAIT', 'LOWZ', 'PS1MD', 'SDSS', 'SNAP', 'SNLS', 'SWIFT']
    surveys = ['SWIFT', 'CFA3K'] #We choose these surveys because SWIFT is the survey with the most Cepheid SNe
                                    #   and we have about as many Hubble Flow CFA3K SNe as we have Hubble Flow SWIFT SNe
    fig_name = 'TwoSurveySampleOfLuminosityOffsetEffect.pdf'
    #surveys = ['CFA1', 'CFA2', 'CFA3S', 'CSP', 'KAIT', 'LOWZ', 'PS1MD', 'SWIFT']
    mu_prior_type = 'gauss'
    base_H0 = 70.0
    base_OmM = 0.3
    base_w0 = -1.0
    ticklabelsize = 13
    labelsize = 15

    fancy_plot = 1
    if fancy_plot:
        # use LaTeX fonts in the plot
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')

    rand_fitter = cfc.CosmicFitter(w_of_funct_str = 'w0', randomize_sn = 1, sn_data_type = 'pantheon_plus', params_to_overplot = ['H0', 'OmM'], overplot_mcmc_steps = overplot_steps, overplot_mcmc_chains = overplot_n_chains, mu_prior_type = mu_prior_type,
                                   canonical_cosmic_vals = {'H0':base_H0, 'L':1.0, 'OmM':base_OmM, 'OmL':1.0 - base_OmM, 'OmR': 0.0001, 'Om0': 1.0, 'w0':base_w0,'wa':0}, force_rand_fit_to_background_cosmology = 0)
    rand_fitter.updateUsedSN(z_lims = z_lims, surveys_to_include  = surveys )

    scalar_params, wOfFunctionParams, muOffsets = rand_fitter.getAllFitParametersFromListOfParamsToVary([base_H0, base_OmM], ['H0','OmM'])
    H0, OmegaM, OmegaLambda, OmegaR, Omega0 = scalar_params
    print ('[H0, OmegaM, OmegaLambda, OmegaR, Omega0] = ' + str([H0, OmegaM, OmegaLambda, OmegaR, Omega0]))
    print ('[rand_fitter.H0, rand_fitter.OmegaM, rand_fitter.OmegaLambda, rand_fitter.OmegaR , rand_fitter.Omega0] = ' + str([rand_fitter.H0, rand_fitter.OmegaM, rand_fitter.OmegaLambda, rand_fitter.OmegaR , rand_fitter.Omega0]))
    print ('[wOfFunctionParams, rand_fitter.default_w_params] = ' + str([wOfFunctionParams, rand_fitter.default_w_params] ))
    loadedWOfFunction = lambda zs: rand_fitter.w_of_funct(zs, wOfFunctionParams[:])
    mu_offsets_by_survey_dict = rand_fitter.getMuOffsetsBySurvey(muOffsets)
    print ('mu_offsets_by_survey_dict = ' + str(mu_offsets_by_survey_dict ))
    #base_calc_mus_nonCeph = rand_fitter.getMusForCosmology(rand_fitter.nonCeph_zs, [rand_fitter.H0, rand_fitter.OmegaM, rand_fitter.OmegaLambda, rand_fitter.OmegaR , rand_fitter.Omega0], wOfFunction = lambda zs: rand_fitter.w_of_funct(zs, rand_fitter.default_w_params) )
    #base_calc_mus_Ceph = rand_fitter.getMusForCosmology([rand_fitter.sorted_zs[index] for index in rand_fitter.cepheid_indeces], [rand_fitter.H0, rand_fitter.OmegaM, rand_fitter.OmegaLambda, rand_fitter.OmegaR , rand_fitter.Omega0], wOfFunction = lambda zs: rand_fitter.w_of_funct(zs, rand_fitter.default_w_params))
    #base_calc_mus_all = rand_fitter.getMusForCosmology(rand_fitter.sorted_zs, [rand_fitter.H0, rand_fitter.OmegaM, rand_fitter.OmegaLambda, rand_fitter.OmegaR , rand_fitter.Omega0], wOfFunction = lambda zs: rand_fitter.w_of_funct(zs, rand_fitter.default_w_params) )
    f, axarr = plt.subplots(4,1, squeeze = False, figsize = [14, 8])
    ylims = [-1.5, 1.5]
    offset_dicts = [{'SWIFT':0.0, 'CFA3K':0.0}, {'SWIFT':0.0, 'CFA3K':1.0}, {'SWIFT':1.0, 'CFA3K':0.0}, {'SWIFT':1.0, 'CFA3K':1.0}]


    for i in range(len(offset_dicts)):
        mu_offsets_by_survey_dict = offset_dicts[i]
        mu_offsets_nonCeph = [mu_offsets_by_survey_dict [survey] for survey in rand_fitter.nonCeph_surveys]
        ceph_zs = [rand_fitter.sorted_zs[index] for index in rand_fitter.cepheid_indeces]
        ceph_mus = [rand_fitter.sorted_mus[index] for index in rand_fitter.cepheid_indeces]
        ceph_surveys = [rand_fitter.sorted_surveys[index] for index in rand_fitter.cepheid_indeces]
        ceph_muErrs = [rand_fitter.sorted_muErrs[index] for index in rand_fitter.cepheid_indeces]
        mu_offsets_Ceph = [mu_offsets_by_survey_dict [survey] for survey in ceph_surveys ]
        calc_mus_nonCeph = rand_fitter.getPredictedMus(rand_fitter.nonCeph_zs, [rand_fitter.sorted_muErrs[index] for index in rand_fitter.cepheid_indeces], [rand_fitter.sorted_surveys[index] for index in rand_fitter.cepheid_indeces], [rand_fitter.H0, rand_fitter.OmegaM, rand_fitter.OmegaLambda, rand_fitter.OmegaR , rand_fitter.Omega0], mu_offsets_by_survey_dict, rand_fitter.cepheid_indeces, loadedWOfFunction,  verbose = 1  )
        calc_mus_nonCeph = np.array(calc_mus_nonCeph) + np.array(mu_offsets_nonCeph)
        calc_mus_Ceph = rand_fitter.getPredictedMus(ceph_zs, [rand_fitter.sorted_muErrs[index] for index in rand_fitter.cepheid_indeces], [rand_fitter.sorted_surveys[index] for index in rand_fitter.cepheid_indeces], [rand_fitter.H0, rand_fitter.OmegaM, rand_fitter.OmegaLambda, rand_fitter.OmegaR , rand_fitter.Omega0], mu_offsets_by_survey_dict, rand_fitter.cepheid_indeces, loadedWOfFunction, verbose = 1  )
        calc_mus_Ceph = np.array(calc_mus_Ceph) + np.array(mu_offsets_Ceph)
        luminosity_mu_offset = -rand_fitter.determineLuminosityOffsetFromCepheids([rand_fitter.sorted_muErrs[index] for index in rand_fitter.cepheid_indeces], [rand_fitter.sorted_surveys[index] for index in rand_fitter.cepheid_indeces], mu_offsets_by_survey_dict, rand_fitter.cepheid_indeces )

        ceph_paired_sn = axarr[i,0].scatter([rand_fitter.nonCeph_zs[i] for i in range(len(rand_fitter.nonCeph_surveys)) if rand_fitter.nonCeph_surveys[i] == 'SWIFT'], [(np.array(calc_mus_nonCeph) - np.array(rand_fitter.nonCeph_mus))[i] for i in range(len(rand_fitter.nonCeph_surveys)) if rand_fitter.nonCeph_surveys[i] == 'SWIFT'], c =[rand_fitter.nonCeph_plotColors[i] for i in range(len(rand_fitter.nonCeph_surveys)) if rand_fitter.nonCeph_surveys[i] == 'SWIFT'] )
        noCeph_paired_sn = axarr[i,0].scatter([rand_fitter.nonCeph_zs[i] for i in range(len(rand_fitter.nonCeph_surveys)) if rand_fitter.nonCeph_surveys[i] == 'CFA3K'], [(np.array(calc_mus_nonCeph) - np.array(rand_fitter.nonCeph_mus))[i] for i in range(len(rand_fitter.nonCeph_surveys)) if rand_fitter.nonCeph_surveys[i] == 'CFA3K'], c =[rand_fitter.nonCeph_plotColors[i] for i in range(len(rand_fitter.nonCeph_surveys)) if rand_fitter.nonCeph_surveys[i] == 'CFA3K'] )
    #axarr[0,0].scatter(rand_fitter.nonCeph_zs, np.array(rand_fitter.nonCeph_mus) - np.array(calc_mus), c ='k', marker = '.'  )
        #axarr[i,0].scatter(ceph_zs, np.array(calc_mus_Ceph) - np.array(ceph_mus), c ='k', marker = 'x' )
        lum_correction_offset = axarr[i,0].axhline(luminosity_mu_offset, c = 'k', linestyle = '--', dashes = (5,1))
        ceph_paired_offset = axarr[i,0].axhline(mu_offsets_by_survey_dict['SWIFT'], c = colors_by_survey['SWIFT'], linestyle = '--', dashes = (5,7))
        noCeph_paired_offset = axarr[i,0].axhline(mu_offsets_by_survey_dict['CFA3K'], c = colors_by_survey['CFA3K'], linestyle = '--', dashes = (5,13))
        if i == len(offset_dicts) - 1:
            axarr[i,0].set_xlabel(r'$z$', fontsize = labelsize)
        else:
            axarr[i,0].set_xlabel('')
            axarr[i,0].set_xticks([])
        axarr[i,0].set_ylabel(r'$\mu - \mu_{\mathrm{model}}$ (mag)', fontsize = labelsize)
        axarr[i,0].set_ylim(ylims)
        if i == 0:
            axarr[i,0].legend([ceph_paired_sn, noCeph_paired_sn, ceph_paired_offset, noCeph_paired_offset, lum_correction_offset],
                              ['HF SNe in survey A (with Cepheid-paired SNe)', 'HF SNe in survey B (no Cepheid-paired SNe)', r'Forced $\Delta \mu$ in survey A', r'Forced $\Delta \mu$ in survey B', r'Additional $\Delta \mu$ from SNe Luminusity correction'],
                              ncol = 3)
        axarr[i,0].tick_params(axis='x', labelsize= ticklabelsize )
        axarr[i,0].tick_params(axis='y', labelsize= ticklabelsize )

    #axarr[0,0].scatter(rand_fitter.sorted_zs, np.array(rand_fitter.sorted_mus) - np.array(base_calc_mus_all), c ='k', marker = '.' )
    #axarr[0,1].scatter(rand_fitter.nonCeph_zs, np.array(rand_fitter.nonCeph_mus) , c =rand_fitter.nonCeph_plotColors )
    #axarr[0,1].scatter([rand_fitter.sorted_zs[i] for i in rand_fitter.cepheid_indeces], np.array([rand_fitter.sorted_mus[i] for i in rand_fitter.cepheid_indeces]), c ='k', marker = 'x' )

    #axarr[0,1].scatter(rand_fitter.sorted_zs, rand_fitter.sorted_mus , c ='k', marker = '.' )
    #axarr[0,0].scatter(rand_fitter.nonCeph_zs, np.array(rand_fitter.nonCeph_mus) - np.array(base_calc_mus_nonCephA), c =rand_fitter.nonCeph_plotColors )
    #axarr[0,0].scatter([rand_fitter.sorted_zs[i] for i in rand_fitter.cepheid_indeces], np.array([rand_fitter.nonCeph_mus[i] for i in rand_fitter.cepheid_indeces]) - np.array(base_calc_mus_CephA), c ='k', marker = 'x' )
    #axarr[0,1].scatter(rand_fitter.nonCeph_zs, np.array(rand_fitter.nonCeph_mus) , c =rand_fitter.nonCeph_plotColors )
    #axarr[0,1].scatter([rand_fitter.sorted_zs[i] for i in rand_fitter.cepheid_indeces], np.array([rand_fitter.nonCeph_mus[i] for i in rand_fitter.cepheid_indeces]) , c ='k', marker = 'x' )

    plt.savefig(plot_dir + fig_name)
    plt.show()

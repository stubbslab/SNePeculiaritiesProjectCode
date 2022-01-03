import cantrips as can
import CosmicFitterClass as cfc
import matplotlib.pyplot as plt
import numpy as np
import sys

if __name__=="__main__":
    """
    $ python makePlotOfCosmicParamConfidenceVsGaussianPriorWidth.py 1 1 2 4 5 6 7 8 9 11 12 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30
    $ python makePlotOfCosmicParamConfidenceVsGaussianPriorWidth.py 1 2
    # missing chains: 3, 10, 13 - 300, 11 - 300
    """
    cl_args = sys.argv[1:]
    #sequence_ids = cl_args[1:]
    real_sequence_ids = cl_args[:]
    #sequence_ids = [int(id) for id in sequence_ids]
    real_sequence_ids = [int(id) for id in real_sequence_ids]
    show_all_chains = 0

    #OneD_posterior_file_name = 'H0_OmM_w_UncertaintyVsMuOffsetPrior_zLim_Full_' + '_'.join([str(id) for id in real_sequence_ids]) + '.pdf'
    OneD_posterior_file_name = 'H0_OmM_UncertaintyVsMuOffsetPrior_zLim_Full_' + '_'.join([str(id) for id in real_sequence_ids]) + '.pdf'
    save_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/variableMuFits/plots/PosteriorWidths/'
    target_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/variableMuFits/mcmcResults/'
    #posterior_data_files = ['PosteriorFitsDifferentOffsets_zLim_Full_' + str(i) + '.txt' for i in range(sequence_start, sequence_end + 1)]
    mu_offsets = [0.001, 0.003, 0.005, 0.007, 0.009, 0.011, 0.013, 0.015, 0.017, 0.02, 0.04, 0.06, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2]
    mu_offsets = [0.001, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6]
    mu_offsets = [0.001, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6]
    #mu_offsets = [0.001, 0.06, 0.1, 0.2]
    OmM_lims, w0_lims = [[0.01, 0.52], [-1.7, -0.4]]
    #surveys = ['CANDELS', 'CFA1',  'CFA2', 'CFA3K', 'CFA3S', 'CFA4p1', 'CFA4p2', 'CSP', 'DES', 'HST', 'KAIT', 'LOWZ', 'PS1MD', 'SDSS', 'SNAP', 'SNLS', 'SWIFT']
    surveys = ['SWIFT', 'CFA2', 'KAIT', 'CFA1', 'LOWZ', 'CFA4p2', 'CFA3S', 'CSP' ,'CFA4p1', 'CFA3K', 'PS1MD', 'SDSS', 'DES', 'SNLS', 'HST', 'SNAP', 'CANDELS']
    #surveys = ['CANDELS', 'CFA1',  'CFA2', 'CFA3K', 'CFA3S']
    #rand_mcmc_data_files_dict = {muSig:['GaussPriorWidth' + str(int(muSig * 1000)) + 'mMags_zLim_Full_MCMC_pantheon_plus_RAND' + 'H0_OmM_w0' + '_mu' + '_mu'.join(surveys) + '_NS10000_NC40_' + str(i) + '.txt' for i in sequence_ids] for muSig in mu_offsets}
    #real_mcmc_data_files_dict = {muSig:['GaussPriorWidth' + str(int(muSig * 1000)) + 'mMags_zLim_Full_MCMC_pantheon_plus_REAL' + 'H0_OmM_w0' + '_mu' + '_mu'.join(surveys) + '_NS10000_NC50_Real' + str(real_sequence_id) + '.txt' for real_sequence_id in real_sequence_ids] for muSig in mu_offsets}
    real_mcmc_data_files_dict = {muSig:['GaussPriorWidth' + str(int(muSig * 1000)) + 'mMags_zLim_Full_MCMC_pantheon_plus_REAL' + 'H0_OmM' + '_mu' + '_mu'.join(surveys) + '_NS10000_NC50_' + str(real_sequence_id) + '.txt' for real_sequence_id in real_sequence_ids] for muSig in mu_offsets}
    f, axarr = plt.subplots(2, 1, figsize = (5,4))
    labelsize = 10
    #f, axarr = plt.subplots(1, 1, figsize = (5,3))
    #H0_ax = axarr
    H0_ax = axarr[0]
    Om_ax = axarr[1]
    #w0_ax = axarr[2]
    #cfc.makePlotOfPosteriorWidths(rand_mcmc_data_files_dict, 0, ax = H0_ax, ref_param_vals_for_other_ax = [67.4, 0.5, r"$H_0$ tension (# of $\sigma$s)"], xlabel = '$\Delta \mu_S$ Gaussian prior width (mags)', ylabel = r'Increase in $H_0$' + '\n'  + r'posterior $\sigma$ (km/s/Mpc)', data_load_dir = target_dir, plot_color = 'k', plot_marker = '.', show_fig = 0, save_fig = 0)
    #cfc.makePlotOfPosteriorWidths(rand_mcmc_data_files_dict, 1, ax = Om_ax, xlabel = '$\Delta \mu_S$ Gaussian prior width (mags)', ylabel = r'Increase in $\Omega_M$' + '\n'  + r'posterior $\sigma$', data_load_dir = target_dir, plot_color = 'k', plot_marker = '.', show_fig = 0, save_fig = 0)
    #cfc.makePlotOfPosteriorWidths(rand_mcmc_data_files_dict, 2, ax = w0_ax, xlabel = '$\Delta \mu_S$ Gaussian prior width (mags)', ylabel = r'Increase in $w$' + '\n'  + r'posterior $\sigma$', data_load_dir = target_dir, plot_color = 'k', plot_marker = '.', show_fig = 0, save_fig = 0)
    #cfc.makePlotOfPosteriorWidths(rand_mcmc_data_files_dict, 0, ax = H0_ax, xlabel = '$\Delta \mu_S$ Gaussian prior width (mags)', ylabel = r'$H_0$' + r' posterior $\sigma$ (km/s/Mpc)', data_load_dir = target_dir, plot_color = 'k', plot_marker = '.', show_fig = 0, save_fig = 0, zero_uncertainties_at_zero_prior = 0, labelsize = labelsize)
    print ('real_mcmc_data_files_dict = ' + str(real_mcmc_data_files_dict))
    #We need to add the pre-inter survey offset uncertainties in quadrature from Planck and SH0ES.
    cfc.makePlotOfPosteriorWidths(real_mcmc_data_files_dict, 0, ax = H0_ax, ref_param_vals_for_other_ax = [73.2, 67.27, np.sqrt(0.6 ** 2.0 + 1.3 ** 2.0), r"Planck vs SH0ES HT"], xlabel = '$\Delta \mu_S$ Gaussian prior width (mags)', ylabel = r'Increase in $H_0$' + '\n' + r' posterior $\sigma$ (km/s/Mpc)', data_load_dir = target_dir, plot_color = 'k', plot_marker = 'x', show_fig = 0, save_fig = 0, zero_uncertainties_at_zero_prior = 1, labelsize = labelsize, show_all_chain_posteriors = show_all_chains)
    #cfc.makePlotOfPosteriorWidths(rand_mcmc_data_files_dict, 1, ax = Om_ax, xlabel = '$\Delta \mu_S$ Gaussian prior width (mags)', ylabel = r'$\Omega_M$' + r' posterior $\sigma$', data_load_dir = target_dir, plot_color = 'k', plot_marker = '.', show_fig = 0, save_fig = 0, zero_uncertainties_at_zero_prior = 0, labelsize = labelsize)
    cfc.makePlotOfPosteriorWidths(real_mcmc_data_files_dict, 1, ax = Om_ax, xlabel = '$\Delta \mu_S$ Gaussian prior width (mags)', ylabel = r'Increase in $\Omega_M$' + '\n' + r' posterior $\sigma$', data_load_dir = target_dir, plot_color = 'k', plot_marker = 'x', show_fig = 0, save_fig = 0, zero_uncertainties_at_zero_prior = 1, labelsize = labelsize, show_all_chain_posteriors = show_all_chains)
    #cfc.makePlotOfPosteriorWidths(rand_mcmc_data_files_dict, 2, ax = w0_ax, xlabel = '$\Delta \mu_S$ Gaussian prior width (mags)', ylabel = r'$w$' + r' posterior $\sigma$', data_load_dir = target_dir, plot_color = 'k', plot_marker = '.', show_fig = 0, save_fig = 0, zero_uncertainties_at_zero_prior = 0, labelsize = labelsize)
    #cfc.makePlotOfPosteriorWidths(real_mcmc_data_files_dict, 2, ax = w0_ax, xlabel = '$\Delta \mu_S$ Gaussian prior width (mags)', ylabel = r'Increase in $w$' + '\n' + r' posterior $\sigma$', data_load_dir = target_dir, plot_color = 'k', plot_marker = 'x', show_fig = 0, save_fig = 0, zero_uncertainties_at_zero_prior = 1, labelsize = labelsize, show_all_chain_posteriors = show_all_chains)
    plt.subplots_adjust (left = 0.20, right = 0.85)
    plt.savefig(save_dir + OneD_posterior_file_name)
    plt.close('all')
    

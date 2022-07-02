import cantrips as can
import CosmicFitterClass as cfc
import matplotlib.pyplot as plt
import numpy as np
import sys

if __name__=="__main__":
    """
    $ python makePlotOfCosmicParamConfidenceVsGaussianPriorWidth.py Test1 1 2 4 5 6 7 8 9 11 12 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30
    $ python makePlotOfCosmicParamConfidenceVsGaussianPriorWidth.py Revised1 Revised2 Revised3 Revised4 Revised5 Revised6 Revised7 Revised8 Revised9 Revised10 Revised11
    # missing chains: 3, 10, 13 - 300, 11 - 300
    """
    cl_args = sys.argv[1:]
    #sequence_ids = cl_args[1:]
    real_sequence_ids = cl_args[:]
    #sequence_ids = [int(id) for id in sequence_ids]
    real_sequence_ids = [id for id in real_sequence_ids]
    show_all_chains = 0
    dir_base = '/Users/sashabrownsberger/Documents/Harvard/physics/'

    #OneD_posterior_file_name = 'H0_OmM_w_UncertaintyVsMuOffsetPrior_zLim_Full_' + '_'.join([str(id) for id in real_sequence_ids]) + '.pdf'
    OneD_posterior_sig_file_name = 'H0sigma_OmMwarea_sMuOffsetPrior_zLim_Full_' + '_'.join([str(id) for id in real_sequence_ids]) + '.pdf'
    OneD_posterior_centr_file_name = 'H0shift_OmMshift_wshift_VsMuOffsetPrior_zLim_Full_' + '_'.join([str(id) for id in real_sequence_ids]) + '.pdf'
    OmM_w0_contour_file_name = 'OmM_and_w0_contoursVsMuOffsetPrior_zLim_Full_' + '_'.join([str(id) for id in real_sequence_ids]) + '.pdf'
    save_dir = dir_base + 'stubbs/variableMuFits/plots/PosteriorWidths/'
    target_dir = 'stubbs/variableMuFits/mcmcResults/'
    #posterior_data_files = ['PosteriorFitsDifferentOffsets_zLim_Full_' + str(i) + '.txt' for i in range(sequence_start, sequence_end + 1)]
    mu_offsets = [0.001, 0.003, 0.005, 0.007, 0.009, 0.011, 0.013, 0.015, 0.017, 0.02, 0.04, 0.06, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2]
    mu_offsets = [0.001, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6]
    mu_offsets = [0.001, 0.012, 0.018, 0.03, 0.06, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6]
    mu_offsets = [0.001, 0.012, 0.03, 0.06, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6]
    #mu_offsets = [0.001, 0.025, 0.1, 0.25, 0.5]
    #mu_offsets = [0.002, 0.005, 0.008, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.25, 0.3] #Chris says info past this isn't useful
    #mu_offsets = [0.002, 0.005, 0.008, 0.01, 0.015, 0.02, 0.025, 0.03,  0.2, 0.3]
    #mu_offsets = [0.002, 0.008, 0.015, 0.02, 0.03,  0.2, 0.3]
    best_fit_mu_offsets = [0.03]
    n_offsets_to_lin_fit = 10
    bins = 50


    fancy_plot = 1
    if fancy_plot:
        # use LaTeX fonts in the plot
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')

    #mu_offsets = [0.002, 0.1, 0.5]
    #mu_offsets = [0.001, 0.04, 0.08, 0.14, 0.2]
    #mu_offsets = [0.1]
    OmM_lims, w0_lims = [[0.01, 0.52], [-1.7, -0.4]]
    OmM_lims, w0_lims = [[-0.6, 0.6], [-2.1, 1.1]]
    #OmM_lims, w0_lims = [[-3.5, -0.25], [-0.1, 0.75]]
    #surveys = ['CANDELS', 'CFA1',  'CFA2', 'CFA3K', 'CFA3S', 'CFA4p1', 'CFA4p2', 'CSP', 'DES', 'HST', 'KAIT', 'LOWZ', 'PS1MD', 'SDSS', 'SNAP', 'SNLS', 'SWIFT']
    surveys = ['SWIFT', 'ASASSN', 'CFA2', 'CFA1', 'KAITM', 'LOWZ', 'CFA4p2', 'KAIT', 'CFA3S', 'CSP', 'CFA3K', 'CFA4p1', 'PS1MD', 'SDSS', 'DES', 'SNLS', 'HST', 'SNAP', 'CANDELS']
    #surveys = ['CANDELS', 'CFA1',  'CFA2', 'CFA3K', 'CFA3S']
    #rand_mcmc_data_files_dict = {muSig:['GaussPriorWidth' + str(int(muSig * 1000)) + 'mMags_zLim_Full_MCMC_pantheon_plus_RAND' + 'H0_OmM_w0' + '_mu' + '_mu'.join(surveys) + '_NS10000_NC40_' + str(i) + '.txt' for i in sequence_ids] for muSig in mu_offsets}
    #real_mcmc_data_files_dict = {muSig:['GaussPriorWidth' + str(int(muSig * 1000)) + 'mMags_zLim_Full_MCMC_pantheon_plus_REAL' + 'H0_OmM_w0' + '_mu' + '_mu'.join(surveys) + '_NS10000_NC50_Real' + str(real_sequence_id) + '.txt' for real_sequence_id in real_sequence_ids] for muSig in mu_offsets}
    real_mcmc_data_files_dict = {muSig:['GaussPriorWidth' + str(int(muSig * 1000)) + 'mMags_zLim_Full_MCMC_pantheon_plus_REAL' + 'H0_OmM_w0' + '_mu' + '_mu'.join(surveys) + '_NS10000_NC50_' + str(real_sequence_id) + '.txt' for real_sequence_id in real_sequence_ids] for muSig in mu_offsets}
    f_sig, axarr_sig = plt.subplots(1,2, figsize = (9,4), sharex = True)
    f_centr, axarr_centr = plt.subplots(1,3, figsize = (13.5,4), sharex = True)
    labelsize = 11
    sup_title_size = 16
    n_burn_in = 20000
    x_scaling = 1 / 0.006


    #f, axarr = plt.subplots(2, 1, figsize = (5,3))
    #H0_ax = axarr
    H0_sig_ax = axarr_sig[0]
    OmM_w0_sig_ax = axarr_sig[1]
    #w0_sig_ax = axarr[2, 0]
    H0_mean_ax = axarr_centr[0]
    OmM_mean_ax = axarr_centr[1]
    w0_mean_ax = axarr_centr[2]
    #Om_ax = axarr[1]
    #w0_ax = axarr[2]
    #cfc.makePlotOfPosteriorWidths(rand_mcmc_data_files_dict, 0, ax = H0_ax, ref_param_vals_for_other_ax = [67.4, 0.5, r"$H_0$ tension (# of $\sigma$s)"], xlabel = '$\Delta \mu_S$ Gaussian prior width (mags)', ylabel = r'Increase in $H_0$' + '\n'  + r'posterior $\sigma$ (km/s/Mpc)', data_load_dir = target_dir, plot_color = 'k', plot_marker = '.', show_fig = 0, save_fig = 0)
    #cfc.makePlotOfPosteriorWidths(rand_mcmc_data_files_dict, 1, ax = Om_ax, xlabel = '$\Delta \mu_S$ Gaussian prior width (mags)', ylabel = r'Increase in $\Omega_M$' + '\n'  + r'posterior $\sigma$', data_load_dir = target_dir, plot_color = 'k', plot_marker = '.', show_fig = 0, save_fig = 0)
    #cfc.makePlotOfPosteriorWidths(rand_mcmc_data_files_dict, 2, ax = w0_ax, xlabel = '$\Delta \mu_S$ Gaussian prior width (mags)', ylabel = r'Increase in $w$' + '\n'  + r'posterior $\sigma$', data_load_dir = target_dir, plot_color = 'k', plot_marker = '.', show_fig = 0, save_fig = 0)
    #cfc.makePlotOfPosteriorWidths(rand_mcmc_data_files_dict, 0, ax = H0_ax, xlabel = '$\Delta \mu_S$ Gaussian prior width (mags)', ylabel = r'$H_0$' + r' posterior $\sigma$ (km/s/Mpc)', data_load_dir = target_dir, plot_color = 'k', plot_marker = '.', show_fig = 0, save_fig = 0, zero_uncertainties_at_zero_prior = 0, labelsize = labelsize)
    print ('real_mcmc_data_files_dict = ' + str(real_mcmc_data_files_dict))
    #We need to add the pre-inter survey offset uncertainties in quadrature from Planck and SH0ES.
    #xlabel = 'Inter-survey offset prior, \n width $\sigma_{P,S}$ (mags)'
    xlabel = 'Cross-survey photometric calibration \n uncertainty scaling, $N_{P, \sigma}$'
    sig_sup_title = 'The Additional Cosmic Parameter Measurement Uncertainty Introduced as the \n Uncertainty in the Pantheon+ Cross-Survey Photometric Zeropoint Increases'
    centr_sup_title = 'The Shift in the Cosmic Parameter Measurement  as the Uncertainty in the Pantheon+ Cross-Survey Photometric Zeropoint Increases'
    f_sig.suptitle(sig_sup_title, fontsize = sup_title_size)
    f_centr.suptitle(centr_sup_title, fontsize = sup_title_size)
    cfc.makePlotOfPosteriorWidths(real_mcmc_data_files_dict, 0, ax = H0_sig_ax, ref_param_vals_for_other_ax = [73.4, 67.27, np.sqrt(0.6 ** 2.0 + 1.22 ** 2.0), r"Planck vs SH0ES HT"], xlabel = xlabel, ylabel = r'Additional $H_0$' + '\n' + r'posterior $\sigma$ (km s$^{-1}$ Mpc$^{-1}$)', data_load_dir = target_dir, plot_color = 'k', plot_marker = 'x', show_fig = 0, save_fig = 0, zero_uncertainties_at_zero_prior = 1, labelsize = labelsize, show_all_chain_posteriors = show_all_chains, bins = bins,  n_points_to_lin_fit = n_offsets_to_lin_fit, n_burn_in = n_burn_in,
    #in_plot_lin_text = ['$\\frac{d \sigma_{H_0}}{d \\sigma_{P,S}} =$ ', r' $\frac{\mathrm{km}}{\mathrm{s} \ \mathrm{Mpc}} \frac{1}{\ 25\ \mathrm{mmag}}$'],
    in_plot_lin_text = ['Asymptotic $\sigma_{H_0}$ Degredation: ', r' $\frac{\mathrm{km}} {\mathrm{s} \ \mathrm{Mpc}} $'],
    #in_plot_lin_text = ['', ''],
    in_plot_loc = [0.14, 0.25], in_plot_text_size = labelsize, dir_base = dir_base, x_scaling = x_scaling)
    H0_sig_ax.text(0.5, 0.1, "Worse cross-survey calibration", ha="center", va="center", transform = H0_sig_ax.transAxes, fontsize = labelsize,  bbox=dict(boxstyle="rarrow,pad=0.3", fc="cyan", ec="b", lw=2))
    cfc.makePlotOfPosteriorWidths(real_mcmc_data_files_dict, 0, ax = H0_mean_ax, xlabel = xlabel, ylabel = r'Shift in $H_0$ centroid, $<H_0>$ (km s$^{-1}$ Mpc$^{-1}$)', data_load_dir = target_dir, plot_color = 'k', plot_marker = 'x', show_fig = 0, save_fig = 0, zero_uncertainties_at_zero_prior = 1, labelsize = labelsize, show_all_chain_posteriors = show_all_chains, plot_sigmas = 0, n_points_to_lin_fit = n_offsets_to_lin_fit, bins = bins,
    #in_plot_lin_text = ['$\\frac{d <H_0>}{d \\sigma_{P,S}} =$ ', r' $\frac{\mathrm{km}}{\mathrm{s} \ \mathrm{Mpc}} \frac{1}{25 \ \mathrm{ mmag}}$'],
    in_plot_lin_text = ['Asymptotic $<H_0>$ Shift: ', r' $\frac{\mathrm{km}}{\mathrm{s} \ \mathrm{Mpc}}$'],
    #in_plot_lin_text = ['', ''],
    in_plot_loc = [0.14, 0.25], in_plot_text_size = labelsize , dir_base = dir_base, n_burn_in = n_burn_in, x_scaling = x_scaling)
    H0_mean_ax.text(0.5, 0.1, "Worse cross-survey calibration", ha="center", va="center", transform = H0_mean_ax.transAxes, fontsize = labelsize,  bbox=dict(boxstyle="rarrow,pad=0.3", fc="cyan", ec="b", lw=2))
    #cfc.makePlotOfPosteriorWidths(rand_mcmc_data_files_dict, 1, ax = Om_ax, xlabel = '$\Delta \mu_S$ Gaussian prior width (mags)', ylabel = r'$\Omega_M$' + r' posterior $\sigma$', data_load_dir = target_dir, plot_color = 'k', plot_marker = '.', show_fig = 0, save_fig = 0, zero_uncertainties_at_zero_prior = 0, labelsize = labelsize)
    #cfc.makePlotOfPosteriorWidths(real_mcmc_data_files_dict, 1, ax = Om_ax, xlabel = '$\Delta \mu_S$ Gaussian prior width (mags)', ylabel = r'Increase in $\Omega_M$' + '\n' + r' posterior $\sigma$', data_load_dir = target_dir, plot_color = 'k', plot_marker = 'x', show_fig = 0, save_fig = 0, zero_uncertainties_at_zero_prior = 1, labelsize = labelsize, show_all_chain_posteriors = show_all_chains)
    cfc.makePlotOfPosteriorWidths(real_mcmc_data_files_dict, 1, ax = OmM_mean_ax, xlabel = xlabel, ylabel = r'Shift in $\Omega_M$ centroid, $<\Omega_M>$', data_load_dir = target_dir, plot_color = 'k', plot_marker = 'x', show_fig = 0, save_fig = 0, zero_uncertainties_at_zero_prior = 1, labelsize = labelsize, show_all_chain_posteriors = show_all_chains, plot_sigmas = 0, n_points_to_lin_fit = n_offsets_to_lin_fit, bins = bins,
    #in_plot_lin_text = ['$\\frac{d <\\Omega_M>}{d \\sigma_{P,S}} =$ ', r' $\frac{1}{25 \ \mathrm{ mmag}}$'],
    in_plot_lin_text = ['Asymptotic $<\\Omega_M>$ Shift: ', ' '],
    #in_plot_lin_text = ['', ''],
    in_plot_loc = [0.14, 0.25], in_plot_text_size = labelsize, dir_base = dir_base, n_burn_in = n_burn_in, x_scaling = x_scaling)
    OmM_mean_ax.text(0.5, 0.1, "Worse cross-survey calibration", ha="center", va="center", transform = OmM_mean_ax.transAxes, fontsize = labelsize,  bbox=dict(boxstyle="rarrow,pad=0.3", fc="cyan", ec="b", lw=2))
    #cfc.makePlotOfPosteriorWidths(rand_mcmc_data_files_dict, 2, ax = w0_ax, xlabel = '$\Delta \mu_S$ Gaussian prior width (mags)', ylabel = r'$w$' + r' posterior $\sigma$', data_load_dir = target_dir, plot_color = 'k', plot_marker = '.', show_fig = 0, save_fig = 0, zero_uncertainties_at_zero_prior = 0, labelsize = labelsize)
    #cfc.makePlotOfPosteriorWidths(real_mcmc_data_files_dict, 2, ax = w0_sig_ax, xlabel = '$\Delta \mu_S$ Gaussian prior width (mags)', ylabel = r'Increase in $w$ posterior $\sigma$', data_load_dir = target_dir, plot_color = 'k', plot_marker = 'x', show_fig = 0, save_fig = 0, zero_uncertainties_at_zero_prior = 1, labelsize = labelsize, show_all_chain_posteriors = show_all_chains, n_points_to_lin_fit = n_offsets_to_lin_fit, in_plot_lin_text = ['$\\frac{d \sigma_w}{d \\sigma_{P,S}} = $', r' mmag$^{-1}$'])
    cfc.makePlotOfPosteriorWidths(real_mcmc_data_files_dict, 2, ax = w0_mean_ax, xlabel = xlabel, ylabel = r'Shift in $w$ centroid, $<w>$', data_load_dir = target_dir, plot_color = 'k', plot_marker = 'x', show_fig = 0, save_fig = 0, zero_uncertainties_at_zero_prior = 1, labelsize = labelsize, show_all_chain_posteriors = show_all_chains, plot_sigmas = 0, n_points_to_lin_fit = n_offsets_to_lin_fit, bins = bins,
    #in_plot_lin_text = ['$\\frac{d <w>}{d \\sigma_{P,S}} =$ ', r' $\frac{1}{25 \ \mathrm{mmag}}$'],
    in_plot_lin_text = ['Asymptotic $<w>$ Shift: ', ' '],
    #in_plot_lin_text = ['', ''],
    in_plot_loc = [0.14, 0.75], in_plot_text_size = labelsize, dir_base = dir_base, n_burn_in = n_burn_in, x_scaling = x_scaling)
    w0_mean_ax.text(0.5, 0.9, "Worse cross-survey calibration", ha="center", va="center", transform = w0_mean_ax.transAxes, fontsize = labelsize,  bbox=dict(boxstyle="rarrow,pad=0.3", fc="cyan", ec="b", lw=2))
    #plt.subplots_adjust (left = 0.20, right = 0.85)
    #plt.savefig(save_dir + OneD_posterior_file_name)
    #plt.close('all')
    full_mcmc_data_files_list = [real_mcmc_data_files_dict[mu_offset] for mu_offset in mu_offsets]
    area_contour_num = 0
    OmM_w0_mesh_axarr = cfc.makePlotOfPosteriorContours(full_mcmc_data_files_list, [2, 3], OmM_lims, w0_lims, xlabel = '$\Delta \Omega_M$', ylabel = r'$\Delta w$', plot_text = [r'$\Delta \mu_S$ priors of ' +  str(int(mu_offset * 1000)) + ' mmag' + (' ' if mu_offset == 0.001 else 's') for mu_offset in mu_offsets], contour_area_ax = OmM_w0_sig_ax, area_contour_num = area_contour_num, contour_area_x_vals = mu_offsets, area_plot_color = 'k', area_plot_marker = 'x', bins = bins, n_points_to_lin_fit = n_offsets_to_lin_fit,
    in_plot_lin_text = [r'Asymptotic $\sigma_{w \times \Omega_M}$ Degredation: ', '\%'],
    #in_plot_lin_text = ['', ''],
    in_plot_text_size = labelsize, in_plot_loc = [0.14, 0.25], dir_base = dir_base, n_burn_in = n_burn_in, x_scaling = x_scaling)
    OmM_w0_sig_ax.text(0.5, 0.1, "Worse cross-survey calibration", ha="center", va="center", transform = OmM_w0_sig_ax.transAxes, fontsize = labelsize,  bbox=dict(boxstyle="rarrow,pad=0.3", fc="cyan", ec="b", lw=2))
    OmM_w0_sig_ax.set_xlabel(xlabel, fontsize = labelsize)
    OmM_w0_sig_ax.set_ylabel(r'Additional $\Omega_M \times w$ ' + str(area_contour_num + 1) + r'$\sigma$' + '\n' +  'posterior area (\%)', fontsize = labelsize)
    H0_sig_ax.axvline(best_fit_mu_offsets[0] * x_scaling, c = 'k', linestyle = '--')
    H0_mean_ax.axvline(best_fit_mu_offsets[0] * x_scaling, c = 'k', linestyle = '--')
    OmM_mean_ax.axvline(best_fit_mu_offsets[0] * x_scaling, c = 'k', linestyle = '--')
    #w0_sig_ax.axvline(best_fit_mu_offsets[0], c = 'k', linestyle = '--')
    w0_mean_ax.axvline(best_fit_mu_offsets[0] * x_scaling, c = 'k', linestyle = '--')
    OmM_w0_sig_ax.axvline(best_fit_mu_offsets[0] * x_scaling, c = 'k', linestyle = '--')
    #OmM_w0_mesh_axarr = cfc.makePlotOfPosteriorContours(full_mcmc_data_files, [2, 3], OmM_lims, w0_lims, plot_text = None)
    #f_sig.subplots_adjust (wspace = 1, hspace = 0.1,  bottom = 0.1, left = 0.1, right = 0.70)
    plt.savefig(save_dir + OmM_w0_contour_file_name)
    plt.tight_layout()
    plt.figure(f_sig.number)
    plt.tight_layout()
    #f_centr.subplots_adjust (wspace = 0.25, hspace = 0.1, bottom = 0.1, right= 2.0)
    plt.savefig(save_dir + OneD_posterior_sig_file_name)
    plt.figure(f_centr.number)
    plt.tight_layout()
    plt.savefig(save_dir + OneD_posterior_centr_file_name)
    plt.close('all')

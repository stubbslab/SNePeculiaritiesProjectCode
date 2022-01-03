import CosmicFitterClass as cfc
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import cantrips as can
import scipy
import sys
import matplotlib.patches as mpatches

def param_fit_funct_full(zs, H0, OmM, w, OmL, OmR, Om0):
    zs = list(zs)
    w_of_funct = lambda w_zs: fitter_to_plot.w_of_funct(w_zs, [w])
    funct_mus = cfc.getMusForCosmology(zs, [H0, OmM, OmL, OmR, Om0], lambda w_zs: fitter_to_plot.w_of_funct(w_zs, [w]) )
    #print ('zs = ' + str(zs))
    return funct_mus

if __name__ == "__main__":
    overplot_steps = 1000
    overplot_n_chains = 2
    z_cut_lims = [0.00, 2.5]
    z_lims = [0.005, 2.5]
    surveys = ['SWIFT', 'ASASSN', 'CFA1', 'CFA2', 'LOWZ', 'KAITM', 'CFA4p2', 'KAIT', 'CFA3S', 'CSP', 'CFA3K', 'CFA4p1', 'PS1MD', 'SDSS', 'DES', 'SNLS', 'HST', 'SNAP', 'CANDELS']
    ticklabelsize = 18
    labelsize = 24
    point_transparency = 0.5
    fitter_to_plot = cfc.CosmicFitter(w_of_funct_str = 'w0', randomize_sn = 0, sn_data_type = 'pantheon_plus', params_to_overplot = ['H0'], overplot_mcmc_steps = overplot_steps, overplot_mcmc_chains = overplot_n_chains)
    fitter_to_plot.updateUsedSN(z_lims = z_cut_lims, surveys_to_include = surveys )

    fancy_plot = 1
    if fancy_plot:
        # use LaTeX fonts in the plot
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')

    #sorted_zs, sorted_mus, sorted_muErrs, colors = [fitter_to_plot.sorted_zs, fitter_to_plot.sorted_mus, fitter_to_plot.sorted_muErrs, fitter_to_plot.sorted_colors]
    save_plot_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/variableMuFits/plots/'
    plot_file_name = 'PantheonPlus_realMuResidsVsZ_boxes.pdf'

    #canon_mus = fitter_to_plot.getMusForCosmology(fitter_to_plot.nonCeph_zs, [fitter_to_plot.H0, fitter_to_plot.OmegaM, fitter_to_plot.OmegaLambda, fitter_to_plot.OmegaR , fitter_to_plot.Omega0], wOfFunction = lambda zs: fitter_to_plot.w_of_funct(zs, fitter_to_plot.default_w_params), )
    wOfFunction = lambda zs: fitter_to_plot.w_of_funct(zs, fitter_to_plot.default_w_params)
    param_fit_funct = lambda zs, H0, OmM, w: param_fit_funct_full(zs, H0, OmM, w, fitter_to_plot.OmegaLambda, fitter_to_plot.OmegaR , fitter_to_plot.Omega0)
    zs_to_fit = fitter_to_plot.nonCeph_zs
    mus_to_fit = fitter_to_plot.nonCeph_mus
    muErrs_to_fit = fitter_to_plot.nonCeph_mus
    canon_mus = cfc.getMusForCosmology(fitter_to_plot.nonCeph_zs, [fitter_to_plot.H0, fitter_to_plot.OmegaM, fitter_to_plot.OmegaLambda, fitter_to_plot.OmegaR , fitter_to_plot.Omega0], wOfFunction )
    init_guess = [fitter_to_plot.H0, fitter_to_plot.OmegaM, fitter_to_plot.default_w_params[0]]
    start_mus = param_fit_funct(zs_to_fit, *init_guess)
    #print ('[len(zs_to_fit), len(mus_to_fit), len(canon_mus), len(start_mus)] = ' + str([len(zs_to_fit), len(mus_to_fit), len(canon_mus), len(start_mus)] ))
    best_fit_params = scipy.optimize.curve_fit(param_fit_funct, zs_to_fit, mus_to_fit, p0 = init_guess, sigma = muErrs_to_fit)[0]
    #print ('best_fit_params = ' + str(best_fit_params))
    wOfFunction = lambda zs: fitter_to_plot.w_of_funct(zs, [best_fit_params[2]])
    canon_mus = cfc.getMusForCosmology(fitter_to_plot.nonCeph_zs, [best_fit_params[0], best_fit_params[1], fitter_to_plot.OmegaLambda, fitter_to_plot.OmegaR , fitter_to_plot.Omega0], wOfFunction )
    fill_zs = np.logspace(-3, np.log10(2.5), 221).tolist()
    #fill_mus = fitter_to_plot.getMusForCosmology(fill_zs, [fitter_to_plot.H0, fitter_to_plot.OmegaM, fitter_to_plot.OmegaLambda, fitter_to_plot.OmegaR , fitter_to_plot.Omega0], wOfFunction = lambda zs: fitter_to_plot.w_of_funct(zs, fitter_to_plot.default_w_params), )
    fill_mus = cfc.getMusForCosmology(fill_zs, [fitter_to_plot.H0, fitter_to_plot.OmegaM, fitter_to_plot.OmegaLambda, fitter_to_plot.OmegaR , fitter_to_plot.Omega0], wOfFunction )
    #print ('[len(fill_zs), len(fill_mus)] = ' + str([len(fill_zs), len(fill_mus)] ))
    read_in_surveys = fitter_to_plot.surveys
    surveys = ['SWIFT', 'ASASSN', 'CFA1', 'CFA2', 'LOWZ', 'KAITM', 'CFA4p2', 'KAIT', 'CFA3S', 'CSP', 'CFA3K', 'CFA4p1', 'PS1MD', 'SDSS', 'DES', 'SNLS', 'HST', 'SNAP', 'CANDELS']
    for survey in read_in_surveys:
        if not(survey in surveys): print ('Survey ' + str(survey) + ' not in expected list of surveys')
    for survey in surveys:
        if not(survey in read_in_surveys): print ('Survey ' + str(survey) + ' not expected, but not found in fitter')
    surveys = [survey for survey in surveys if survey in read_in_surveys]
    #print ('surveys = '  +str(surveys) )
    f, axarr = plt.subplots(1,1, figsize = [11, 7], sharex = True)
    scats = []
    #ref_cosmo_plot = axarr[0].plot([-0.01] + fill_zs, [0.0] + fill_mus.tolist(), c = 'k')[0]
    bins = np.logspace(np.log10(z_lims[0]),np.log10(z_lims[1]), 25)
    all_hist_handles = []
    #full_hist = axarr[1].hist(fitter_to_plot.nonCeph_zs, bins = bins, edgecolor = 'k', fill = False, histtype = 'step', linestyle = '--')
    survey_boxes = []
    for survey in surveys:
        zs = [fitter_to_plot.nonCeph_zs[i] for i in range(len(fitter_to_plot.nonCeph_zs)) if fitter_to_plot.nonCeph_surveys[i] == survey]
        z_range = [min(zs), max(zs)]
        muDiffs = [fitter_to_plot.nonCeph_mus[i] - canon_mus[i] for i in range(len(fitter_to_plot.nonCeph_mus)) if fitter_to_plot.nonCeph_surveys[i] == survey]
        muErrs = [fitter_to_plot.nonCeph_muErrs[i] for i in range(len(fitter_to_plot.nonCeph_muErrs)) if fitter_to_plot.nonCeph_surveys[i] == survey]
        muWeights = [err ** -2.0 for err in muErrs]
        mu_wmean = np.average(muDiffs, weights = muWeights)
        #mu_wstd =  np.sqrt(np.sum((np.array(muWeights) * np.array(muErrs))** 2.0)) / np.sum(muWeights)
        mu_wstd =  np.sqrt(np.sum((np.array(muWeights) * np.array(muErrs))** 2.0)) / np.sum(muWeights)
        colors = [fitter_to_plot.nonCeph_plotColors[i] for i in range(len(fitter_to_plot.nonCeph_plotColors)) if fitter_to_plot.nonCeph_surveys[i] == survey]
        #scats = scats + [axarr[0].scatter(zs, muDiffs, c = colors, marker = '.', alpha = point_transparency)]
        #axarr[0].errorbar(zs,  muDiffs, yerr = muErrs, color = colors, fmt = 'none', alpha = point_transparency)
        survey_box = Rectangle((z_range[0], mu_wmean - mu_wstd), z_range[1] - z_range[0], 2.0 * mu_wstd, ec = colors[0], alpha = 1, facecolor='none')
        survey_boxes = survey_boxes + [survey_box]
        axarr.add_patch(survey_box)
        #new_hist = axarr[1].hist(zs, bins = bins, edgecolor = colors[0], histtype = 'step', alpha = 0.75)
        all_hist_handles = all_hist_handles + [Rectangle((0,0),1,1,facecolor='none',ec=colors[0], alpha = 0.75)]
        #axarr[0].plot([-1.0] + fill_zs + [2.1], [weighted_mean] + [weighted_mean for z in fill_zs] + [weighted_mean], c = colors[0])
        #print (survey + ': \Delta \mu = ' + str(can.round_to_n(mu_wmean, 3)) + ' \pm ' + str(can.round_to_n(mu_wstd, 3)))
        #if survey == 'CFA1':
        #    print ('[muDiffs, muWeights] = ' + str([muDiffs, muWeights]))
    all_hist_handles = all_hist_handles + [Rectangle((0,0),1,1,color='white',ec='k', linestyle = '--')]
    #axarr[1].legend(all_hist_handles, surveys + ['All'], ncol = 9, fontsize = 8 )

    axarr.set_ylabel(r'$\mu$ (mag)', fontsize = labelsize )
    ref_cosmo_plot = axarr.plot([-1.0] + fill_zs + [2.1], [0.0] + [0.0 for z in fill_zs] + [0.0], c = 'k')[0]
    axarr.legend(survey_boxes + [ref_cosmo_plot], surveys + ['SH0ES+Pantheon+ ' +  r'best-fit $\Lambda$CDM'], ncol = 5, fontsize = 11, loc = 'lower left' )
    axarr.tick_params(labelsize=ticklabelsize)
    #axarr[1].scatter(fitter_to_plot.nonCeph_zs, np.array(fitter_to_plot.nonCeph_mus) - np.array(canon_mus), c = fitter_to_plot.nonCeph_plotColors, marker = '.')
    #axarr[1].errorbar(fitter_to_plot.nonCeph_zs,  np.array(fitter_to_plot.nonCeph_mus) - np.array(canon_mus), yerr = fitter_to_plot.nonCeph_muErrs, color = fitter_to_plot.nonCeph_plotColors, fmt = 'none')
    #axarr[1].set_xlabel(r'$z$', fontsize = labelsize)
    axarr.set_xlabel(r'$z$', fontsize = labelsize)
    axarr.set_ylabel(r'$\mu-\mu_{\mathrm{model}}$ (mag)', fontsize = labelsize)
    axarr.set_xlim(z_lims)
    axarr.set_ylim([-0.25, 0.25])
    #axarr[1].set_xlim(z_lims)
    #axarr[1].set_ylim(-1.05, 1.05)
    #axarr[0].set_ylim(-1.2, 1.2)
    #axarr[1].set_ylim(0.5, 1000)
    #axarr[0].set_ylim ([min(fitter_to_plot.nonCeph_mus) - (max(fitter_to_plot.nonCeph_mus) - min(fitter_to_plot.nonCeph_mus)) * 0.05, max(fitter_to_plot.nonCeph_mus) + (max(fitter_to_plot.nonCeph_mus) - min(fitter_to_plot.nonCeph_mus)) * 0.05])
    #axarr[1].set_ylabel(r'Number of SNe Ia', fontsize = labelsize )
    #axarr[1].set_ylabel(r'$\mu-\mu_{model}$ (mag)', fontsize = labelsize)
    #axarr[1].tick_params(labelsize=ticklabelsize)
    axarr.set_xscale('log')
    #axarr[1].set_xscale('log')
    #axarr[1].set_yscale('log')
    #plt.text(0.02, 0.11, r"$\}$", fontsize=50, rotation = 90, horizontalalignment = 'center',  stretch = 1)
    axarr.annotate(r'Low-$z$ surveys have SNe in both' + '\n calibrator and Hubble Flow sample.', xy = (0.022, 0.17), fontsize = 14, ha='center', va='bottom',
            bbox=dict(boxstyle='square', fc='white'),
            arrowprops=dict(arrowstyle='-[, widthB=10.0, lengthB=1.5', lw=2.0))
    CFA1_arrow_coords = [[0.13, -0.11], [0.15, -0.11]]
    arrow = mpatches.Arrow(*CFA1_arrow_coords[0], -0.07, 0.0, width = 0.1, facecolor = 'grey', edgecolor = 'k')
    axarr.add_patch(arrow)
    #axarr.text(0.12, -0.17, r'Because CFA1 is notably' + '\n' + 'offset from the best' + '\n'+ r'Pantheon+ $\Lambda$CDM fit, $\Delta \mu_{CFA1}$' + '\n' + 'is statistically distinct' + '\n' + 'from $0$ (see Figure 2).', fontsize = 14)
    axarr.annotate(r'Because the CFA1 data is the' + '\n' + 'most offset from the best fit' + '\n'+ r'Pantheon+ $\Lambda$CDM fit, $\Delta \mu_{CFA1}$' + '\n' + r'is the $\Delta \mu_S$ that is most offset' + '\n'  + 'from $0$ (see Table 1).', fontsize = 14, xy = (0.14, -0.11), ha='left', va='center',
            bbox=dict(boxstyle='square', fc='white'))
    #axarr.text(0.12, 0.9, r'Inter-survey offset $\Delta \mu_S$ vertically' + '\n' +  'shifts all SNe in survey $S$', fontsize = 14, horizontalalignment = 'center')
    up_arrow = mpatches.Arrow(0.25, 0.09, 0.0, 0.05, width = 0.1, facecolor = 'k', edgecolor = 'k')
    axarr.add_patch(up_arrow)
    down_arrow = mpatches.Arrow(0.25, 0.09, 0.0, -0.05, width = 0.1, facecolor = 'k', edgecolor = 'k')
    axarr.add_patch(down_arrow)
    axarr.annotate(r'Inter-survey offset $\Delta \mu_S$ vertically' + '\n' +  'shifts all SNe in survey $S$.', xy = (0.25, 0.09), fontsize = 14, ha='center', va='center', bbox=dict(boxstyle='square', fc='white'))
    f.subplots_adjust(hspace=0.05)
    plt.tight_layout()
    plt.savefig(save_plot_dir + plot_file_name)
    plt.show()

import CosmicFitterClass as cfc
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    overplot_steps = 1000
    overplot_n_chains = 2
    z_lims = [0.005, 2.1]
    ticklabelsize = 12
    labelsize = 18
    fitter_to_plot = cfc.CosmicFitter(w_of_funct_str = 'w0', randomize_sn = 1, sn_data_type = 'pantheon_plus', params_to_overplot = ['H0'], overplot_mcmc_steps = overplot_steps, overplot_mcmc_chains = overplot_n_chains)
    #sorted_zs, sorted_mus, sorted_muErrs, colors = [fitter_to_plot.sorted_zs, fitter_to_plot.sorted_mus, fitter_to_plot.sorted_muErrs, fitter_to_plot.sorted_colors]
    save_plot_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/variableMuFits/plots/'
    plot_file_name = 'PantheonPlus_randomizeMuVsZ.pdf'

    #canon_mus = fitter_to_plot.getMusForCosmology(fitter_to_plot.nonCeph_zs, [fitter_to_plot.H0, fitter_to_plot.OmegaM, fitter_to_plot.OmegaLambda, fitter_to_plot.OmegaR , fitter_to_plot.Omega0], wOfFunction = lambda zs: fitter_to_plot.w_of_funct(zs, fitter_to_plot.default_w_params), )
    wOfFunction = lambda zs: fitter_to_plot.w_of_funct(zs, fitter_to_plot.default_w_params)
    canon_mus = cfc.getMusForCosmology(fitter_to_plot.nonCeph_zs, [fitter_to_plot.H0, fitter_to_plot.OmegaM, fitter_to_plot.OmegaLambda, fitter_to_plot.OmegaR , fitter_to_plot.Omega0], wOfFunction )
    fill_zs = np.logspace(-3, np.log10(2.5), 221).tolist()
    #fill_mus = fitter_to_plot.getMusForCosmology(fill_zs, [fitter_to_plot.H0, fitter_to_plot.OmegaM, fitter_to_plot.OmegaLambda, fitter_to_plot.OmegaR , fitter_to_plot.Omega0], wOfFunction = lambda zs: fitter_to_plot.w_of_funct(zs, fitter_to_plot.default_w_params), )
    fill_mus = cfc.getMusForCosmology(fill_zs, [fitter_to_plot.H0, fitter_to_plot.OmegaM, fitter_to_plot.OmegaLambda, fitter_to_plot.OmegaR , fitter_to_plot.Omega0], wOfFunction )
    print ('[len(fill_zs), len(fill_mus)] = ' + str([len(fill_zs), len(fill_mus)] ))
    surveys = np.unique(fitter_to_plot.nonCeph_surveys).tolist()
    print ('surveys = '  +str(surveys) )
    f, axarr = plt.subplots(2,1, figsize = [10, 6], sharex = True)
    scats = []
    ref_cosmo_plot = axarr[0].plot([-0.01] + fill_zs, [0.0] + fill_mus.tolist(), c = 'k')[0]
    for survey in surveys:
        zs = [fitter_to_plot.nonCeph_zs[i] for i in range(len(fitter_to_plot.nonCeph_zs)) if fitter_to_plot.nonCeph_surveys[i] == survey]
        mus = [fitter_to_plot.nonCeph_mus[i] for i in range(len(fitter_to_plot.nonCeph_mus)) if fitter_to_plot.nonCeph_surveys[i] == survey]
        muErrs = [fitter_to_plot.nonCeph_muErrs[i] for i in range(len(fitter_to_plot.nonCeph_muErrs)) if fitter_to_plot.nonCeph_surveys[i] == survey]
        colors = [fitter_to_plot.nonCeph_plotColors[i] for i in range(len(fitter_to_plot.nonCeph_plotColors)) if fitter_to_plot.nonCeph_surveys[i] == survey]
        scats = scats + [axarr[0].scatter(zs, mus, c = colors, marker = '.')]
        axarr[0].errorbar(zs,  mus, yerr = muErrs, color = colors, fmt = 'none')

    axarr[0].set_ylabel(r'$\mu$ (mag)', fontsize = labelsize )
    axarr[0].legend(scats + [ref_cosmo_plot], surveys + ['Pantheon+ ' +  r'$\Lambda$CDM'], ncol = 3, fontsize = 9)
    axarr[0].tick_params(labelsize=ticklabelsize)
    axarr[1].plot([-1.0] + fill_zs + [2.1], [0.0] + [0.0 for z in fill_zs] + [0.0], c = 'k')
    axarr[1].scatter(fitter_to_plot.nonCeph_zs, np.array(fitter_to_plot.nonCeph_mus) - np.array(canon_mus), c = fitter_to_plot.nonCeph_plotColors, marker = '.')
    axarr[1].errorbar(fitter_to_plot.nonCeph_zs,  np.array(fitter_to_plot.nonCeph_mus) - np.array(canon_mus), yerr = fitter_to_plot.nonCeph_muErrs, color = fitter_to_plot.nonCeph_plotColors, fmt = 'none')
    axarr[1].set_xlabel(r'$z$', fontsize = labelsize)
    axarr[1].set_ylabel(r'$\mu-\mu_{model}$ (mag)', fontsize = labelsize)
    axarr[0].set_xlim(z_lims)
    axarr[1].set_xlim(z_lims)
    axarr[1].set_ylim(-1.05, 1.05)
    axarr[0].set_ylim ([min(fitter_to_plot.nonCeph_mus) - (max(fitter_to_plot.nonCeph_mus) - min(fitter_to_plot.nonCeph_mus)) * 0.05, max(fitter_to_plot.nonCeph_mus) + (max(fitter_to_plot.nonCeph_mus) - min(fitter_to_plot.nonCeph_mus)) * 0.05])
    #axarr[1].set_ylabel(r'$\Delta \mu$ (mag)', fontsize = labelsize )
    axarr[1].set_ylabel(r'$\mu-\mu_{model}$ (mag)', fontsize = labelsize)
    axarr[1].tick_params(labelsize=ticklabelsize)
    axarr[0].set_xscale('log')
    axarr[1].set_xscale('log')
    f.subplots_adjust(hspace=0.05)
    plt.tight_layout()
    plt.savefig(save_plot_dir + plot_file_name)
    plt.show()

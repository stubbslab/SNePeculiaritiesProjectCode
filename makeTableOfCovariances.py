import cantrips as can
import CosmicFitterClass as cfc
import matplotlib.pyplot as plt
import binData as bd
import numpy as np
import random
import sys
import scipy
import SNDataArchive as snda


if __name__=="__main__":
    """
    Makes a table formatted specifically to paste into the latex table and prints it to the command line.
    $ python makeTableOfCovariances.py 0.2 1
    #  missing chains - 15
    """
    cl_args = sys.argv[1:]
    test_mu_offset, sequence_ids = [cl_args[0], cl_args[1:]]
    test_mu_offset, sequence_ids = [float(test_mu_offset), [id for id in sequence_ids]]
    print ('[test_mu_offset, sequence_ids] = ' + str([test_mu_offset, sequence_ids]))

    dir_base = '/Users/sashabrownsberger/Documents/Harvard/physics/'

    fitter_for_sn = cfc.CosmicFitter(w_of_funct_str = 'w0', randomize_sn = 0, sn_data_type = 'pantheon_plus', params_to_overplot = None, dir_base = dir_base)

    print('fitter_for_sn.surveys  = ' + str(fitter_for_sn.surveys))
    n_sig_fig_errs = 2

    save_dir = dir_base + 'stubbs/variableMuFits/plots/PosteriorWidths/'
    target_dir = dir_base + 'stubbs/variableMuFits/mcmcResults/'
    surveys = ['SWIFT', 'ASASSN', 'CFA2', 'CFA1', 'KAITM', 'LOWZ', 'CFA4p2', 'KAIT', 'CFA3S', 'CSP', 'CFA3K', 'CFA4p1', 'PS1MD', 'SDSS', 'DES', 'SNLS', 'HST', 'SNAP', 'CANDELS']
    covariance_files = ['GaussPriorWidth' + str(int(test_mu_offset * 1000)) + 'mMags_zLim_Full_PosteriorStats_pantheon_plus_REALH0_OmM_w0_mu' + '_mu'.join(surveys) + '_NS10000_NC50_' + str(i) + '.txt' for i in sequence_ids]
    raw_mcmc_files = ['GaussPriorWidth' + str(int(test_mu_offset * 1000)) + 'mMags_zLim_Full_MCMC_pantheon_plus_REALH0_OmM_w0_mu' + '_mu'.join(surveys) + '_NS10000_NC50_' + str(i) + '.txt' for i in sequence_ids]
    raw_mcmc_data = can.readInColumnsToList(target_dir + raw_mcmc_files[0], n_ignore = 1, delimiter = ', ', convert_to_float = 1, n_ignore_end = 1)
    median_survey_offsets = [np.median(col) for col in raw_mcmc_data[3:]]
    std_survey_offsets = [np.std(col) for col in raw_mcmc_data[3:]]
    print ('[len(median_survey_offsets), len(surveys)] = ' + str([len(median_survey_offsets), len(surveys)]))
    covariance_matrices = [ can.readInColumnsToList(target_dir + file, n_ignore = 0, delimiter = ', ', convert_to_float = 0) for file in covariance_files ][0]
    covariance_matrix_labels = [ col[0] for col in covariance_matrices ]
    covariance_matrices = [can.readInColumnsToList(target_dir + file, n_ignore = 1, delimiter = ', ', convert_to_float = 1) for file in covariance_files]
    #print ('covariance_matrices[0] = ' + str(covariance_matrices[0]))
    #print ('covariance_matrices[-1] = ' + str(covariance_matrices[-1]))
    mean_covariance_matrix = np.mean(covariance_matrices, axis = 0)
    std_covariance_matrix = np.std(covariance_matrices, axis = 0)

    #print ('std_covariance_matrix = ' + str(std_covariance_matrix))

    #extra_print_elems = [   '$(\\cdot) \\times H_0$ (km s$^{-1}$ Mpc$^{-1}$) & ',         '$(\\cdot) \\times \Omega_M$ & ',                     '$(\\cdot) \\times w$ & '
    truncated_labels = [label[2:] for label in covariance_matrix_labels]
    surveys = truncated_labels[2:]
    print ('surveys = ' + str(surveys))

    cepheid_indeces = fitter_for_sn.cepheid_indeces
    n_cepheids_by_survey = {survey:0 for survey in surveys}
    n_HF_by_survey = {survey:0 for survey in surveys}
    for k in range(len(fitter_for_sn.sorted_surveys)):
        survey = fitter_for_sn.sorted_surveys[k]
        if k in cepheid_indeces:
            n_cepheids_by_survey[survey] = n_cepheids_by_survey[survey]  + 1
        else:
            n_HF_by_survey[survey] = n_HF_by_survey[survey]  + 1
    print ('n_cepheids_by_survey = ' + str(n_cepheids_by_survey))
    print ('n_HF_by_survey = ' + str(n_HF_by_survey))


    extra_print_elems = ['$\Delta \mu_{' + '{}'.format(label) + '}$' for label in truncated_labels]
    print ('covariance_matrix_labels = ' + str(covariance_matrix_labels))
    print ('extra_print_elems = ' + str(extra_print_elems))

    #top_row = '$\\Delta \\mu_S$ & $(\\cdot) \\times H_0$ (km s$^{-1}$ Mpc$^{-1}$)' + ' & ' + '$(\\cdot) \\times \Omega_M$' + ' & ' + '$(\\cdot) \\times w$' + ' & ' + '$N_{cal}$' + ' & ' + '$N_{HF}$' + ' \\\\'
    top_row = '$ \\Delta \\mu_S $ & $< \\Delta \\mu_S>$(mmag) & $\\frac{d H_0} { d (\\cdot)} ( \\frac{\\textrm{km}}{\\textrm{s Mpc 100 mmag }} )$  & $\\frac{d \\Omega_M } { d (\\cdot)} (\\frac{1}{\\textrm{100 mmag}})$ & $\\frac{d w } { d (\\cdot)}  (\\frac{1}{\\textrm{100 mmag}})$ & $N_{cal}$ & $N_{HF}$ & $N_{cal} / N_{HF}$ \\\\'
    print (top_row)
    print ('\\hline')
    x_vals_to_plot = []
    y_vals_to_plot = []
    colors = []
    legends = []
    data_points_text = []
    color_dict = snda.SNDataArchive().survey_color_map
    for j in range(3, len(mean_covariance_matrix[0])):
        survey = truncated_labels[j]
        mean_row = mean_covariance_matrix[j]
        offset_variance = mean_row[j]
        #linear_fits = [scipy.optimize.curve_fit(lambda xs, slope, intercept: slope * xs + intercept, np.array(raw_mcmc_data[j]), np.array(raw_mcmc_data[i]), p0 = [0, 0]) for i in range(3)]
        n_bootstraps = 100
        boot_linear_fits = [[] for k in range(n_bootstraps )]
        for k in range(n_bootstraps):
            boot_indeces = [random.randrange(len(raw_mcmc_data[0])) for k in range(len(raw_mcmc_data[0]))]
            boot_xs = [raw_mcmc_data[j][index] for index in boot_indeces]
            boot_ys = [[raw_mcmc_data[i][index] for index in boot_indeces] for i in range(3)]
            linear_fits = [np.polyfit(np.array(boot_xs), np.array(boot_ys[i]), 1)[0] for i in range(3)]
            boot_linear_fits[k] = linear_fits
        boot_mean_linear_fits = [np.mean([linear_fits[i] for linear_fits in boot_linear_fits]) for i in range(3)]
        boot_std_linear_fits = [np.std([linear_fits[i] for linear_fits in boot_linear_fits]) for i in range(3)]
        true_linear_fits = [np.polyfit(np.array(raw_mcmc_data[j]), np.array(raw_mcmc_data[i]), 1) for i in range(3)]
        """
        print ('true_linear_fits[0] = ' + str(true_linear_fits[0]))
        print ('boot_mean_linear_fits[0] = ' + str(boot_mean_linear_fits[0]))
        print ('boot_std_linear_fits[0] = ' + str(boot_std_linear_fits[0]))
        plt.scatter(np.array(raw_mcmc_data[j]), np.array(raw_mcmc_data[0]), marker = '.')
        plot_xs = [np.min(raw_mcmc_data[j]), np.max(raw_mcmc_data[j])]
        plt.plot(plot_xs, np.poly1d(true_linear_fits[0])(plot_xs), c = 'r')
        plt.show()
        print ('true_linear_fits = ' + str(true_linear_fits))
        print ('boot_mean_linear_fits = ' + str(boot_mean_linear_fits))
        print ('boot_std_linear_fits = ' + str(boot_std_linear_fits))
        """
        #print ('offset_variance = ' + str(offset_variance))
        #std_row = std_covariance_matrix[i]
        #print ('std_row = ' + str(std_row))
        #print ([[mean_row[j], std_row[j]] for j in range(len(mean_row))])
        #print ([[np.log10(abs(mean_row[j])),np.log10(abs( std_row[j]))] for j in range(len(mean_row))])
        #print ([[int(np.log10(abs(mean_row[j]))), int(np.log10(abs( std_row[j])))] for j in range(len(mean_row))])
        #round_n_sig_figs = [int(np.log10(abs(mean_row[j]))) - int(np.log10(abs(std_row[j]))) for j in range(len(mean_row)) ]
        round_n_sig_figs = [1 for i in range(3)]
        #print ('round_n_sig_figs = ' + str(round_n_sig_figs))
        #row_to_print = extra_print_elems[i] + '$' + ('$ & $'.join([str(can.round_to_n(mean_row[j] * (1000 if j > 2 else 1), n_sig_fig_errs + round_n_sig_figs[j])) + r' \pm ' + str(can.round_to_n(std_row[j] * (1000 if j > 2 else 1), n_sig_fig_errs) ) + ('$ \\\\ \n $ ' if j % 5 == 5-1 else '') for j in range(len(mean_row))])) [0:-2] + '\n'
        i = 2
        x_vals_to_plot  = x_vals_to_plot + [n_cepheids_by_survey[survey] / n_HF_by_survey[survey]]
        y_vals_to_plot  = y_vals_to_plot + [true_linear_fits[0][0] * 1000]
        legends = legends + [survey]
        colors = colors + [color_dict[survey]]
        data_points_text = data_points_text + [extra_print_elems[j]]

        row_to_print = extra_print_elems[j] + ' & $' + str(int(median_survey_offsets[j - 3])) + ' \\pm ' + str(int(std_survey_offsets[j - 3])) + '$ & $' + ('$ & $'.join([str(can.round_to_n(true_linear_fits[i][0] * 1000 * 0.100, n_sig_fig_errs + round_n_sig_figs[i])) + ' \\pm ' + str(can.round_to_n(boot_std_linear_fits[i] * 1000 * 0.100, n_sig_fig_errs + round_n_sig_figs[i] - 1))  for i in range(3)]))  + '$ & ' + str(n_cepheids_by_survey[survey]) + ' & ' + str(n_HF_by_survey[survey]) + ' & ' + str(can.round_to_n(n_cepheids_by_survey[survey] / n_HF_by_survey[survey], 3)) + ' \\\\'
        print (row_to_print)
    f, ax = plt.subplots(1,1, figsize = [10,5])
    scats = [ax.scatter(x_vals_to_plot[i], y_vals_to_plot[i], c = colors[i], marker = 'x', s = 50) for i in range(len(x_vals_to_plot))]
    ax.set_xlabel(r'$N_{cal} / N_{HF}$', fontsize = 12)
    xticks = [can.round_to_n(tick, 2) for tick in ax.get_xticks()]
    ax.set_xticklabels(xticks, fontsize = 10)
    ax.set_ylabel(r'$d H_0/d (\cdot)}$', fontsize = 12)
    yticks = ax.get_yticks()
    ax.set_yticklabels(yticks, fontsize = 10)
    ax.legend(scats, legends, ncol = 5)
    #[plt.annotate(r"%s" %data_points_text[i] , (x_vals_to_plot[i], y_vals_to_plot[i])) for i in range(len(data_points_text))]
    plt.show()
    n_bins = 100
    fitter_levels = can.niceReverse([0.0] + (1.0 - np.exp(-(np.arange(1.0, 3.1, 1.0) ** 2.0 )/ 2.0)).tolist())
    truth_vals = [70., 0.3, -1] #H0, OmM, w

import cantrips as can
import CosmicFitterClass as cfc
import matplotlib.pyplot as plt
import binData as bd
import numpy as np
import scipy.special as special
import sys
import corner

if __name__=="__main__":
    """
    $ python makePlotOfCosmicParamContoursWithAndWithoutOffsets.py 0.001 0.025 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 #wierd chains - 2, 3, 7, 11, 13, 17, 18,
    #  missing chains - 15
    """
    cl_args = sys.argv[1:]
    base_mu_offset, test_mu_offset, sequence_ids = [cl_args[0], cl_args[1], cl_args[2:]]
    base_mu_offset, test_mu_offset, sequence_ids = [float(base_mu_offset), float(test_mu_offset), [int(id) for id in sequence_ids]]
    print ('[base_mu_offset, test_mu_offset, sequence_ids] = ' + str([base_mu_offset, test_mu_offset, sequence_ids]))

    contour_file_name = 'H0_OmM_w0_contoursWithAndWithout' + str(int(test_mu_offset * 1000)) + 'mMagMuOffset_' + '_'.join([str(id) for id in sequence_ids]) + '.pdf'
    full_fig_file_name = 'FullCornerPlot_' + str(int(test_mu_offset * 1000)) + '_vs_' + str(int(base_mu_offset * 1000)) + 'mMags.pdf'
    save_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/variableMuFits/plots/PosteriorWidths/'
    target_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/variableMuFits/mcmcResults/'
    show_fig = 1
    n_bins = 50
    fitter_levels = can.niceReverse([0.0] + (1.0 - np.exp(-(np.arange(1.0, 3.1, 1.0) ** 2.0 )/ 2.0)).tolist())
    truth_vals = [70., 0.3, -1] #H0, OmM, w
    #base_mu_offset = 0.001
    #test_mu_offset = 0.06
    H0_lims, OmM_lims, w0_lims = [[60.5, 79.5], [0.001, 1.0], [-14.0, -0.0]]
    #lims_std_bounds = [[-6.0, 6.0], [-5.5, 6.0], [-13.5, 6.0] ]
    lims_std_bounds = [[-6.0, 6.0], [-5.5, 6.0], [-7.5, 5.5] ]
    #tick_std_labels = [[-5, -2.5, 0, 2.5, 5], [-5, -2.5, 0, 2.5, 5], [-12, -8, -4, 0, 4 ]]
    tick_std_labels = [[-5, -2.5, 0, 2.5, 5], [-5, -2.5, 0, 2.5, 5], [ -6, -3, 0, 3 ]]
    param_ticks_set = [[68, 69, 70, 71, 72], [0.1, 0.2, 0.3, 0.4, 0.5], [-1.8, -1.4,  -1.0, -0.6]]
    param_ticks_set = [[], [], []]
    param_tick_labels_set = [[], [], []]
    param_lims_set = [H0_lims, OmM_lims, w0_lims]
    param_numbers_in_chain = [1, 2, 3]
    param_labels = ['$\Delta H_0$ (km s$^{-1}$ Mpc$^{-1}$)',  r'$\Delta \Omega_M$', r'$\Delta w$']
    plot_titles = ['$\Delta H_0 = $',  r'$\Delta \Omega_M = $', r'$\Delta w = $']
    ticklabelsize = 14
    plot_elem_size = 3
    n_params = len(param_numbers_in_chain)
    labelsize = 20
    titlesize = 18

    fancy_plot = 1
    if fancy_plot:
        # use LaTeX fonts in the plot
        plt.rc('text', usetex=True)
        plt.rc('font', family='serif')

    surveys = ['SWIFT', 'ASASSN', 'CFA1', 'CFA2', 'LOWZ', 'KAITM', 'CFA4p2', 'KAIT', 'CFA3S', 'CSP', 'CFA3K', 'CFA4p1', 'PS1MD', 'SDSS', 'DES', 'SNLS', 'HST', 'SNAP', 'CANDELS']
    #surveys = ['CANDELS', 'CFA1',  'CFA2', 'CFA3K', 'CFA3S',  ]
    #surveys = ['CANDELS', 'CFA1',  'CFA2', 'CFA3K', 'CFA3S']
    base_full_mcmc_data_files = ['GaussPriorWidth' + str(int(base_mu_offset * 1000)) + 'mMags_zLim_Full_MCMC_pantheon_plus_REAL' + 'H0_OmM_w0' + '_mu' + '_mu'.join(surveys) + '_NS10000_NC50_Cov' + str(i) + '.txt' for i in sequence_ids]
    test_full_mcmc_data_files = ['GaussPriorWidth' + str(int(test_mu_offset * 1000)) + 'mMags_zLim_Full_MCMC_pantheon_plus_REAL' + 'H0_OmM_w0' + '_mu' + '_mu'.join(surveys) + '_NS10000_NC50_Cov' + str(i) + '.txt' for i in sequence_ids]

    #First, make a big corner plot, normalized to the reference plot
    base_data = can.readInColumnsToList(base_full_mcmc_data_files[0], file_dir = target_dir, n_ignore = 1, delimiter = ', ', convert_to_float = 1)
    base_data = [np.array(col[0:-1]) - col[-1] for col in base_data]
    base_data_median = [np.median(col) for col in base_data]
    test_data = can.readInColumnsToList(test_full_mcmc_data_files[0], file_dir = target_dir, n_ignore = 1, delimiter = ', ', convert_to_float = 1)
    test_data = [np.array(col[0:-1]) - col[-1] for col in test_data]
    test_data[0], test_data[1], test_data[2] = [np.array(test_data[0]) - base_data_median[0], np.array(test_data[1]) - base_data_median[1], np.array(test_data[2]) - base_data_median[2]]
    labels = [r'$\Delta H_0$ (km s$^{-1}$ Mpcs$^{-1}$)', r'$\Delta \Omega_M$', r'$\Delta w$'] + [r'$\Delta \mu$' + ''.join([r'$_{}$'.format(char) for char in survey] ) for survey in surveys]


    make_big_corner_plot = 1

    if make_big_corner_plot:
        full_corner_truth_vals = [None, None, None] + [0 for i in range(len(base_data) - 3)]
        big_corner_labelsize = 23
        big_corner_titlesize = big_corner_labelsize - 2
        print ('Making big corner plot... ')
        fig_corner = corner.corner( np.transpose(np.array(test_data[3:] + test_data[0:3])), bins = n_bins, labels=labels[3:] + ['$\Delta H_0$'] + labels[1:3] , truths = full_corner_truth_vals[3:] + full_corner_truth_vals[0:3], levels = fitter_levels,  label_kwargs = {'fontsize':big_corner_labelsize}, title_kwargs  = {'fontsize':big_corner_titlesize, 'pad':9}, show_titles=True, titles = ['' for label in labels])
        #fig_corner.subplots_adjust(top = 20)
        plt.savefig(save_dir + full_fig_file_name)
        plt.close('all')
        print ('Big corner plot made. ')

    #Now make a plot with both surveys stacked on top of each other
    f, axarr = plt.subplots(n_params, n_params, figsize = (plot_elem_size * n_params, plot_elem_size * n_params))
    for i in range(n_params):
        print ('Working on i = ' + str(i))
        param_lims = param_lims_set[i]
        param_ticks = param_ticks_set[i]
        print ('param_lims = ' + str(param_lims))
        print ('param_ticks = ' + str(param_ticks))
        param_col = param_numbers_in_chain[i] - 1
        ax_col, ax_row = [i, i]
        ax = axarr[ax_col, ax_row]
        #bin_x_edges = np.linspace(*param_lims, n_bins + 1)
        base_hists = [[] for mcmc_file in base_full_mcmc_data_files ]
        test_hists = [[] for mcmc_file in test_full_mcmc_data_files ]
        for j in range(len(test_full_mcmc_data_files)):
            base_mcmc_file = base_full_mcmc_data_files[j]
            test_mcmc_file = test_full_mcmc_data_files[j]
            base_data = [float(elem) for elem in can.readInColumnsToList(base_mcmc_file, file_dir = target_dir, n_ignore = 1, delimiter = ', ')[param_col]]
            base_data_offset = base_data[-1]
            base_data = base_data[0:-1]
            test_data = [float(elem) for elem in can.readInColumnsToList(test_mcmc_file, file_dir = target_dir, n_ignore = 1, delimiter = ', ')[param_col]]
            test_data_offset = test_data[-1]
            relative_offset = base_data_offset - test_data_offset
            test_data = test_data[0:-1]
            test_data = [elem + relative_offset for elem in test_data]
            base_data_median, base_data_std = [np.median(base_data), np.std(base_data)]
            print ('base_data_median = ' + str(base_data_median))
            all_data_median, all_data_std = [np.median(test_data + base_data), np.std(test_data + base_data)]
            #param_lims_set[i] = [all_data_median - 5.0 * all_data_std,  all_data_median + 5.0 * all_data_std]
            param_lims_set[i] = [base_data_median + lims_std_bounds[i][0] * base_data_std,  base_data_median + lims_std_bounds[i][1] * base_data_std]
            param_ticks_set[i] = [n_std * base_data_std + base_data_median for n_std in tick_std_labels[i]]
            param_tick_labels_set[i] = [can.round_to_n(tick - base_data_median, 2) for tick in param_ticks_set[i]]
            #param_ticks_set[i] = [can.round_to_n(tick, 2) for tick in param_ticks_set[i]]
            print ('param_lims_set[i] = ' + str(param_lims_set[i]))
            print ('param_ticks_set[i] = ' + str(param_ticks_set[i]))
            #param_lims_set[i] = [param_lims_set[i][0] + base_data_offset, param_lims_set[i][1] + base_data_offset]
            bin_x_edges = np.linspace(*param_lims_set[i], n_bins + 1)
            #param_ticks_set[i] = [can.round_to_n(tick + base_data_offset, 2) for tick in param_ticks_set[i]]
            truth_vals[i] = truth_vals[i] + base_data_offset
            #bin_x_edges = [edge + base_data_offset for edge in bin_x_edges]
            base_binned_data = bd.binData(base_data, [1 for elem in base_data], bin_borders = bin_x_edges, trim = 0)[1]
            #test_binned_data = bd.binData([elem for elem in test_data], [1 for elem in base_data], bin_borders = bin_x_edges, trim = 0)[1]
            test_binned_data = bd.binData(test_data, [1 for elem in base_data], bin_borders = bin_x_edges, trim = 0)[1]
            #print ('base_binned_data = ' + str(base_binned_data))
            #print ('test_binned_data = ' + str(test_binned_data))
            base_hists[j] = np.array([0.0 if np.isnan(elem) else elem for elem in base_binned_data[0]]) * np.array(base_binned_data[2])
            test_hists[j] = np.array([0.0 if np.isnan(elem) else elem for elem in test_binned_data[0]]) * np.array(test_binned_data[2])

            print ('[test_data_offset, base_data_offset, relative_offset] = ' + str([test_data_offset, base_data_offset, relative_offset]))
            base_data_quantiles = np.quantile(base_data, [((special.erf(-1 / np.sqrt(2.0)) + 1 ) / 2 ), 0.5, ((special.erf(1 / np.sqrt(2.0)) + 1 ) / 2 )]) - base_data_median
            base_data_quantiles = [can.round_to_n(quant, 2) for quant in base_data_quantiles]
            test_data_quantiles = np.quantile(test_data, [((special.erf(-1 / np.sqrt(2.0)) + 1 ) / 2 ), 0.5, ((special.erf(1 / np.sqrt(2.0)) + 1 ) / 2 )]) - base_data_median
            test_data_quantiles = [can.round_to_n(quant, 2) for quant in test_data_quantiles]
            #plot_titles[i] = plot_titles[i] + r'$E^{\alpha}_{\beta}$' + r' (no $\Delta \mu_S$)' + '\n' + plot_titles[i] + r'$E^{\alpha}_{\beta}$' + ' (' + str(int(test_mu_offset * 1000)) + r'mMag $\Delta \mu_S$ prior)'
            print ('base_data_quantiles = ' + str(base_data_quantiles))
            print ('test_data_quantiles = ' + str(test_data_quantiles))
            plot_titles[i] = (plot_titles[i]
                              + "${" + str(base_data_quantiles[1]) + "}^{" + str(can.round_to_n(base_data_quantiles[1] - base_data_quantiles[0],2)) + "}_{" + str(can.round_to_n(base_data_quantiles[2] - base_data_quantiles[1],2)) + "}$" + ' (Base)'
                              + '\n' + plot_titles[i]
                              + "${" + str(test_data_quantiles[1]) + "}^{" + str(can.round_to_n(test_data_quantiles[1] - test_data_quantiles[0],2)) + "}_{" + str(can.round_to_n(test_data_quantiles[2] - test_data_quantiles[1],2)) + "}$" + r' ($\Delta \mu_S$)' )
            #plot_titles[i] = plot_titles[i] + r'$E^{\alpha}_{\beta}$'

        param_lims = param_lims_set[i]
        param_ticks = param_ticks_set[i]
        param_tick_labels = param_tick_labels_set[i]
        print ('param_lims = ' + str(param_lims))
        print ('param_ticks = ' + str(param_ticks))
        base_med_hist = np.median(base_hists, axis = 0)
        test_med_hist = np.median(test_hists, axis = 0)
        ax.fill_between(bin_x_edges, np.concatenate(([0],test_med_hist)), step="pre", color = 'white', edgecolor = 'k')
        ax.fill_between(bin_x_edges, np.concatenate(([0],base_med_hist)), step="pre", color = 'r', edgecolor = 'k', alpha = 0.25)
        #ax.axvline(truth_vals[i], c = 'blue', alpha = 0.5, linestyle = '--')
        #ax.axvline()
        ax.set_yticklabels([])
        ax.set_yticks([])
        ax.set_xlim(param_lims)
        ax.set_xticks(param_ticks)
        if ax_row == n_params-1:
            ax.set_xlabel(param_labels[i], fontsize = labelsize)
            #x_param_tick_labels = [can.round_to_n(tick - base_data_median, 2) for tick in param_ticks]
            x_param_tick_labels = param_tick_labels
            print ('Here 3')
            print ('param_ticks = ' + str(param_ticks))
            print ('x_param_tick_labels = ' + str(x_param_tick_labels))
            ax.set_xticklabels(x_param_tick_labels, fontsize = ticklabelsize)
        else:
            ax.set_xticklabels([])
        ax.set_title(plot_titles[i], fontsize = titlesize)

        local_ylims = ax.get_ylim()
        ax.set_ylim([0.0, local_ylims[1]])
    for i in range(n_params - 1):
        x_param_lims = param_lims_set[i]
        x_param_ticks = param_ticks_set[i]
        x_param_tick_labels = param_tick_labels_set[i]
        x_param_col = param_numbers_in_chain[i] - 1
        ax_col = i
        for j in range(i+1, n_params):
            ax_row = j
            print ('[i, j] = ' + str([i, j]) + ' => [ax_col, ax_row] = ' + str([ax_col, ax_row] ))
            y_param_lims = param_lims_set[j]
            y_param_ticks = param_ticks_set[j]
            y_param_tick_labels = param_tick_labels_set[j]
            y_param_col = param_numbers_in_chain[j] - 1
            ax = axarr[ax_row, ax_col]
            base_mcmc_mesh_stack = [[] for k in range(len(base_full_mcmc_data_files))]
            test_mcmc_mesh_stack = [[] for k in range(len(test_full_mcmc_data_files))]
            for k in range(len(test_full_mcmc_data_files)):
                base_mcmc_file = base_full_mcmc_data_files[k]
                test_mcmc_file = test_full_mcmc_data_files[k]
                base_x_data = [float(elem) for elem in can.readInColumnsToList(base_mcmc_file, file_dir = target_dir, n_ignore = 1, delimiter = ', ')[x_param_col]]
                base_x_offset = base_x_data[-1]
                base_x_data = base_x_data[0:-1]
                base_y_data = [float(elem) for elem in can.readInColumnsToList(base_mcmc_file, file_dir = target_dir, n_ignore = 1, delimiter = ', ')[y_param_col]]
                base_y_offset = base_y_data[-1]
                base_y_data = base_y_data[0:-1]
                test_x_data = [float(elem) for elem in can.readInColumnsToList(test_mcmc_file, file_dir = target_dir, n_ignore = 1, delimiter = ', ')[x_param_col]]
                test_x_offset = test_x_data[-1]
                x_relative_offset = base_x_offset - test_x_offset
                test_x_data = [elem + x_relative_offset for elem in test_x_data[0:-1]]
                test_y_data = [float(elem) for elem in can.readInColumnsToList(test_mcmc_file, file_dir = target_dir, n_ignore = 1, delimiter = ', ')[y_param_col]]
                test_y_offset = test_y_data[-1]
                y_relative_offset = base_y_offset - test_y_offset
                test_y_data = [elem + y_relative_offset for elem in test_y_data[0:-1]]
                base_mcmc_samples_to_plot = [base_x_data, base_y_data]
                test_mcmc_samples_to_plot = [test_x_data, test_y_data]
                base_x_mesh, base_y_mesh, base_mcmc_mesh = cfc.getContourMeshFromMCMCChainComponents(base_mcmc_samples_to_plot, x_param_lims, y_param_lims, n_bins )
                test_x_mesh, test_y_mesh, test_mcmc_mesh = cfc.getContourMeshFromMCMCChainComponents(test_mcmc_samples_to_plot, x_param_lims, y_param_lims, n_bins )
                base_mcmc_mesh_stack[k] = base_mcmc_mesh
                test_mcmc_mesh_stack[k] = test_mcmc_mesh
            med_base_mcmc_mesh = np.mean(base_mcmc_mesh_stack, axis = 0)
            med_test_mcmc_mesh = np.mean(test_mcmc_mesh_stack, axis = 0)
            sorted_base_mcmc_samples = np.flip(np.sort(med_base_mcmc_mesh.flatten()))
            sorted_test_mcmc_samples = np.flip(np.sort(med_test_mcmc_mesh.flatten()))
            base_levels = cfc.determineContourLevelsFromMCMCPosterior(sorted_base_mcmc_samples, fitter_levels, samples_already_sorted = 1)
            test_levels = cfc.determineContourLevelsFromMCMCPosterior(sorted_test_mcmc_samples, fitter_levels, samples_already_sorted = 1)
            ax.contourf(base_x_mesh, base_y_mesh, med_base_mcmc_mesh, levels = base_levels, cmap = 'plasma', alpha = 0.5)
            ax.contour(test_x_mesh, test_y_mesh, med_test_mcmc_mesh, levels = test_levels, colors = 'k')
            #ax.axvline(truth_vals[i], c = 'blue', alpha = 0.5, linestyle = '--')
            #ax.axhline(truth_vals[j], c = 'blue', alpha = 0.5, linestyle = '--')
            ax.set_xlim(x_param_lims)
            ax.set_ylim(y_param_lims)
            ax.set_xticks(x_param_ticks)
            ax.set_yticks(y_param_ticks)
            if ax_col == 0:
                print ('Here 1')
                ax.set_ylabel(param_labels[j], fontsize = labelsize)
                print ('y_param_ticks = ' + str(y_param_ticks))
                #y_param_tick_labels = [can.round_to_n(tick - y_param_ticks[len(y_param_ticks) // 2], 2) for tick in y_param_ticks]
                print ('y_param_tick_labels = ' + str(y_param_tick_labels))
                ax.set_yticklabels(y_param_tick_labels, fontsize = ticklabelsize)
            else:
                ax.set_yticklabels([])
            if ax_row == n_params-1:
                print ('Here 2')
                ax.set_xlabel(param_labels[i], fontsize = labelsize)
                print ('x_param_ticks = ' + str(x_param_ticks))
                #x_param_tick_labels = [can.round_to_n(tick - x_param_ticks[len(x_param_ticks) // 2],2) for tick in x_param_ticks]
                print ('x_param_tick_labels = ' + str(x_param_tick_labels))
                ax.set_xticklabels(x_param_tick_labels, fontsize = ticklabelsize)
            else:
                ax.set_xticklabels([])
    for i in range(1, n_params):
        ax_col = i
        for j in range(0, i):
            ax_row = j
            ax = axarr[ax_row, ax_col]
            ax.axis('off')
            #ax.set_xticklabels([])
            #ax.set_xticks([])
            #ax.set_yticklabels([])
            #ax.set_yticks([])

    plt.subplots_adjust (wspace=0.1, hspace=0.1)
    plt.savefig(save_dir + contour_file_name )
    if show_fig:
        plt.show()

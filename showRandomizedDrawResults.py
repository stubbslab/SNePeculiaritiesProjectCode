import matplotlib.pyplot as plt
import numpy as np
import math
from loadSN import loadSN 
from plotGroupedHistograms import plotGroupedHistograms
from loadRandomizedDrawResults import loadRandomizedDrawResults
from calculateSplineWithSpecificKnots import getKnotXPositions
from DirectoryArchive import DirectoryArchive
from AstronomicalParameterArchive import AstronomicalParameterArchive
from CosmologicalParameterArchive import CosmologicalParameterArchive

#At present, cannot make separate plots per observed field 
def showFittedMuResidualCurves(axarr, data_sets_set, column_to_display, surveys_to_fit_set, functions_to_fit_set, best_fit_params_set_of_sets, best_fit_vals_set, ind_var, 
                               pull_extinctions = 1, results_type = '', xlabels= '', ylabels = '', titles = '', n_knots_set = [], spline_min_sep_set = [],
                               x_scale = 'log', legend = 1, legend_labels_set = [['survey']], data_type = 'real', single_shift = 0):

    astro_archive = AstronomicalParameterArchive()
    s_to_yr = astro_archive.getSecondToYear()
    cosmo_archive = CosmologicalParameterArchive()
    age_of_universe = cosmo_archive.getAgeOfUniverse(units = 'yr')[0] 
    
    data_set_dict = {}
    for data_set in np.unique(data_sets_set):
        print 'data_type = ' + str(data_type) 
        data_set_dict[data_set] = loadSN(data_set, pull_extinctions = pull_extinctions, data_type = data_type)
    print 'len(best_fit_params_set_of_sets) = ' + str(len(best_fit_params_set_of_sets)) 
    for i in range(len(data_sets_set)):
        print 'i = ' + str(i) 
        data_set = data_sets_set[i]
        all_sn = data_set_dict[data_set]
        functions_to_fit = functions_to_fit_set[i]
        best_fit_params_set = best_fit_params_set_of_sets[i]
        best_fit_vals = best_fit_vals_set[i]
        best_fit_params = best_fit_params_set[np.argmin(best_fit_vals)]
        median_fit_params = best_fit_params_set[[int(elem) for elem in best_fit_vals].index(int(sorted(best_fit_vals)[len(best_fit_vals) / 2]))]
        worst_fit_params = best_fit_params_set[np.argmax(best_fit_vals)]
        surveys_to_fit = surveys_to_fit_set[i]
        observation_scatter_plots = []
        sample_spline_plots = []
        if surveys_to_fit[0].lower() == 'all':
            surveys_to_fit = np.unique([sn['survey'] for sn in all_sn]).tolist()
            
        min_x = np.inf
        max_x = - np.inf
        for j in range(len(surveys_to_fit)):
            survey = surveys_to_fit[j]
            sn_for_survey = [sn for sn in all_sn if sn['survey'] == survey]
            xs = [sn[ind_var] for sn in sn_for_survey]
            if ind_var is 't':
                xs = [-1.0 * 10 ** (-6.0) * (t * s_to_yr - age_of_universe) for t in xs]
            #Times should be expressed in lookback time
            #if ind_var in ['t']:
            #    xs = [-1.0 * 10 ** (-6.0) * (t * s_to_yr - age_of_universe)  for t in xs]
            new_min = min(xs)
            new_max = max(xs)
            if new_min < min_x: min_x = new_min
            if new_max > max_x: max_x = new_max
        xs_for_function = np.linspace(min_x, max_x, 100)
        for j in range(len(surveys_to_fit)):
            survey = surveys_to_fit[j]
            function_to_fit = functions_to_fit[j]
            sn_for_survey = [sn for sn in all_sn if sn['survey'] == survey]
            color = sn_for_survey[0]['color']
            xs = [sn[ind_var] for sn in sn_for_survey]
            if ind_var is 't':
                xs = [-1.0 * 10 ** (-6.0) * (t * s_to_yr - age_of_universe) for t in xs]
            ys = [sn['muDiff'] for sn in sn_for_survey] 
            yerrs = [sn['muErr'] for sn in sn_for_survey]
            if x_scale is 'log': 
                axarr[i, column_to_display].set_xscale('log', nonposx='clip')
            axarr[i, column_to_display].set_xlim([min_x - (max_x - min_x) / 20.0, max_x + (max_x - min_x) / 20.0])
            observation_scatter_plots = observation_scatter_plots + [axarr[i, column_to_display].scatter(xs, ys, c = color)] 
            axarr[i, column_to_display].errorbar(xs, ys, yerr = yerrs, ecolor = color, fmt = None)
            
            #spline should look like string that looks like: 'technique_spline'
            if 'spline' in results_type.lower() and len(n_knots_set) >  0:
                if single_shift:
                    shift = best_fit_params[0]
                else: 
                    shift = best_fit_params[j] 
                spline_technique = results_type.split('_')[0]
                n_knots = n_knots_set[i]
                knot_args = best_fit_params[len(best_fit_params) - 2 * n_knots:]
                if len(spline_min_sep_set) > 0: min_sep = spline_min_sep_set[i]
                else: min_sep = 0.0
                knot_xs, knot_ys = getKnotXPositions(knot_args, spline_technique, min_sep)
                knot_ys = [y + shift for y in knot_ys]
                sample_spline_plots = sample_spline_plots + [axarr[i, column_to_display].plot(np.array(sorted(xs_for_function.tolist() + knot_xs)), function_to_fit(np.array(sorted(xs_for_function.tolist() + knot_xs)), best_fit_params), c='r')]
                sample_spline_plots = sample_spline_plots + [axarr[i, column_to_display].plot(np.array(sorted(xs_for_function.tolist() + knot_xs)), function_to_fit(np.array(sorted(xs_for_function.tolist() + knot_xs)), median_fit_params), c='c')]
                sample_spline_plots = sample_spline_plots + [axarr[i, column_to_display].plot(np.array(sorted(xs_for_function.tolist() + knot_xs)), function_to_fit(np.array(sorted(xs_for_function.tolist() + knot_xs)), worst_fit_params), c='b')]
                
                #axarr[i, column_to_display].scatter(knot_xs, knot_ys, c = 'r', marker = 'x')
                
            #Otherwise, just plot from data
            
            else: axarr[i, column_to_display].plot(xs_for_function, function_to_fit(xs_for_function, best_fit_params), c = color)
    
            if type(xlabels) is str: xlabel = xlabels
            else: xlabel = xlabels[i]
            if type(ylabels) is str: ylabel = ylabels
            else: ylabel = ylabels[i]
            if type(titles) is str: title = titles
            else: title = titles[i]
            if i == 0:
                axarr[i, column_to_display].set_title(title)
            if i == len(best_fit_params_set[0]) - 1:
                axarr[i, column_to_display].set_xlabel(xlabel)
            axarr[i, column_to_display].set_ylabel(ylabel)
        if legend:
            if len(legend_labels_set) is 1: legend_labels = legend_labels_set[0]
            else: legend_labels = legend_labels_set[i]
            if legend_labels[0].lower() in ['survey','surveys']:
                legend_labels = surveys_to_fit
            axarr[i, column_to_display].legend(observation_scatter_plots + [plot[0] for plot in sample_spline_plots], legend_labels + ['Best Fit Spline', 'Average Fit Spline', 'Worst Fit Spline'], loc = 'upper right', ncol = 7, prop = {'size':6})
        
    

def showRandomizedDrawResults(data_sets, results_files, surveys_to_fit, functions_to_fit, ind_var = 'z', plot_rand_data = 1, plot_true_data = 0, 
                              show = 0, save = 0, figsize = None, results_type = 'standard_spline', fig_name = 'randomized_draw_results_test.png',
                              hist_x_labels = 'Normalized sum of squares', hist_l_ylabels = 'N draws in bin', hist_titles = 'Randomized draw results histogram',
                              bins = 40, hist_bounds = None, hist_r_ylabels = '', data_type = 'real', single_shift = 0, 
                              data_x_labels = 'z', data_y_labels = 'mu residuals', data_titles = 'Residuals data and best fit curve', n_knots_set = [], spline_min_sep_set = [], legend_labels_set = [['survey']]):
    dir_archive = DirectoryArchive()
    plot_dir = dir_archive.getPlotDirectory()
    
    results_arrays = []
    param_sets = []
    true_values = []
    rand_results = []
    for i in range(len(results_files)):
        results_file = results_files[i]
        loaded_results = loadRandomizedDrawResults(results_file)
        results_arrays = results_arrays + [loaded_results[key] for key in loaded_results.keys()] 
    for i in range(len(results_arrays)):
        results_array = results_arrays[i]
        if 'spline' in results_type:
            param_sets = param_sets + [results_array[1]]
            true_values = true_values + [results_array[2]]
            rand_results = rand_results + [results_array[3]]

    if (plot_rand_data and plot_true_data): 
        f, axarr = plt.subplots(len(rand_results), 3, figsize = figsize, squeeze = False)
        axarr_rand_hist_index = 2
        axarr_true_hist_index = 1
    elif plot_rand_data:
        f, axarr = plt.subplots(len(rand_results), 2, figsize = figsize, squeeze = False)
        axarr_rand_hist_index = 1
        axarr_true_hist_index = -1
    elif plot_true_data:
        f, axarr = plt.subplots(len(rand_results), 2, figsize = figsize, squeeze = False)
        axarr_rand_hist_index = -1
        axarr_true_hist_index = 1
    else:
        f, axarr = plt.subplots(len(rand_results), 1, figsize = figsize, squeeze = False)
        axarr_rand_hist_index = -1
        axarr_true_hist_index = -1
    showFittedMuResidualCurves(axarr, data_sets, 0, surveys_to_fit, functions_to_fit, param_sets, true_values, ind_var, 
                               pull_extinctions = 1, results_type = results_type, xlabels = data_x_labels, ylabels = data_y_labels, titles = data_titles,
                               n_knots_set = n_knots_set, spline_min_sep_set = spline_min_sep_set, legend_labels_set = legend_labels_set, data_type = data_type, single_shift = single_shift )
    if plot_rand_data: 
        plotGroupedHistograms(rand_results, axarr = axarr, axarr_hist_index = axarr_rand_hist_index, 
                              labels = [{'set_to_label':i, 'val': np.median(true_values[i]) , 'text':'true val', 'specs':None} for i in range(len(true_values))], add_arrow = 0, 
                              bins = bins, xlabels = hist_x_labels, left_ylabels = hist_l_ylabels, right_ylabels = hist_r_ylabels, titles = hist_titles, bounds = hist_bounds  )
    if plot_true_data: 
        plotGroupedHistograms(true_values, axarr = axarr, axarr_hist_index = axarr_true_hist_index, 
                              labels = [{'set_to_label':i, 'val': np.median(rand_results[i]), 'text':'true val', 'specs':None} for i in range(len(true_values))], add_arrow = 0, 
                              bins = bins, xlabels = hist_x_labels, left_ylabels = hist_l_ylabels, right_ylabels = hist_r_ylabels, titles = hist_titles, bounds = hist_bounds  )

    print 'Finished loading histograms. '

    if save:
        print 'Saving figure to ' + plot_dir + fig_name 
        plt.savefig(plot_dir + fig_name)
        
    if show:
        plt.show()

    #print 'results_arrays = ' + str(results_arrays) 

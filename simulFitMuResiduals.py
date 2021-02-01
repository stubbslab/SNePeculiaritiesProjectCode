import math
import time 
import numpy as np
from simultaneousFit import getSimultaneousFitParams
from SDSSArchive import SDSSArchive
from DirectoryArchive import DirectoryArchive
from PanStars1Archive import PanStars1Archive
from AstronomicalParameterArchive import AstronomicalParameterArchive
from CosmologicalParameterArchive import CosmologicalParameterArchive 
import matplotlib.pyplot as plt 
from loadSN import loadSN
#from SimultaneousFitter import SimultaneousFitter
from checkIfInCyclicRange import checkIfInCyclicRange
from binData import binData
from randomSimulationFunctions import randomSortNSigma

#Assume that all data follow same functional form (fit_function).
#They have some parameters that can be different
# (ones that are allowed to be different between surveys are specified in diff_param_indeces_by_survey). 
#The initial values for parameters are the same for each data set. 
def simulFitMuResiduals(data_set, fit_function, n_fit_params, data_set_type = 'real', guess_params = [], n_steps_per_param = 10, diff_param_indeces_by_survey=[], diff_param_indeces_by_field = [], param_bounds = (),
                        ind_variable = 'z', save = 0, show = 1, surveys_to_fit = ['PS1MD', 'SDSS','SNLS'], fit_information = {'funct':'none'}, data_to_fit = 'residuals', 
                        show_fit_label = 1, fit_label_start = 'fit_', ind_var_range = [0.0001, 1.1], res_range = [-0.8, 0.8], archive_to_use = 'none' , pull_extinctions = 1, plt_spline_points = 0,
                        n_bins = 20, bin_scheme = 'bin_size', binned_res_range = [-0.4, 0.4], zero_ind_var = 0, 
                        rand_arrange = 0, space_sampling_method = 'random', show_binned = 0 ):

    start = time.time() 
    dir_archive = DirectoryArchive()
    astro_archive = AstronomicalParameterArchive()
    cosmo_archive = CosmologicalParameterArchive() 
    plot_dir = dir_archive.getPlotDirectory() 

    PSArch = PanStars1Archive()
    SDSSArch = SDSSArchive()

    if archive_to_use.lower() == 'none':
        archive = None
        fields = {0:astro_archive.getFullSky()}
    elif archive_to_use.lower() == 'sdss':
        archive = SDSSArch
        fields = archive.fields 
    else:
        archive = PSArch
        fields = archive.fields

    s_to_yr = astro_archive.getSecondToYear()
    ageOfUniv = cosmo_archive.getAgeOfUniverse(units = 'yr')[0] 
    
    

    all_sns = loadSN(data_set, surveys_of_interest = ['all'], pull_extinctions = pull_extinctions, data_type = data_set_type)
    
    if surveys_to_fit[0].lower() == 'all':
        surveys_to_fit = np.unique([sn['survey'] for sn in all_sns]).tolist()
    
    sn_by_survey = [[sn for sn in all_sns if sn['survey'] == survey] for survey in surveys_to_fit ]

    
    colors = [sn_set[0]['color'] for sn_set in sn_by_survey] 

    sn_in_survey_by_field = []
    flat_survey_groupings = []
    for sn_set in sn_by_survey:
        sn_by_field = {}
        for key in fields:
            field = fields[key]
            #sn_by_field[key] = [sn for sn in sn_set if (sn['RA']>field[0] and sn['RA'] < field[1] and sn['Dec'] > field[2] and sn['Dec'] < field[3] )]
            sn_by_field[key] = [sn for sn in sn_set if ( checkIfInCyclicRange(sn['RA'], [field[0], field[1]], 360.0) and  checkIfInCyclicRange(sn['Dec'], [field[2], field[3]], 180.0, cyclic_start = -90.0) ) ]
            flat_survey_groupings = flat_survey_groupings + [sn_by_field[key]]
        sn_in_survey_by_field = sn_in_survey_by_field + [sn_by_field]
        flat_survey_groupings 
 
    if ind_variable.lower()  in ['z', 'zs', 'redshift']:
        x_arrays = [[sn['z'] for sn in sn_set] for sn_set in flat_survey_groupings ]
    elif ind_variable.lower() in ['ext', 'exts', 'extinction']:
        x_arrays = [[sn['extinction'] for sn in sn_set] for sn_set in flat_survey_groupings ]
    elif ind_variable.lower() in ['t', 'ts', 'time']:
        time_arrays_in_seconds = [[sn['t'] for sn in sn_set] for sn_set in flat_survey_groupings ]
        #Want to plot in terms of lookback time in Myrs
        x_arrays = [[-1.0 * 10 ** (-6.0) *  (t * s_to_yr - ageOfUniv) for t in time_array] for time_array in time_arrays_in_seconds]
    if  data_to_fit == 'mu' or data_to_fit == 'mus' or data_to_fit == 'distance modulus':
        y_arrays = [[sn['mu'] for sn in sn_set ] for sn_set in flat_survey_groupings]
    else:
        y_arrays = [[sn['muDiff'] for sn in sn_set] for sn_set in flat_survey_groupings ]
        
    y_err_arrays = [[sn['muErr'] for sn in sn_set] for sn_set in flat_survey_groupings ]
    
    if zero_ind_var:
        print 'Subtracting weighted mean, by surveys... '
        mean_sub_y_arrays = []
        for i in range(len(y_arrays)):
            survey = surveys_to_fit[i]
            ys = y_arrays[i]
            y_errs = y_err_arrays[i]
            weights = [1.0 / (err**2.0) for err in y_errs]
            w_mean = sum([ys[i] * weights[i] for i in range(len(weights))]) / sum(weights)
            mean_sub_y_arrays = mean_sub_y_arrays + [[y - w_mean for y in ys]]
        y_arrays = mean_sub_y_arrays
        
    #Allow for doing measurement on randomly sorted data, rather than actual data
    if rand_arrange:
        rand_x_arrays = []
        rand_y_arrays = []
        rand_y_err_arrays = []
        for i in range(len(x_arrays)):
            print 'Randomizing data from survey ' + surveys_to_fit[i] 
            x_array = x_arrays[i]
            y_array = y_arrays[i]
            y_err_array = y_err_arrays[i]
            weighted_mean = sum( [y_array[i] * 1 / ((y_err_array[i]) ** 2.0) for i in range(len(y_array)) ] ) / sum([1 / (err ** 2.0) for err in y_err_array]) 
            rand_x, rand_y, rand_y_err = randomSortNSigma(x_array, y_array, y_err_array, [weighted_mean for y_val in y_array])
            rand_x_arrays = rand_x_arrays + [rand_x]
            rand_y_arrays = rand_y_arrays + [rand_y]
            rand_y_err_arrays = rand_y_err_arrays + [rand_y_err]
        x_arrays = rand_x_arrays
        y_arrays = rand_y_arrays
        y_err_arrays = rand_y_err_arrays
    
    n_versions_of_params = [len(surveys_to_fit) * len(fields) if (i in diff_param_indeces_by_survey and i in diff_param_indeces_by_field)
                               else len(surveys_to_fit)  if (i in diff_param_indeces_by_survey)
                               else len(fields) if (i in diff_param_indeces_by_field)
                               else 1 for i in range(n_fit_params) ]
    n_free_parameters = sum(n_versions_of_params)

    if (type (n_steps_per_param) is int or type(n_steps_per_param) is float) and space_sampling_method.lower() == 'grid':
        n_steps_per_param = [n_steps_per_param for bound in param_bounds[0]]
    
    param_indeces = [[] for i in range(len(flat_survey_groupings))]
    full_guess_params = []
    n_prior_indeces = 0

    if len(param_bounds) == 0:
        param_bounds = ([-np.inf for i in range(n_fit_params) ], [np.inf for i in range(n_fit_params)])
    full_param_bounds = [param_bounds[0] if (type(param_bounds[0]) is float or type(param_bounds[0]) is int or param_bounds[0] is np.inf) else [],
                         param_bounds[1] if (type(param_bounds[1]) is float or type(param_bounds[1]) is int or param_bounds[1] is np.inf) else [] ]
    
    full_guess_parameters = []
    full_n_steps_per_param = []
    for i in range(n_fit_params):
        base_indeces = []
        for j in range(len(surveys_to_fit)):
            for k in range(len(fields)):
                base_indeces = base_indeces + [(j * len(fields) + k if (i in diff_param_indeces_by_survey and i in diff_param_indeces_by_field)
                                                else j if (i in diff_param_indeces_by_survey)
                                                else k if (i in diff_param_indeces_by_field)
                                                else 0)]
        param_indeces = [param_indeces[j] + [n_prior_indeces + base_indeces[j]] for j in range(len(param_indeces)) ]
        n_prior_indeces = n_prior_indeces + len(np.unique(base_indeces))
        if len(guess_params) != 0:
            #print 'guess_params = ' + str(guess_params) 
            full_guess_parameters = full_guess_parameters + (np.zeros(len(fields) * len(surveys_to_fit) if (i in diff_param_indeces_by_survey and i in diff_param_indeces_by_field)
                                                                       else len(surveys_to_fit) if (i in diff_param_indeces_by_survey)
                                                                       else len(fields) if (i in diff_param_indeces_by_field)
                                                                       else 1) + guess_params[i] ).tolist()
            
        full_param_bounds[0] = (full_param_bounds[0] if (type(full_param_bounds[0]) is float or type(full_param_bounds[0]) is int or full_param_bounds[0] is np.inf)
                                else full_param_bounds[0] + (np.zeros(len(fields) * len(surveys_to_fit) if (i in diff_param_indeces_by_survey and i in diff_param_indeces_by_field)
                                                                       else len(surveys_to_fit) if (i in diff_param_indeces_by_survey)
                                                                       else len(fields) if (i in diff_param_indeces_by_field)
                                                                       else 1) + param_bounds[0][i] ).tolist()
                                )
        full_param_bounds[1] = (full_param_bounds[1] if (type(full_param_bounds[1]) is float or type(full_param_bounds[1]) is int or full_param_bounds[1] is np.inf)
                                else full_param_bounds[1] + (np.zeros(len(fields) * len(surveys_to_fit) if (i in diff_param_indeces_by_survey and i in diff_param_indeces_by_field)
                                                                       else len(surveys_to_fit) if (i in diff_param_indeces_by_survey)
                                                                       else len(fields) if (i in diff_param_indeces_by_field)
                                                                       else 1) + param_bounds[1][i] ).tolist()
                                )
        if space_sampling_method.lower() == 'grid':
            full_n_steps_per_param = full_n_steps_per_param + (np.zeros(len(fields) * len(surveys_to_fit) if (i in diff_param_indeces_by_survey and i in diff_param_indeces_by_field)
                                                                           else len(surveys_to_fit) if (i in diff_param_indeces_by_survey)
                                                                           else len(fields) if (i in diff_param_indeces_by_field)
                                                                           else 1) + n_steps_per_param[i] ).tolist()
        else:
            full_n_steps_per_param = n_steps_per_param 
        

    full_param_bounds = (full_param_bounds[0], full_param_bounds[1])
    best_fit, best_fit_errors = getSimultaneousFitParams(x_arrays, y_arrays, [fit_function for sn_set in flat_survey_groupings], param_indeces, guess_params = full_guess_parameters, y_err_arrays = y_err_arrays, param_bounds=full_param_bounds, n_steps = n_steps_per_param, space_sampling_method = space_sampling_method)
    #fitter = SimultaneousFitter(x_arrays, y_data_arrays, [fit_function for sn_set in flat_survey_groupings], param_indeces, [guess_parameters for sn_set in flat_survey_groupings], y_err_arrays = y_err_arrays)
    #print 'best_fit = ' + str(best_fit)
    #print 'best_fit_errors = ' + str(best_fit_errors) 

    chi_sqr = 0.0
    n_sns_fitted = 0

    for field in fields:
        if show:
            if show_binned: 
                f, axarr = plt.subplots(2, sharex = True, figsize = (17.5,9.0), squeeze = 0)
            else:
                f, axarr = plt.subplots(1, sharex = True, figsize = (17.5,9.0), squeeze = 0)
            survey_plots = []
            fitted_plots = []
        #print 'Displaying field ' + str(field)
        fit_position_index = 0
        for i in range(len(surveys_to_fit)):
            sns_in_field_in_survey = sn_in_survey_by_field[i][field]
            zs = [sn['z'] for sn in sns_in_field_in_survey]
            mus = [sn['mu'] for sn in sns_in_field_in_survey]
            muDiffs = [sn['muDiff'] for sn in sns_in_field_in_survey]
            if data_to_fit == 'mus' or data_to_fit == 'distance modulus':
                y_data = mus
            else:
                y_data = muDiffs 
            muErrs = [sn['muErr'] for sn in sn_in_survey_by_field[i][field]]
            x_data = x_arrays[field + i * len(fields)]
            y_data = y_arrays[field + i * len(fields)]
            y_errs = y_err_arrays[field + i * len(fields)]
            #print 'x_data = ' + str(x_data)
            #print 'Doing fit.'
            single_set_param_indeces = param_indeces[field + i * len(fields)]
            single_set_params =  [best_fit[j] for j in single_set_param_indeces ]
            #fitters[field][i].generateFit(zs, y_data, muErrs)
            #print 'For survey ' + surveys_to_fit[i] + ' and field ' + str(field) + ', fit params are: '
            #print single_set_params 
            #print 'For survey ' + surveys_to_fit[i] + ' and field ' + str(field) + ', fit param indeces are: '
            #print single_set_param_indeces 

            #Now compute the contribution of this fit to the chi_square
            #print 'single_set_params = ' + str(single_set_params) 
            #chi_sqr = chi_sqr + sum([ (y_data[j] - fit_function(x_data[j], *single_set_params))**2.0 / y_errs[j]**2.0 for j in range(len(x_data)) ])
            chi_sqr = chi_sqr + sum( (np.array(y_data) - fit_function(x_data, *single_set_params)) ** 2.0 / np.array(y_errs) ** 2.0) 
            n_sns_fitted = n_sns_fitted + len(y_data)
            #print 'n_sns_fitted is now ' + str(n_sns_fitted) 
            #print 'chi_sqr is now = ' + str(chi_sqr)
            #Now plot, if we want to
            if show:
                color = colors[i]
                survey_plots = survey_plots + [axarr[0,0].scatter(x_data, y_data, c = color) ]
                axarr[0,0].errorbar(x_data, y_data, yerr = y_errs, ecolor = color, fmt = None)
            
                if (len(x_data) > 1 and show_binned):
                    z_bin_centers, binned_Sn_data = binData(x_data, y_data, y_errs = y_errs, n_bins = n_bins, bin_scheme = bin_scheme)
        
                    axarr[1,0].scatter(z_bin_centers, binned_Sn_data[0], c = color) 
                    axarr[1,0].errorbar(z_bin_centers, binned_Sn_data[0], yerr = binned_Sn_data[1], fmt = None, ecolor = color)
                ind_var_step = 1.0
                extra_points_for_fit = np.arange(ind_var_range[0],ind_var_range[1], ind_var_step).tolist()
                fitted_plots = fitted_plots + [ axarr[0,0].plot( sorted(x_data + extra_points_for_fit ), fit_function(sorted(x_data + extra_points_for_fit ), *single_set_params), c =color ) ]
                if show_binned: 
                    fitted_plots = fitted_plots + [ axarr[1,0].plot( sorted(x_data + extra_points_for_fit ), fit_function(sorted(x_data + extra_points_for_fit ), *single_set_params), c =color ) ]
                if show_fit_label:
                    fit_string = fit_label_start + '[' + ', '.join([str(np.around(param,5)) for param in single_set_params]) + ']'
                    axarr[0,0].text(0.0, res_range[0]+ 0.06 + 0.14 * (fit_position_index) ,fit_string, color = color)
                    if show_binned: 
                        axarr[1,0].text(0.0, binned_res_range[0] + 0.03 + 0.07 * (fit_position_index) ,fit_string, color = color)
                    fit_position_index = fit_position_index + 1

        if show:
            if plt_spline_points:
                shift = single_set_params[0]
                n_knots = (len(single_set_params) - 1 ) / 2
                knot_xs = single_set_params[1:n_knots+1]
                knot_ys = [shift + knot_y for knot_y in single_set_params[n_knots + 1:]]
                axarr[0,0].scatter(knot_xs, knot_ys, c = 'r', marker = 'x')
                if show_binned: 
                    axarr[1,0].scatter(knot_xs, knot_ys, c = 'r') 
            axarr[0,0].set_xlim(ind_var_range[0], ind_var_range[1])
            axarr[0,0].set_ylim(res_range[0], res_range[1])
            axarr[0,0].set_ylabel('mu residual')
            if show_binned: 
                axarr[1,0].set_xlim(ind_var_range[0], ind_var_range[1])
                axarr[1,0].set_ylim(binned_res_range[0], binned_res_range[1])
                axarr[1,0].set_xlabel('z')
                if data_to_fit == 'mus' or data_to_fit == 'distance modulus':
                    axarr[1,0].set_ylabel('mu residual') 
                else: 
                    axarr[1,0].set_ylabel('Binned mu residual')
            

            plt.suptitle('mu residuals in PS1MD field ' + str(field)) 

            plt.legend(survey_plots, surveys_to_fit, loc = 'upper right', prop = {'size':8})

            if save: plt.savefig(plot_dir + 'SN_residuals_v_' + ind_variable + '_' + 'surveys_' + '_'.join(surveys_to_fit) + '_archive_' + archive_to_use + '_field' + str(field) + '_fit_' + fit_label_start + 'bin_' + bin_scheme + str(n_bins) + '.png') 
        
            plt.show() 

    plt.close('all')
    #print 'best_fit = ' + str(best_fit) 

    #print 'chi_square before reduction is = ' + str(chi_sqr) 
    #print 'total of ' + str(n_sns_fitted) + ' SN considered.'
    #print 'Considered a total of ' + str(len(best_fit)) + ' best fit params.' 
    reduced_chi_square = chi_sqr / float(n_sns_fitted - len(best_fit))

    #Now that we have generated the fitted plots, we want to measure the reduced chi square 

    
    #return reduced_chi_square
    end = time.time()
    print 'Took ' + str(end - start) + 's' 
    #return [chi_sqr, best_fit, best_fit_errors, x_arrays, y_arrays, y_err_arrays]
    return [chi_sqr, best_fit,best_fit_errors]

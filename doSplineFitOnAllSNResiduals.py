#We want to fit a spline of a specified number of points to the SN residual data

from simulFitMuResiduals import simulFitMuResiduals
import numpy as np
from loadSN import loadSN
from calculateSplineWithSpecificKnots import getIterableSplineFunction
from DirectoryArchive import DirectoryArchive 
from AstronomicalParameterArchive import AstronomicalParameterArchive
from CosmologicalParameterArchive import CosmologicalParameterArchive


if __name__ == '__main__':
    #parameters at present are: [data_set, data_set_type, ind_variable, surveys to fit, archive to determine fields, n_burn_in_steps, n_true_value_sampling, number_random_draws,
    #                            spline_type,    number of knots, minimum separation between adjacent knots,
    #                            min knot y value, max knot y value,
    #                            min offset value, max offset value,
    #                            whether shift can vary by survey, whether shift can vary by field
    #                            whether knot xs can vary by survey, whether knot ys can vary by survey,
    #                            whether knot xs can vary by field, whether knot ys can vary by field
    #                            ]
    parameter_sets = [#[1, 'real', 'extinction', ['PS1MD'], 'none', 1000, 100, 'recursive_groups', 4, 0.01, -0.2, 0.2, -0.15, 0.15, 0, 0, 0, 0, 0, 0],
                      #[1, 'real', 'extinction', ['PS1MD'], 'none', 1000, 100, 'recursive_groups', 2, 0.01, -0.2, 0.2, -0.15, 0.15, 0, 0, 0, 0, 0, 0],
                      #[1, 'real' , 'extinction', ['PS1MD'], 'none', 1000, 100, 'recursive_groups', 3, 0.01, -0.2, 0.2, -0.15, 0.15, 0, 0, 0, 0, 0, 0],
                      #[1, 'real' , 'extinction', ['PS1MD'], 'none', 1000, 100, 'standard', 4, 0.01, -0.2, 0.2, -0.15, 0.15, 0, 0, 0, 0, 0, 0],
                      [1, 'real' , 't', ['all'], 'none', 2000, 200, 200,  'standard', 4, 1.0, -0.5, 0.5, -0.15, 0.15, 0, 0, 0, 0, 0, 0],
                      [1, 'real' , 't', ['all'], 'none', 2000, 200, 200,  'standard', 5, 1.0, -0.5, 0.5, -0.15, 0.15, 0, 0, 0, 0, 0, 0],
                      [1, 'real' , 't', ['all'], 'none', 2000, 200, 200,  'standard', 6, 1.0, -0.5, 0.5, -0.15, 0.15, 0, 0, 0, 0, 0, 0],
                      [1, 'real' , 't', ['all'], 'none', 2000, 200, 200,  'standard', 7, 1.0, -0.5, 0.5, -0.15, 0.15, 0, 0, 0, 0, 0, 0],
                      [1, 'real' , 't', ['all'], 'none', 2000, 200, 200,  'standard', 8, 1.0, -0.5, 0.5, -0.15, 0.15, 0, 0, 0, 0, 0, 0],
                      [1, 'real' , 't', ['all'], 'none', 2000, 200, 200,  'standard', 9, 1.0, -0.5, 0.5, -0.15, 0.15, 0, 0, 0, 0, 0, 0],
                      [1, 'real' , 't', ['all'], 'none', 2000, 200, 200,  'standard', 10, 1.0, -0.5, 0.5, -0.15, 0.15, 0, 0, 0, 0, 0, 0],
                      [1, 'real' , 't', ['all'], 'none', 2000, 200, 200,  'standard', 11, 1.0, -0.5, 0.5, -0.15, 0.15, 0, 0, 0, 0, 0, 0],
                      [1, 'real' , 't', ['all'], 'none', 2000, 200, 200,  'standard', 12, 1.0, -0.5, 0.5, -0.15, 0.15, 0, 0, 0, 0, 0, 0] 
                     ]
    astro_archive = AstronomicalParameterArchive()
    s_to_yr = astro_archive.getSecondToYear()
    cosmo_archive = CosmologicalParameterArchive()
    age_of_universe = cosmo_archive.getAgeOfUniverse(units = 'yr')[0]
    saverun = 1
    run_denoter = 'spline_vs_realData_in_t_single_shift_true_and_rand_data_1'
    show_best_fit = 0
    data_set_index = 0
    data_set_type_index = data_set_index + 1
    ind_var_index = data_set_type_index + 1
    surveys_index = ind_var_index + 1
    fields_index = surveys_index + 1
    n_burn_index = fields_index + 1
    n_true_sample_index = n_burn_index + 1
    n_random_draws_index = n_true_sample_index + 1
    spline_index = n_random_draws_index +1
    n_knots_index = spline_index + 1
    min_sep_index = n_knots_index + 1
    min_knot_y_index = min_sep_index + 1
    max_knot_y_index = min_knot_y_index + 1
    min_shift_index = max_knot_y_index + 1
    max_shift_index = min_shift_index + 1
    shift_by_survey_index = max_shift_index + 1
    shift_by_field_index = shift_by_survey_index + 1
    x_by_survey_index = shift_by_field_index + 1
    y_by_survey_index = x_by_survey_index + 1
    x_by_field_index = y_by_survey_index + 1
    y_by_field_index = x_by_field_index + 1

    sig_sqr_results = {}

    set_index = 0
    for parameter_set in parameter_sets:
        print 'Starting to compute series with parameter_set: '
        print parameter_set 
        data_set = parameter_set[data_set_index]
        data_set_type = parameter_set[data_set_type_index]
        ind_var = parameter_set[ind_var_index]
        surveys_to_fit = parameter_set[surveys_index]
        field_archive = parameter_set[fields_index]
        n_burn_steps = parameter_set[n_burn_index]
        n_true_samples = parameter_set[n_true_sample_index]
        n_random_draws = parameter_set[n_random_draws_index] 
        spline_type = parameter_set[spline_index]
        n_knots = parameter_set[n_knots_index]
        min_sep = parameter_set[min_sep_index] 
        min_knot_y = parameter_set[min_knot_y_index]
        max_knot_y = parameter_set[max_knot_y_index] 
        min_shift = parameter_set[min_shift_index]
        max_shift = parameter_set[max_shift_index]
        can_shift_vary_by_survey = parameter_set[shift_by_survey_index]
        can_shift_vary_by_field = parameter_set[shift_by_field_index]
        can_x_vary_by_survey = parameter_set[x_by_survey_index]
        can_y_vary_by_survey = parameter_set[y_by_survey_index] 
        can_x_vary_by_field = parameter_set[x_by_field_index]
        can_y_vary_by_field = parameter_set[y_by_field_index]

        if ind_var.lower in ['ext','ext','extinction','extinctions']:
            pull_extinctions = 1
        else:
            pull_extinctions = 0
        print 'pull_extinctions = ' + str(pull_extinctions) 
        all_sn = loadSN(data_set, pull_extinctions = pull_extinctions, data_type = data_set_type  )
        xs = []
        res = []
        if surveys_to_fit[0] == 'all':
            xs = [sn[ind_var] for sn in all_sn]
            res = [sn['muDiff'] for sn in all_sn]
        else:
            for survey in surveys_to_fit:
                xs = xs + [sn[ind_var] for sn in all_sn if sn['survey'] == survey ]
                res = res + [sn['muDiff'] for sn in all_sn if sn['survey'] == survey]
        if ind_var is 't':
            xs = [-1.0 * 10 ** (-6.0) * (t * s_to_yr - age_of_universe) for t in xs]
        min_x = min(xs)
        max_x = max(xs)
        min_res = min(res)
        max_res = max(res)         

        #I decided that it would be better to allow the limits of the results to set the limits of the y-variables. 
        # If one wishes to set bounds by hand, comment these lines out
        min_knot_y = min_res
        max_knot_y = max_res 
        min_shift = min_res
        max_shift = max_res 

        lower_bounds = [min_shift] + [min_x + min_sep * i for i in range(n_knots)] + [min_knot_y for i in range(n_knots)]
        upper_bounds = [max_shift] + [max_x - min_sep * (n_knots - 1 - i) for i in range(n_knots)] + [max_knot_y for i in range(n_knots)]

        diff_param_indeces_by_survey = [0] if can_shift_vary_by_survey else []
        diff_param_indeces_by_field = [0] if can_shift_vary_by_field else []
        diff_param_indeces_by_survey = diff_param_indeces_by_survey + [i+1 for i in range(n_knots)] if can_x_vary_by_survey else diff_param_indeces_by_survey
        diff_param_indeces_by_field = diff_param_indeces_by_field + [i+1 for i in range(n_knots)] if can_x_vary_by_field else diff_param_indeces_by_field
        diff_param_indeces_by_survey = diff_param_indeces_by_survey + [i + n_knots + 1 for i in range(n_knots)] if can_y_vary_by_survey else diff_param_indeces_by_survey
        diff_param_indeces_by_field = diff_param_indeces_by_field + [i + n_knots + 1 for i in range(n_knots)] if can_y_vary_by_field else diff_param_indeces_by_field

        spline_funct = getIterableSplineFunction(min_x, max_x, min_sep, spline_type)
        true_data_chi_sqr_set = []
        true_data_best_fit_params_set = []
        for true_sample in range(n_true_samples):
            print 'Working on true sample ' + str(true_sample) 
            true_data_chi_sqr,true_data_best_fit_params = simulFitMuResiduals(data_set, spline_funct, 2 * n_knots + 1,
                                                                              data_set_type = data_set_type, n_steps_per_param = n_burn_steps, surveys_to_fit = surveys_to_fit, ind_variable= ind_var, 
                                                                              param_bounds = (lower_bounds, upper_bounds),
                                                                              diff_param_indeces_by_field = diff_param_indeces_by_field,
                                                                              diff_param_indeces_by_survey = diff_param_indeces_by_survey, rand_arrange = 0,
                                                                              show = 0, pull_extinctions = pull_extinctions, zero_ind_var = 1)[0:2]
            true_data_chi_sqr_set = true_data_chi_sqr_set + [true_data_chi_sqr]
            true_data_best_fit_params_set = true_data_best_fit_params_set + [true_data_best_fit_params]
        random_draw_chi_sqr = []
        for random_sample in range(n_random_draws):
            print 'Working on random draw ' + str(random_sample) 
            random_draw_chi_sqr = random_draw_chi_sqr + [simulFitMuResiduals(data_set, spline_funct, 2 * n_knots + 1,
                                                                              data_set_type = data_set_type, n_steps_per_param = n_burn_steps, surveys_to_fit = surveys_to_fit, ind_variable= ind_var, 
                                                                             param_bounds = (lower_bounds, upper_bounds), 
                                                                             diff_param_indeces_by_field = diff_param_indeces_by_field, 
                                                                             diff_param_indeces_by_survey = diff_param_indeces_by_survey, rand_arrange = 1,
                                                                             show = 0, pull_extinctions = pull_extinctions, zero_ind_var = 1)[0]]
        sig_sqr_results[set_index] = [n_knots, true_data_best_fit_params_set, true_data_chi_sqr_set, random_draw_chi_sqr]
        set_index = set_index + 1

    print 'sig_sqr_results = '
    print sig_sqr_results 
        
    if saverun:
        dir_archive = DirectoryArchive() 
        save_dir = dir_archive.getRandomizedDrawResultsDir()
        file_name = 'sequenceOfSplineFitsOf' + ind_var + 'vsMuResidualFromRandomizedDraws_run' + run_denoter + '.npy'
        print 'Saving file to ' + file_name 
        np.save(save_dir + file_name, sig_sqr_results) 

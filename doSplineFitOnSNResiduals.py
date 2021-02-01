#We want to fit a spline of a specified number of points to the SN residual data

from simulFitMuResiduals import simulFitMuResiduals
import numpy as np
from loadSN import loadSN
from calculateSplineWithSpecificKnots import getIterableSplineFunction
from DirectoryArchive import DirectoryArchive 


if __name__ == '__main__':
    #parameters at present are: [data_set, surveys to fit, archive to determine fields, n_burn_in_steps, number_random_draws,
    #                            spline_type,    number of knots, minimum separation between adjacent knots,
    #                            min knot y value, max knot y value,
    #                            min offset value, max offset value,
    #                            whether shift can vary by survey, whether shift can vary by field
    #                            whether knot xs can vary by survey, whether knot ys can vary by survey,
    #                            whether knot xs can vary by field, whether knot ys can vary by field
    #                            ]
    parameter_sets = [[1, ['PS1MD'], 'none', 1000, 3, 'recursive_groups', 5, 0.01, -0.2, 0.2, -0.15, 0.15, 0, 0, 0, 0, 0, 0],
                      [1, ['PS1MD'], 'none', 1000, 4, 'recursive_groups', 6, 0.01, -0.2, 0.2, -0.15, 0.15, 0, 0, 0, 0, 0, 0] ]
    saverun = 1
    run_denoter = '_spline_test'
    show_best_fit = 0
    data_set_index = 0
    surveys_index = data_set_index + 1
    fields_index = surveys_index + 1
    n_burn_index = fields_index + 1
    n_random_draws_index = n_burn_index + 1
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
        surveys_to_fit = parameter_set[surveys_index]
        field_archive = parameter_set[fields_index]
        n_burn_steps = parameter_set[n_burn_index]
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

        all_sn = loadSN(data_set)
        zs = []
        if surveys_to_fit[0] == 'all':
            zs = [sn['z'] for sn in all_sn]
        else:
            for survey in surveys_to_fit:
                zs = zs + [sn['z'] for sn in all_sn if sn['survey'] == survey ]
        min_z = min(zs)
        max_z = max(zs)
        
        lower_bounds = [min_shift] + [min_z + min_sep * i for i in range(n_knots)] + [min_knot_y for i in range(n_knots)]
        upper_bounds = [max_shift] + [max_z + min_sep * (n_knots - 1 - i) for i in range(n_knots)] + [max_knot_y for i in range(n_knots)]

        diff_param_indeces_by_survey = [0] if can_shift_vary_by_field else []
        diff_param_indeces_by_field = [0] if can_shift_vary_by_field else []
        diff_param_indeces_by_survey = diff_param_indeces_by_survey + [i+1 for i in range(n_knots)] if can_x_vary_by_survey else diff_param_indeces_by_survey
        diff_param_indeces_by_field = diff_param_indeces_by_field + [i+1 for i in range(n_knots)] if can_x_vary_by_field else diff_param_indeces_by_field
        diff_param_indeces_by_survey = diff_param_indeces_by_survey + [i + n_knots + 1 for i in range(n_knots)] if can_y_vary_by_survey else diff_param_indeces_by_survey
        diff_param_indeces_by_field = diff_param_indeces_by_field + [i + n_knots + 1 for i in range(n_knots)] if can_y_vary_by_field else diff_param_indeces_by_field

        spline_funct = getIterableSplineFunction(min_z, max_z, min_sep, spline_type)

        true_data_chi_sqr = simulFitMuResiduals(data_set, spline_funct, 2 * n_knots + 1, n_steps_per_param = n_burn_steps, surveys_to_fit = surveys_to_fit, param_bounds = (lower_bounds, upper_bounds),
                                                diff_param_indeces_by_field = diff_param_indeces_by_field, diff_param_indeces_by_survey = diff_param_indeces_by_survey, show = 0)[0:2]
        random_draw_chi_sqr = []
        for i in range(n_random_draws):
            random_draw_chi_sqr = random_draw_chi_sqr + [simulFitMuResiduals(data_set, spline_funct, 2 * n_knots + 1, n_steps_per_param = n_burn_steps, surveys_to_fit = surveys_to_fit,
                                                                             param_bounds = (lower_bounds, upper_bounds),
                                                                             diff_param_indeces_by_field = diff_param_indeces_by_field, 
                                                                             diff_param_indeces_by_survey = diff_param_indeces_by_survey, rand_arrange = 1, show = 0)[0:2]]
        sig_sqr_results[set_index] = [true_data_chi_sqr, random_draw_chi_sqr]
        set_index = set_index + 1
        print 'Finished parameter set ' + str(set_index) 
        
    if saverun:
        dir_archive = DirectoryArchive() 
        save_dir = dir_archive.getRandomizedDrawResultsDir()
        file_name = 'sequenceOfSplineFitsFromRandomizedDraws_run' + run_denoter + '.npy'
        print 'Saving file to ' + file_name 
        np.save(save_dir + file_name, sig_sqr_results) 

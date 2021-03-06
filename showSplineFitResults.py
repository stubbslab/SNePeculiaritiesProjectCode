import numpy as np
from loadRandomizedDrawResults import loadRandomizedDrawResults


def showSplineFitResults(data_set, results_file,
                         results_indeces_to_show = [0], set_of_surveys_to_fit = [['all']], archives_to_use = ['none'],
                         spline_types = ['standard'], minimum_knot_seps = [0.01],
                         ind_vars = ['z'], dep_vars = ['residuals'], 
                         shifts_variabiltity_by_survey = [0], shifts_variabiltity_by_field = [0],
                         knot_xs_variability_by_survey = [0], knot_xs_variability_by_field = [0],
                         knot_ys_variability_by_survey = [0], knot_ys_variability_by_field = [0],
                         plot_fields_separately = 1, plot_surveys_separately = 0, 
     ):
    full_results = loadRandomizedDrawResult(results_file)

    n_knots_index = 0
    best_fit_vals_index = 1
    true_sum_of_sqrs_index = 2
    rand_sum_of_sqrs_index = 3

    for results_index in results_indeces_to_show:
        results_to_show = full_results[results_index]
        surveys_to_fit = set_of_surveys_to_fit[results_index]
        archive_to_use = archives_to_use[results_index]
        spline_type = spline_types[results_index]
        minimum_knot_sep = minimum_knot_seps[results_index]
        ind_var = ind_vars[results_index]
        dep_var = dep_vars[results_index]
        shift_variabiltity_by_survey = shifts_variabiltity_by_survey[results_index]
        shift_variabiltity_by_field = shifts_variabiltity_by_field[results_index]
        knot_x_variabiltity_by_survey = knot_xs_variabiltity_by_survey[results_index]
        knot_x_variabiltity_by_field = knot_xs_variabiltity_by_field[results_index]
        knot_y_variabiltity_by_survey = knot_ys_variabiltity_by_survey[results_index]
        knot_y_variabiltity_by_field = knot_ys_variabiltity_by_field[results_index]
        
        n_knots = results_to_show[n_knots_index]
        best_fit_vals = results_to_show[best_fit_vals_index]
        true_sum_of_sqrs = results_to_show[true_sum_of_sqrs_index]
        rand_sum_of_sqrs = results_to_show[rand_sum_of_sqrs_index]

        

        
        

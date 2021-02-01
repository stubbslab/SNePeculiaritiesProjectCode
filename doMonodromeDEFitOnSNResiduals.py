#We want to fit a spline of a specified number of points to the SN residual data

from simulFitMuResiduals import simulFitMuResiduals
import numpy as np
from loadSN import loadSN
from cantrips import insertListElements
from calculateSplineWithSpecificKnots import getIterableSplineFunction
from computeMuForDifferentCosmologies import computeMuForCosmology 
from CosmologicalParameterArchive import CosmologicalParameterArchive 
from DirectoryArchive import DirectoryArchive
from calculateHForMonodromeDE import getMuForMonodromeDE
from calculateHForMonodromeDE import getParamsForFlattenedHForMonodromeDE
from cantrips import combineTwoDicts 

if __name__ == '__main__':
    #parameters at present are: [set_name,
    #                            data_set, surveys to fit, archive to determine fields, n_burn_in_steps, number_random_draws,
    #                            init_z, flip_z, normalize_to_H0, minimize_average_from_canonical, 
    #                            vary_A,     fixed_A_val,     A_bounds,     vary_A_by_survey,     vary_A_by_field, 
    #                            vary_alpha, fixed_alpha_val, alpha_bounds, vary_alpha_by_survey, vary_alpha_by_field, 
    #                            vary_nu,    fixed_nu_val,    nu_bounds     vary_nu_by_survey,    vary_nu_by_field,
    #                            vary_phi0,  fixed_phi0_val,  phi0_bounds   vary_phi0_by_survey,  vary_phi0_by_field, 
    #                            vary_shift, fixed_shift_val, shift_bounds, vary_shift_by_survey, vary_shift_by_field
    #                            ]
    parameter_sets = [#['oscil_test_1', 
                      # 1, ['PS1MD'], 'none', 100, 1,
                      # 1000.0, 1, 1, 0, 
                      # 1, 0.0, [0.0,0.1], 0, 0,
                      # 0, 0.2, [0.01, 0.3], 0, 0,
                      # 1, 70.0,[50.0, 300.0], 0, 0,
                      # 0, 0.1,[0.0001, 0.2], 0, 0,
                      # 1, 0.0, [-0.1, 0.1], 0, 0 ],
                      #['no_oscil_test_1',
                      # 1, ['PS1MD'], 'none', 50, 1,
                      # 1000.0, 1, 0, 1, 
                      # 0, 0.0, [0.0,0.1], 0, 0,
                      # 1, 0.2, [0.01, 1.0], 0, 0,
                      # 0, 70.0,[50.0, 300.0], 0, 0,
                      # 1, 0.1,[0.0001, 0.2], 0, 0,
                      # 0, 0.0, [-0.1, 0.1], 0, 0 ],
                      ['just_phi0_test_1',
                       1, ['PS1MD'], 'none', 4, 1,
                       1000.0, 1, 0, 0, 
                       1, 0.05, [0.0,0.1], 0, 0,
                       0, 0.2, [0.01, 1.0], 0, 0,
                       1, 170.0,[50.0, 300.0], 0, 0,
                       1, 0.21,[0.0001, 1.0], 0, 0,
                       0, 0.0, [-0.2, 0.2], 0, 0 ]
                       ] 

    saverun = 1
    run_denoter = '_monodrome_test'
    show_best_fit = 0
    set_name_index = 0
    data_set_index = set_name_index + 1
    surveys_index = data_set_index + 1
    fields_index = surveys_index + 1
    n_burn_index = fields_index + 1
    n_random_draws_index = n_burn_index + 1
    init_z_index = n_random_draws_index +1
    flip_z_index = init_z_index + 1
    normalize_H_index = flip_z_index + 1
    minimize_average_index = normalize_H_index + 1

    vary_param_index_modulus = 1
    fixed_param_index_modulus = vary_param_index_modulus + 1
    param_bounds_index_modulus = fixed_param_index_modulus + 1
    param_by_survey_index_modulus = param_bounds_index_modulus + 1
    param_by_field_index_modulus = param_by_survey_index_modulus + 1
    n_param_specifiers = len([vary_param_index_modulus, fixed_param_index_modulus, param_bounds_index_modulus,param_by_survey_index_modulus, param_by_field_index_modulus])

    param_order = ['A','alpha','nu','phi_0','shift']

    sig_sqr_results = {}

    set_index = 0
    for parameter_set in parameter_sets:
        print 'Starting to compute series with parameter_set: '
        print parameter_set
        set_name = parameter_set[set_name_index]
        data_set = parameter_set[data_set_index] 
        surveys_to_fit = parameter_set[surveys_index]
        field_archive = parameter_set[fields_index]
        n_burn_steps = parameter_set[n_burn_index]
        n_random_draws = parameter_set[n_random_draws_index]
        init_z = parameter_set[init_z_index]
        flip_z = parameter_set[flip_z_index ]
        normalize_H = parameter_set[normalize_H_index]
        minimize_average = parameter_set[minimize_average_index]

        vary_params_array = []
        fixed_params_array = []
        params_bounds_array = []
        params_by_survey_array = []
        params_by_field_array = []

        for i in range(len(param_order)):
            base_index = minimize_average_index + i * n_param_specifiers
            vary_params_array = vary_params_array + [parameter_set[base_index + vary_param_index_modulus]]
            fixed_params_array = fixed_params_array + [parameter_set[base_index + fixed_param_index_modulus]]
            params_bounds_array = params_bounds_array + [parameter_set[base_index + param_bounds_index_modulus]]
            params_by_field_array = params_by_field_array + [parameter_set[base_index + param_by_field_index_modulus]]
            params_by_survey_array = params_by_survey_array + [parameter_set[base_index + param_by_survey_index_modulus]]
            
        all_sn = loadSN(data_set)
        zs = []
        if surveys_to_fit[0] == 'all':
            zs = [sn['z'] for sn in all_sn]
        else:
            for survey in surveys_to_fit:
                zs = zs + [sn['z'] for sn in all_sn if sn['survey'] == survey ]
        min_z = min(zs)
        max_z = max(zs)

        lower_bounds_variable_params = []
        upper_bounds_variable_params = []
        lower_bounds_fixed_params = []
        upper_bounds_fixed_params = []
        diff_param_indeces_by_survey = []
        diff_param_indeces_by_field = []
        fixed_params_dict = {}
        
        for i in range(len(param_order)):
            if vary_params_array[i]:
                lower_bounds_variable_params = lower_bounds_variable_params + [params_bounds_array[i][0]]
                upper_bounds_variable_params = upper_bounds_variable_params + [params_bounds_array[i][1]]
            else:
                lower_bounds_fixed_params = lower_bounds_fixed_params + [params_bounds_array[i][0]]
                upper_bounds_fixed_params = upper_bounds_fixed_params + [params_bounds_array[i][1]]
                fixed_params_dict[i] = fixed_params_array[i]
            if params_by_survey_array[i]:
                diff_param_indeces_by_survey = diff_param_indeces_by_survey + [i]
            if params_by_field_array[i]:
                diff_param_indeces_by_field = diff_param_indeces_by_field + [i]
                
        n_free_params = len(lower_bounds_variable_params)

        if minimize_average and normalize_H:
            print 'Cannot both normalize H to H0 at z=0 AND minimize the average separation from expected values.  Will NOT normalize H to H0. '
            normalize_H = 0

        fullMonodromeDEResidualFunct = lambda zs, A, alpha, nu, phi_0, shift: ( getMuForMonodromeDE( zs, A, alpha, nu, phi_0 = phi_0, init_z = init_z, flip_z = flip_z, normalize_H_to_H0 = normalize_H )
                                                                                - computeMuForCosmology(zs)
                                                                                + shift )

        if minimize_average:
            n_var_for_H_calc = 4
            bounds_dict = {}
            fixed_params_to_compute_H_dict = {} #all fixed params that go into determining H(z) (among A, alpha, nu, phi0) 
            extra_params_dict = {} #parameters that dont' go directly into determining H(z) (eg, systematic shift)
            print 'fixed_params_dict = ' + str(fixed_params_dict) 
            for i in range(len(lower_bounds_fixed_params)):
                key = sorted(fixed_params_dict.keys())[i]
                print 'key = ' + str(key) 
                if key < n_var_for_H_calc: 
                    bounds_dict[key] = (lower_bounds_fixed_params[i], upper_bounds_fixed_params[i])
                    fixed_params_to_compute_H_dict[key] = fixed_params_dict[key] 
                else:
                    extra_params_dict[key] = fixed_params_dict[key]

            print 'bounds_dict = ' + str(bounds_dict)
            print 'fixed_params_to_compute_H_dict = ' + str(fixed_params_to_compute_H_dict)
            print 'extra_params_dict = ' + str(extra_params_dict) 
                    
            usedMonodromeDEResidualFunct = lambda zs, *args: fullMonodromeDEResidualFunct(zs, *insertListElements(list(args),
                                                                                                                  combineTwoDicts( getParamsForFlattenedHForMonodromeDE(zs, list(args)[0:n_var_for_H_calc - len(fixed_params_to_compute_H_dict)],
                                                                                                                                                                        fixed_params_to_compute_H_dict,
                                                                                                                                                                        bounds_dict = bounds_dict, init_z = init_z, flip_z = flip_z),
                                                                                                                                   extra_params_dict 
                                                                                                                                 )
                                                                                                                  )
                                                                                          )
            #*args are free parameters
        else: 
            usedMonodromeDEResidualFunct = lambda zs, *args: fullMonodromeDEResidualFunct(zs, *insertListElements(list(args), fixed_params_dict ))
        
        true_data_chi_sqr,true_data_best_fit_params = simulFitMuResiduals(data_set, usedMonodromeDEResidualFunct, n_free_params, n_steps_per_param = n_burn_steps, surveys_to_fit = surveys_to_fit, param_bounds = (lower_bounds_variable_params, upper_bounds_variable_params),
                                                diff_param_indeces_by_field = diff_param_indeces_by_field, diff_param_indeces_by_survey = diff_param_indeces_by_survey, show = 0)[0:2]
        
        #rebuild the fixed_params dict now that I know the best fit values:
        if minimize_average:
            fixed_params_dict = combineTwoDicts( getParamsForFlattenedHForMonodromeDE(zs, true_data_best_fit_params[0:n_var_for_H_calc - len(fixed_params_to_compute_H_dict)],
                                                                                      fixed_params_to_compute_H_dict, bounds_dict = bounds_dict, init_z = init_z, flip_z = flip_z),
                                                 extra_params_dict 
                                               )
        print 'fixed_params_dict = ' + str(fixed_params_dict)
        print 'true_data_best_fit_params = ' + str(true_data_best_fit_params)
        print 'insertListElements(true_data_best_fit_params, fixed_params_dict) = ' + str(insertListElements(true_data_best_fit_params, fixed_params_dict)) 
        true_data_best_fit_params = insertListElements(true_data_best_fit_params, fixed_params_dict)
        print 'true_data_best_fit_params = ' + str(true_data_best_fit_params) 
        random_draw_chi_sqr = []
        for i in range(n_random_draws):
            random_draw_chi_sqr = random_draw_chi_sqr + [simulFitMuResiduals(data_set, usedMonodromeDEResidualFunct, n_free_params, n_steps_per_param = n_burn_steps, surveys_to_fit = surveys_to_fit,
                                                                             param_bounds = (lower_bounds_variable_params, upper_bounds_variable_params),
                                                                             diff_param_indeces_by_field = diff_param_indeces_by_field, 
                                                                             diff_param_indeces_by_survey = diff_param_indeces_by_survey, rand_arrange = 1, show = 0)[0]]
        sig_sqr_results[set_index] = [set_name, true_data_best_fit_params, true_data_chi_sqr, random_draw_chi_sqr]
        set_index = set_index + 1
        print 'Finished parameter set ' + str(set_index) 
        
    if saverun:
        dir_archive = DirectoryArchive() 
        save_dir = dir_archive.getRandomizedDrawResultsDir()
        file_name = 'sequenceOfMonodromeFitsFromRandomizedDraws_run' + run_denoter + '.npy'
        print 'Saving file to ' + file_name 
        np.save(save_dir + file_name, sig_sqr_results) 
    else:
        print 'sig_sqr_results are: '
        print sig_sqr_results 

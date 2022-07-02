
import CosmicFitterClass as cfc
import numpy as np
import matplotlib.pyplot as plt
import cantrips as can
import time
import sys

if __name__=="__main__":
    """
    $ python MakeTableOfCosmicParamContoursFromToyModelAfterMCMCs.py LowZ 0 LowZ 1 MidZ 1 HighZ 1
    """
    cl_args = sys.argv[1:]
    toy_data_file_ids = [cl_args[i * 2] for i in range(len(cl_args) // 2)]
    n_extra_surveys_to_plot = [int(cl_args[i * 2 + 1]) for i in range(len(cl_args) // 2)]
    results_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/variableMuFits/mcmcResults/ToyModel/'
    plot_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/variableMuFits/plots/ToyModel/'
    mu_offsets = [0.001, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.3, 0.5]
    #mu_offsets = [0.001, 0.003, 0.005, 0.007, 0.009, 0.011, 0.013, 0.015, 0.017, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2]
    mu_offsets = [0.001, 0.005, 0.01, 0.03, 0.05, 0.1, 0.15, 0.2]
    mu_offsets = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
    mu_offsets = [0.005, 0.01, 0.025, 0.05]
    #mu_offsets = [ 0.01, 0.025]
    ref_mu_offset = 0.01
    ref_mu_index = mu_offsets.index(ref_mu_offset)
    #mu_offsets = [0.001, 0.01, 0.1]
    z_lims_str = 'zLim_Full'
    #z_lims = [0.0, 0.15]
    #z_lims_str = 'zLim_0p0_0p15'
    cosmicParams_plot_file_name = 'H0OmMw0_Unertainties_VsMuOffsetPrior_' + z_lims_str + '_ToyData.pdf'
    #mu_offsets = [0.0002, 0.2]
    full_steps = 10000
    #full_steps = 2000
    full_n_chains = 50
    OmM_lims, w0_lims = [[0.0, 0.8], [-2.2, -0.3]]
    #full_n_chains = 10
    #surveys = ['CANDELS', 'CFA1',  'CFA2', 'CFA3K', 'CFA3S', 'CFA4p1', 'CFA4p2', 'CSP', 'DES', 'FOUNDATION', 'HST', 'KAIT', 'LOWZ', 'PS1MD', 'SDSS', 'SNAP', 'SNLS', 'SWIFT', 'SWIFTNEW']
    n_iterations = 1

    bins = 50
    fitter_levels = can.niceReverse([0.0] + (1.0 - np.exp(-(np.arange(1.0, 3.1, 1.0) ** 2.0 )/ 2.0)).tolist())

    #all_surveys = ['A','B','C','D','E','F','G','H','I','J', 'H', 'I','J','K','L','M','N','O','P', 'Q']
    #base_surveys = ['CANDELS', 'CFA1', 'CFA2', 'CFA3K', 'CFA3S', 'CFA4p1', 'CFA4p2', 'CSP', 'DES',  'HST', 'KAIT', 'LOWZ', 'PS1MD', 'SDSS', 'SNAP', 'SNLS', 'SWIFT']
    base_surveys = ['SWIFT', 'ASASSN', 'CFA1', 'CFA2', 'LOWZ', 'KAITM', 'CFA4p2', 'KAIT', 'CFA3S', 'CSP', 'CFA3K', 'CFA4p1', 'PS1MD', 'SDSS', 'DES', 'SNLS', 'HST', 'SNAP', 'CANDELS']
    all_surveys_set = [[], ['LowZ'], ['MidZ'], ['HighZ']]
    #all_surveys_set = [[], ['LowZ'], ['MidZ']]
    paper_col_headers = ['Just Pantheon+', '+ Low-$z$ Survey', '+ Mid-$z$ Survey', '+ High-$z$ Survey']

    surveys_files = ['/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/variableMuFits/ArtificialSurveys/ArtificialSNe_muOffset' + str(int(0.01 * 1000)) + '_' + str(toy_data_file_id) + '.csv' for toy_data_file_id in toy_data_file_ids]
    print ('surveys_files = ' + str(surveys_files))
    toy_data_sets = [can.readInColumnsToList(surveys_file, delimiter = ', ', n_ignore = 1) for surveys_file in surveys_files]
    surveys_col = 1
    colors_col = 5
    colors_dicts = [ {toy_data[surveys_col][i]:toy_data[colors_col][i] for i in range(len(toy_data[0]))} for toy_data in toy_data_sets ]
    colors = [[colors_dicts[j][survey] for survey in all_surveys_set[j]] for j in range(len(all_surveys_set))]
    print ('colors = ' + str(colors))
    mu_prior_type = 'gauss'
    all_cols_to_save = []
    #surveys = ['SDSS']
    extra_surveys_to_include_in_sequence = n_extra_surveys_to_plot
    #for survey in surveys_to_plot:
    #    if survey in all_surveys :
    #        extra_surveys_to_include_in_sequence = extra_surveys_to_include_in_sequence + [all_surveys.index(survey) + 1]
    #    else:
    #        extra_surveys_to_include_in_sequence = extra_surveys_to_include_in_sequence  + [0]
    #extra_surveys_to_include_in_sequence = [all_surveys.index(survey)+1 for survey in surveys_to_plot]
    print ('extra_surveys_to_include_in_sequence = ' + str(extra_surveys_to_include_in_sequence))
    f, axarr = plt.subplots(3,1, figsize = (5,6))
    H0_ylims, OmM_ylims, w_ylims = [ [0.21, 0.34], [0.28, 0.62], [0.07, 0.21] ]
    labelsize = 10
    H0_scats = []
    OmM_scats = []
    w_scats = []
    H0_cols_to_print = [['Relative H$_0$ uncertainty when:'] + ['$\sigma_{P,S}  = ' + str(int(offset * 1000)) + '$ mmag' for offset in mu_offsets] ]
    OmM_cols_to_print = [['Relative $\\Omega_M$ uncertainty when:'] + ['$\sigma_{P,S}  = ' + str(int(offset * 1000)) + '$ mmag' for offset in mu_offsets] ]
    w_cols_to_print = [['Relative $w$ uncertainty reduction when:'] + ['$\sigma_{P,S}  = ' + str(int(offset * 1000)) + '$ mmag' for offset in mu_offsets] ]
    OmM_times_w_to_print = [['Relative $\\Omega_M \\times w$ area when:'] + ['$\sigma_{P,S}  = ' + str(int(offset * 1000)) + '$ mmag' for offset in mu_offsets] ]

    print ('extra_surveys_to_include_in_sequence = ' + str(extra_surveys_to_include_in_sequence))
    for j in range(len(extra_surveys_to_include_in_sequence)):
        toy_data_file_id = toy_data_file_ids[j]
        n_extra_surveys = extra_surveys_to_include_in_sequence[j]
        all_surveys = all_surveys_set[j]
        if n_extra_surveys> 0:
            color = colors[j][n_extra_surveys-1]
        else:
            color = 'k'
        surveys_to_include = all_surveys[0:n_extra_surveys]
        print ('surveys_to_include = ' + str(surveys_to_include))
        params_to_vary = ['H0','OmM','w0'] + ['mu' + survey for survey in base_surveys] + ['mu' + survey for survey in surveys_to_include]
        save_files_dict = {mu_offset:['ToySNe_GPW' + str(int(1000 * mu_offset)) + 'mMags_' + z_lims_str + '_' + 'MCMC_' + 'toy_surveys' + '_REAL' + '_'.join(params_to_vary) + '_NS' + str(full_steps) + '_NC' + str(full_n_chains) + '_' + str(toy_data_file_id) + str(iter) + '.txt' for iter in range(1, n_iterations + 1)] for mu_offset in mu_offsets}
        x_vals = list(save_files_dict.keys())
        H0_val_stds_dict = {}
        H0_val_means_dict = {}
        OmM_val_stds_dict = {}
        OmM_val_means_dict = {}
        OmM_times_w_area_dict = {}
        w_val_stds_dict = {}
        w_val_means_dict = {}

        x_vals = list(save_files_dict.keys())
        print ('x_vals = ' + str(x_vals))
        f, axarr = plt.subplots(len(x_vals), 1)
        for k in range(len(x_vals)):
            x_val = x_vals[k]
            mcmc_output_files = save_files_dict [x_val]
            read_in_H0_vals = [[] for i in range(len(mcmc_output_files))]
            read_in_OmM_vals = [[] for i in range(len(mcmc_output_files))]
            read_in_OmM_times_w_vals = [[] for i in range(len(mcmc_output_files))]
            read_in_w_vals = [[] for i in range(len(mcmc_output_files))]
            for i in range(len(mcmc_output_files)):
                mcmc_output_file = mcmc_output_files[i]
                print ('mcmc_output_file = ' + str(mcmc_output_file))
                mcmc_data = can.readInColumnsToList(results_dir + mcmc_output_file, delimiter = ', ', n_ignore = 1, convert_to_float = 1)
                loaded_H0_data = [elem - mcmc_data[0][-1] for elem in mcmc_data[0][:-1]]
                loaded_OmM_data = [elem - mcmc_data[1][-1] for elem in mcmc_data[1][:-1]]
                loaded_w_data = [elem - mcmc_data[2][-1] for elem in mcmc_data[2][:-1]]

                OmM_mesh, w_mesh, mcmc_mesh = cfc.getContourMeshFromMCMCChainComponents([loaded_OmM_data, loaded_w_data], OmM_lims, w0_lims, bins )
                sorted_full_mcmc_mesh = np.flip(np.sort(mcmc_mesh.flatten()))
                full_levels = cfc.determineContourLevelsFromMCMCPosterior(sorted_full_mcmc_mesh, fitter_levels, samples_already_sorted = 1)
                contours = axarr[k].contour(OmM_mesh, w_mesh, mcmc_mesh, levels =  full_levels)

                contour_vertices = [contours.allsegs[-1 - i][0] for i in range(1, len(full_levels)  )] #nth outer contour is -1 - n
                contour_areas = [np.abs(can.measureAreaOfContour(contour)) for contour in contour_vertices]

                H0_std, H0_mean = [np.std(loaded_H0_data), np.mean(loaded_H0_data)]
                OmM_std, OmM_mean = [np.std(loaded_OmM_data), np.mean(loaded_OmM_data)]
                w_std, w_mean = [np.std(loaded_w_data), np.mean(loaded_w_data)]
                read_in_H0_vals[i] = [H0_std, H0_mean]
                read_in_OmM_vals[i] = [OmM_std, OmM_mean]
                read_in_OmM_times_w_vals[i] = [contour_areas[1], 0.0 ]
                read_in_w_vals[i] = [w_std, w_mean]

                print ('read_in_w_vals = ' + str(read_in_w_vals))
                print ('read_in_OmM_times_w_vals = ' + str(read_in_OmM_times_w_vals))
            H0_val_stds_dict[x_val] = [read_in_val[0] for read_in_val in read_in_H0_vals]
            H0_val_means_dict[x_val] = [read_in_val[1] for read_in_val in read_in_H0_vals]
            OmM_val_stds_dict[x_val] = [read_in_val[0] for read_in_val in read_in_OmM_vals]
            OmM_val_means_dict[x_val] = [read_in_val[1] for read_in_val in read_in_OmM_vals]
            w_val_stds_dict[x_val] = [read_in_val[0] for read_in_val in read_in_w_vals]
            w_val_means_dict[x_val] = [read_in_val[1] for read_in_val in read_in_w_vals]
            OmM_times_w_area_dict[x_val] = [read_in_val[0] for read_in_val in read_in_OmM_times_w_vals]
            #print ('H0_val_stds_dict = ' + str(H0_val_stds_dict))
            #print ('H0_val_means_dict = ' + str(H0_val_means_dict))
            print ('OmM_val_stds_dict = ' + str(OmM_val_stds_dict))
            print ('OmM_times_w_area_dict = ' + str(OmM_times_w_area_dict))
            #print ('OmM_val_means_dict = ' + str(OmM_val_means_dict))
            #print ('w_val_stds_dict = ' + str(w_val_stds_dict))
            #print ('w_val_means_dict = ' + str(w_val_means_dict))
        plt.close('all')
        #H0_cols_to_print = H0_cols_to_print + ['+ ' + str(n_extra_surveys) + ' surveys, ' + ', '.join([str(np.mean(H0_val_stds_dict[offset])) for offset in mu_offsets])]
        #OmM_cols_to_print = OmM_cols_to_print + ['+ ' + str(n_extra_surveys) + ' surveys, ' + ', '.join([str(np.mean(OmM_val_stds_dict[offset])) for offset in mu_offsets])]
        #w_cols_to_print = w_cols_to_print + ['+ ' + str(n_extra_surveys) + ' surveys, ' + ', '.join([str(np.mean(w_val_stds_dict[offset])) for offset in mu_offsets])]
        H0_cols_to_print = H0_cols_to_print + [[paper_col_headers[j]] + [np.mean(H0_val_stds_dict[offset]) for offset in mu_offsets] ]
        OmM_cols_to_print = OmM_cols_to_print + [[paper_col_headers[j]] + [np.mean(OmM_val_stds_dict[offset]) for offset in mu_offsets] ]
        OmM_times_w_to_print = OmM_times_w_to_print + [[paper_col_headers[j]] + [np.mean(OmM_times_w_area_dict[offset] ) for offset in mu_offsets] ]
        w_cols_to_print = w_cols_to_print + [[paper_col_headers[j]]+ [np.mean(w_val_stds_dict[offset]) for offset in mu_offsets] ]

    print ('H0_cols_to_print: ')
    print (H0_cols_to_print)
    print ('OmM_cols_to_print: ')
    print (OmM_cols_to_print)
    print ('OmM_times_w_to_print: ')
    print (OmM_times_w_to_print)
    print ('w_cols_to_print: ')
    print (w_cols_to_print)
    ref_H0s = np.array(H0_cols_to_print[1][1:])
    ref_OmMs = np.array(OmM_cols_to_print[1][1:])
    ref_ws = np.array(w_cols_to_print[1][1:])
    ref_H0 = H0_cols_to_print[1][1:][ref_mu_index]
    ref_OmM = OmM_cols_to_print[1][1:][ref_mu_index]
    ref_OmM_w_area = OmM_times_w_to_print[1][1:][ref_mu_index]
    ref_w = w_cols_to_print[1][1:][ref_mu_index]
    for i in range(1, len(H0_cols_to_print)):
        H0_cols_to_print[i][1:] = np.array(H0_cols_to_print[i][1:]) / ref_H0 * 100
        OmM_cols_to_print[i][1:] = np.array(OmM_cols_to_print[i][1:]) / ref_OmM * 100
        OmM_times_w_to_print[i][1:] = np.array(OmM_times_w_to_print[i][1:]) / ref_OmM_w_area * 100
        w_cols_to_print[i][1:] = np.array(w_cols_to_print[i][1:]) / ref_w * 100
    print ('H0_cols_to_print: ')
    print (H0_cols_to_print)
    print ('OmM_cols_to_print: ')
    print (OmM_cols_to_print)
    print ('OmM_times_w_to_print: ')
    print (OmM_times_w_to_print)
    print ('w_cols_to_print: ')
    print (w_cols_to_print)

    H0_rows_to_print = np.transpose(H0_cols_to_print)
    OmM_rows_to_print = np.transpose(OmM_cols_to_print)
    OmM_times_w_rows_to_print = np.transpose(OmM_times_w_to_print)
    w_rows_to_print = np.transpose(w_cols_to_print)
    print ('\\hline')
    for i in range(len(H0_rows_to_print)):
        row = H0_rows_to_print[i]
        #print ('row = ' + str(row))
        if i == 0:
            row_to_print = ' & '.join(str(elem) for elem in row) + ' \\\\'
            print (row_to_print)
            print ('\\hline')
        else:
            row_to_print = row[0] + ' & ' + ' & '.join(str(can.round_to_n(float(elem), 3)) + '\\%' for elem in row[1:]) + ' \\\\'
            print (row_to_print)

    print ('\\hline')
    #for i in range(len(OmM_rows_to_print)):
    #    row = OmM_rows_to_print[i]
    for i in range(len(OmM_times_w_rows_to_print)):
        row = OmM_times_w_rows_to_print[i]
        #print ('row = ' + str(row))
        if i == 0:
            row_to_print = ' & '.join(str(elem) for elem in row) + ' \\\\'
            print (row_to_print)
            print ('\\hline')
        else:
            row_to_print = row[0] + ' & ' + ' & '.join(str(can.round_to_n(float(elem), 3)) + '\\%' for elem in row[1:]) + ' \\\\'
            print (row_to_print)
    print ('\\hline')
    for i in range(len(w_rows_to_print)):
        row = w_rows_to_print[i]
        #print ('row = ' + str(row))
        if i == 0:
            row_to_print = ' & '.join(str(elem) for elem in row) + ' \\\\'
            print (row_to_print)
            print ('\\hline')
        else:
            row_to_print = row[0] + ' & ' + ' & '.join(str(can.round_to_n(float(elem), 3)) + '\\%' for elem in row[1:]) + ' \\\\'
            print (row_to_print)

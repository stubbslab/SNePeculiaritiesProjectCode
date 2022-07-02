import cantrips as can
import numpy as np
import matplotlib.pyplot as plt
import makePlotOfPS1MDFieldsClass as mpfc
import sys

def ConsolidateSNeFieldFits(grid_density, fitter_ids = [1], data_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/SNIsotropyProject/SNeFieldFits/',
                            min_sn = 14, comoving_bin = 100, rand = 0, hemispheres = [0, 1], angle_divs = 30):


        sn_data_type = 'pantheon_plus' #'pantheon_plus' #real - for old Pantheon data
        zHD = 1 # Was 1 for old Pantheon data
        overdensity_param = 200 #Delta, the fraction of background mass density that the average halo mass density must be
        n_good_fits_to_show = 5
        resid_profile_funct = 'exp_void'
        gal_dens_weighting = 0.0
        z_range = [0.0, 0.8]
        #angle_slices = range(1, angle_divs + 1)
        delimiter = ', '
        n_ignore = 4
        n_mcmc_steps = 1500

        init_comoving_bin_guess = 200
        angular_scale_PS1MD = 3.3   #The angular scale of PS1 is 3.3 degree in diameter (see Abstract of: https://iopscience.iop.org/article/10.1088/0004-637X/745/1/42/pdf).  So that times root(2)
        angular_scale_to_bin = angular_scale_PS1MD * 2

        files_to_consolidate = [can.flattenListOfLists([['OverdensityFit' + str(fitter_id) + '_MinNSN_' + str(min_sn) + '_GridDens_' + str(grid_density) + '_BinSize_' + str(comoving_bin) + '_Hemisphere_' + str(hemispheres[i]) + '_RArange_' + str(angle_slice) + 'of' + str(angle_divs) + '_Rand' + str(rand) + '_fits.txt' for angle_slice in range(1, angle_divs+1)] for i in range(len(hemispheres))]) for fitter_id in fitter_ids]
        file_roots = ['OverdensityFit' + str(fitter_id) + '_MinNSN_' + str(min_sn) + '_GridDens_' + str(grid_density) + '_BinSize_' + str(comoving_bin) +  '_Rand' + str(rand) + '_fits'  for fitter_id in fitter_ids]
        files_to_save = [root + '.txt' for root in file_roots]


        best_rows = []

        for i in range(len(files_to_consolidate)):
            print ('angle_divs = ' + str(angle_divs))
            field_plotter = mpfc.PanStarsFieldManager(1, full_sdss_gal_data_file = 'SDSS_fullCoverage_SDSSGals_pzAll.csv', preloaded_sdss_gals = None, gal_dens_weighting = gal_dens_weighting, z_range = z_range, sn_data_type = sn_data_type, zHD = zHD, cutoff_radius_Mpc = init_comoving_bin_guess, NFW_overdensity_param = overdensity_param, resid_from_grav = 0, resid_from_vel = 1, resid_profile_funct = resid_profile_funct, randomize_all_sn = rand, n_mcmc_steps = n_mcmc_steps, hemisphere = hemispheres[0], angle_divs = angle_divs, angle_slice = 0, angular_scale_to_bin = angular_scale_to_bin )
            header = (can.readInRowsToList(files_to_consolidate[i][0], file_dir = data_dir, delimiter = delimiter) [0:n_ignore])
            print ('header = ' + str(header))
            header = [delimiter.join(header_elem) for header_elem in header]
            #print ('header = ' + str(header))
            read_in_data = [can.readInRowsToList(file, file_dir = data_dir, delimiter = delimiter, n_ignore = n_ignore) for file in files_to_consolidate[i]]
            data_to_save = can.flattenListOfLists(read_in_data)
            with open(data_dir + files_to_save[i], 'w') as f:
                print ('Writing consolidated data to ' + str(files_to_save[i]))
                for header_line in header:
                    f.write(header_line + '\n')
                for row in data_to_save:
                    f.write(delimiter.join(row) + '\n')
            print ('Identifying best row...')
            best_row = data_to_save[np.argmin([float(row[5]) / float(row[6]) for row in data_to_save])]
            best_rows = best_rows + [best_row]
            spherical_points = [ [float(row[k]) for k in range(0, 3) ] for row in data_to_save]
            #fit_params =  [ [float(row[k]) for k in range(3, 8) ] for row in data_to_save ]
            #fit_params = [ spherical_points[j] + [float(data_to_save[j][k]) for k in range(3, 5) ] for j in range(len(data_to_save)) ]
            fit_params =  [ [float(row[k]) for k in range(3, 5) ] for row in data_to_save ]
            #fit_results = [ [float(row[k]) for k in range(8, 11) ] for row in data_to_save]
            fit_results = [ [float(row[k]) for k in range(5, 8) ] for row in data_to_save]
            fits = [[fit_params[k]] + fit_results[k] for k in range(len(data_to_save))]


            #print ('spherical_points = ' + str(spherical_points))
            #print ('fit_params = ' + str(fit_params))
            #print ('fit_results = ' + str(fit_results))
            print ('Plotting fit results...')
            field_plotter.makePlotOfOneDFits(spherical_points, fits, save_plot_prefix = file_roots[i], fig_size_unit = 2.5, comoving_bin = comoving_bin)

        return best_rows



if __name__ == "__main__":
    comoving_bin = 150
    rand = 0
    if comoving_bin == 100:
        if rand:
            fitter_ids = [1,2,3,4,5,6,7,8,9,10, 11, 12, 13]
            grid_densities = [24 for id in fitter_ids]
            hemispheres_sets = [ ['both'] for id in fitter_ids]
            angle_divs_set = [60 for id in fitter_ids]
        else:
            grid_densities =  [ 19,       20,    21,      22,       23,       24,        25,       26,       27,       28,       29,        30,     35,      40,     45,      50,     55,       60 ]
            hemispheres_sets = [['both'],[0,1], ['both'], ['both'], ['both'], ['both'], ['both'], ['both'], ['both'], ['both'], ['both'],  [0,1], ['both'], [0,1], ['both'], [0,1], ['both'], ['both']]
            angle_divs_set = [ 30,       30,    30,      30,       30,       30,        30,       30,       30,       30,       30,        30,     30,      30,     30,      30,     30,       30 ]
            fitter_ids = [1]
    elif comoving_bin == 150:
        grid_densities =  [ 25,       27,        30,        30,       32,        35,         40,      50,       60 ]
        hemispheres_sets = [['both'], ['both'],  ['both'],  ['both'],  ['both'], ['both'],  ['both'],  ['both'], ['both']]
        angle_divs_set = [180,         180,        180,        60,        180,      180,       60,       60,       60]
        fitter_ids = [1]

    best_rows = []
    for i in range(len(grid_densities)):
        print ('i = ' + str(i))
        grid_density = grid_densities[i]
        hemispheres = hemispheres_sets[i]
        angle_divs = angle_divs_set [i]
        print ('[grid_density, hemispheres, angle_divs] = ' + str([grid_density, hemispheres, angle_divs]))
        best_row = ConsolidateSNeFieldFits(grid_density, data_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/SNIsotropyProject/SNeFieldFits/',
                                min_sn = 14, comoving_bin = comoving_bin, rand = rand, hemispheres = hemispheres, angle_divs = angle_divs, fitter_ids = fitter_ids)[0]
        best_rows = best_rows + [best_row]
        print ('Ratio of real to null r Chi Sqr = ' + str(float(best_row[5]) / float(best_row[6])))

    for i in range(len(grid_densities)):
        print ('For grid density: ' + str(grid_densities[i]))
        print ('Ratio of real to null r Chi Sqr = ' + str(float(best_rows[i][5]) / float(best_rows[i][6])))

    print ('Done.')

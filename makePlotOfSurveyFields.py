import numpy as np
import math
from PanStars1Archive import PanStars1Archive
from loadSN import loadSN
from binData import binData 
import matplotlib.pyplot as plt
from DirectoryArchive import DirectoryArchive
from FitStorer import FitStorer
from SDSSArchive import SDSSArchive
from checkIfInCyclicRange import checkIfInCyclicRange


def makePlotOfSurveyFields(data_set, n_bins = 20, bin_scheme = 'bin_size', save = 1, show = 0, surveys_to_display = ['PS1MD', 'SDSS','SNLS'], fit_information = {'funct':'none'}, show_fit_label = 1, z_range = [-0.1, 1.1], res_range = [-0.8, 0.8], binned_res_range = [-0.4, 0.4], archive_to_use = 'PS1MD', separate_fields_by_color = 0, separate_fields_by_plot = 1, pull_extinctions = 1):

    dir_archive = DirectoryArchive()
    plot_dir = dir_archive.getPlotDirectory() 
    
    PSArch = PanStars1Archive()
    SDSSArch = SDSSArchive()

    if archive_to_use.lower() == 'sdss':
        archive = SDSSArch
    else:
        archive = PSArch
    fields = archive.fields 

    all_sns = loadSN(1, ['all'], pull_extinctions = pull_extinctions)
    sn_by_survey = [[sn for sn in all_sns if sn['survey'] == survey] for survey in surveys_to_display ]
 
    colors = [sn_set[0]['color'] for sn_set in sn_by_survey]

    sn_in_survey_by_field = []

    for sn_set in sn_by_survey:
        sn_by_field = {}
        for key in fields:
            field = fields[key] 
            #sn_by_field[key] = [sn for sn in sn_set if (sn['RA']>field[0] and sn['RA'] < field[1] and sn['Dec'] > field[2] and sn['Dec'] < field[3] )]
            sn_by_field[key] = [sn for sn in sn_set if ( checkIfInCyclicRange(sn['RA'], [field[0], field[1]], 360.0) and  checkIfInCyclicRange(sn['Dec'], [field[2], field[3]], 180.0, cyclic_start = -90.0) ) ] 
        sn_in_survey_by_field = sn_in_survey_by_field + [sn_by_field]
    if separate_fields_by_color and len(surveys_to_display) == 1:
        field_colors = archive.field_colors
        colors = [field_colors[field] for field in fields]

    fitters = [ [FitStorer(fit_information) for survey in surveys_to_display ] for field in fields]

    if separate_fields_by_color:
        f, axarr = plt.subplots(2, sharex = True, figsize = (10,5))

    n_fields = len(fields) 
    if separate_fields_by_plot:
        f, axarr = plt.subplots( math.ceil(np.sqrt(n_fields  + 1)), math.ceil(np.sqrt(n_fields  + 1))) 
    for field_num in range(n_fields):
        field = fields[field_num] 
        if separate_fields_by_plot: 
            ax = axarr[n_fields / field_num, n_fields % field_num]
        else:
            f, ax = plt.subplots(1, sharex = True, figsize = (10,5))
        print 'Displaying field ' + str(field)
        survey_plots = []
        fitted_plots = []
        fit_position_index = 0
        for i in range(len(surveys_to_display)):
            sns_in_field_in_survey = sn_in_survey_by_field[i][field]
            zs = [sn['z'] for sn in sns_in_field_in_survey]
            muDiffs = [sn['muDiff'] for sn in sns_in_field_in_survey]
            muErrs = [sn['muErr'] for sn in sn_in_survey_by_field[i][field]]
            if not separate_fields_by_color: color = colors[i]
            else: color = colors[field]
            survey_plots = survey_plots + [ax.scatter(zs, muDiffs, c = color) ]
            ax.errorbar(zs, muDiffs, yerr = muErrs, ecolor = color, fmt = None)

            if len(sns_in_field_in_survey) > 0:
                z_bin_centers, binned_Sn_data = binData(zs, muDiffs, y_errs = muErrs, n_bins = n_bins, bin_scheme = bin_scheme)
        
                ax.scatter(z_bin_centers, binned_Sn_data[0], c = color) 
                ax.errorbar(z_bin_centers, binned_Sn_data[0], yerr = binned_Sn_data[1], fmt = None, ecolor = color)
            if fitters[field][i].getDimension() == 1 and len(zs) > 2:
                print 'Doing fit.' 
                fitters[field][i].generateFit(zs, muDiffs, muErrs)
                print 'For survey ' + surveys_to_display[i] + ', fit params are: '
                print (fitters[field][i].fit_params)
            #    #print [ [[ sn['z'] for sn in sns if sn['survey'] == surveys_to_display[i] ], fitters[i].getFitValues([ sn['z'] for sn in sns if sn['survey'] == surveys_to_display[i] ])] for i in range(len(surveys_to_display)) ]
                z_step = 0.001
                extra_points_for_fit = np.arange(z_range[0],z_range[1], z_step).tolist()
                fitted_plots = fitted_plots + [ ax.plot( sorted(zs + extra_points_for_fit ), fitters[field][i].getFitValues(sorted(zs + extra_points_for_fit)), c =color ) ]
                if show_fit_label:
                    fit_string = fitters[field][i].fit_string
                    print 'Showing text: ' + fit_string 
                    ax.text(0.0,-0.5 + 0.1 * (fit_position_index) ,fit_string, color = color)
                    fit_position_index = fit_position_index + 1

            ax.set_xlim(z_range[0], z_range[1])
            axarr[1].set_xlim(z_range[0], z_range[1])

            axarr[0].set_ylim(res_range[0], res_range[1])
            axarr[1].set_ylim(binned_res_range[0], binned_res_range[1])

            axarr[1].set_xlabel('z')
            axarr[1].set_ylabel('Binned mu residual')
            axarr[0].set_ylabel('mu residual')

        plt.suptitle('mu residuals in PS1MD field ' + str(field)) 

        print 'len(survey_plots) = ' + str(len(survey_plots)) 
        plt.legend(survey_plots, surveys_to_display, loc = 'upper left', prop = {'size':8})

        if save: plt.savefig(plot_dir + 'SN_residuals_v_z_' + archive_to_use + '_field' + str(field) + '_bin_' + bin_scheme + str(n_bins) + '_fit_' + fit_information['funct'] + '.png') 
        
        if show and (not separate_fields_by_color or not(len(surveys_to_display) == 1)): plt.show()

    if show and (separate_fields_by_color and len(surveys_to_display) == 1): plt.show()

    plt.close('all') 

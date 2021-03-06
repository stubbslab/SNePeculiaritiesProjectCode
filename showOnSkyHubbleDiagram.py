import numpy as np
import math
from calculateAngularSeparation import calculateAngularSeparation
from SNDataArchive import SNDataArchive 
from RawSNDataStorer import RawSNDataStorer
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from logList import logList
from computeMuForCosmology import computeMuForCosmology
import scipy.integrate as integrate
from binData import binData
from query3dDustMap import query 

#possible plot types: 'on_sky_residual_contour', 'sn_positions', 'residuals_v_z'
def showSNPlot(bin_angular_size, data_set, RA_centers, Dec_centers, plot_type = 'sn_positions', n_levels = 20, surveys_to_display = 'all'):

    xlims = [np.min(RA_centers), np.max(RA_centers)]
    ylims = [np.min(Dec_centers), np.max(Dec_centers)]

    sn_archive = SNDataArchive()
    sn_storer = RawSNDataStorer(data_set)
    sn_ras = sn_storer.getDataArray('RA')
    sn_ras = [float(ra) for ra in sn_ras]
    sn_decs = sn_storer.getDataArray('DECL')
    sn_decs = [float(dec) for dec in sn_decs]
    sn_zs = sn_storer.getDataArray('zHD')
    sn_zs = [float(z) for z in sn_zs]
    sn_mBs = sn_storer.getDataArray('mB')
    sn_mBs = [float(mB) for mB in sn_mBs]
    sn_mBErrs = sn_storer.getDataArray('mBERR')
    sn_mBErrs = [float(mBErr) for mBErr in sn_mBErrs]
    sn_mus = sn_storer.getDataArray('MU')
    sn_mus = [float(mu) for mu in sn_mus]
    sn_muErrs = sn_storer.getDataArray('MUERR') 
    sn_muErrs = [float(muErr) for muErr in sn_muErrs]
    sn_survey = sn_storer.getSurveys()

    #print sn_survey

    #sns = [ {'RA':sn_ras[i], 'Dec':sn_decs[i], 'z':sn_zs[i], 'mu':sn_mus[i], 'muErr':sn_muErrs[i], 'suvey':sn_survey[i]} for i in range(len(sn_zs)) ]

    unique_surveys = np.unique(np.array(sn_survey))
    #print 'unique_surveys = '
    #print unique_surveys 
    #surveys_to_display = surveys_to_display.lower()
    if surveys_to_display == 'all':
        surveys_to_display = unique_surveys 
    
    survey_color_map = sn_archive.getSurveyColorMap() 
    sn_colors = np.array([ survey_color_map[survey] for survey in sn_survey ])

    #print sn_colors 
    memberships = []
    expected_mus = [computeMuForCosmology(z) for z in sn_zs]
    #print 'expected_mus = '
    #print expected_mus

    mu_diffs = [sn_mus[i] - expected_mus[i] for i in range(len(sn_mBs))]
    #print 'mu_diffs = '
    #print mu_diffs
    
    sns = [ {'RA':sn_ras[i], 'Dec':sn_decs[i], 'z':sn_zs[i], 'mu':sn_mus[i], 'muDiff':mu_diffs[i], 'muErr':sn_muErrs[i], 'survey':sn_survey[i], 'color':sn_colors[i]} for i in range(len(sn_zs)) ]

    RAmesh,Decmesh = np.meshgrid(RA_centers,Dec_centers)
    #print 'RAmesh ='
    #print RAmesh
    #print 'Decmesh ='
    #print Decmesh

    #print calculateAngularSeparation(RAmesh, Decmesh, sn_ras[0], sn_decs[0]) 
    angular_sep_meshes = [calculateAngularSeparation(RAmesh, Decmesh, sn_ras[i], sn_decs[i]) for i in range(len(sn_ras)) ]
    

    sn_membership_array = [angular_sep_mesh <= bin_angular_size for angular_sep_mesh in angular_sep_meshes]
    n_sn_array = sum(sn_membership_array)
    sn_mean_diff_array = np.nan_to_num ( sum([ sn_membership_array[i] * mu_diffs[i] for i in range(len(mu_diffs)) ]) / n_sn_array )  
    sn_mean_diff_err_array = np.nan_to_num ( np.sqrt(sum([ sn_membership_array[i] * (sn_muErrs[i]) ** 2.0 for i in range(len(mu_diffs)) ])) / n_sn_array )

    #normal_distribution = lambda x, mu, sigma: 1.0 / np.sqrt(2.0 * math.pi * sigma ** 2.0) * np.exp( - (x - mu)**2.0 / (2.0 * sigma ** 2.0) )
    #normal_probability = lambda x, mu, sigma: integrate.quad(lambda y: 1.0 / np.sqrt(2.0 * math.pi * sigma ** 2.0) * np.exp( - (y - mu)**2.0 / (2.0 * sigma ** 2.0) ), 0.0, 100  )
    n_sigma = np.abs(sn_mean_diff_array) / sn_mean_diff_err_array
    n_sigma = np.nan_to_num(n_sigma) 
    #print 'n_sigma = '
    #print n_sigma

    #print 'sn_mean_diff_array = '
    #print sn_mean_diff_array 
    #print 'sn_mean_diff_err_array = '
    #print sn_mean_diff_err_array

    #print 'sn_mean_diff_array = '
    #print sn_mean_diff_array

    #if show_sn or show_contours or show_residuals_vs_z:
    #    fig = plt.figure(figsize = [9.0,9.0])

    if plot_type == 'sn_positions':
        plt.xlim(xlims)
        plt.ylim(ylims)
        #print [survey_color_map[survey] for survey in surveys_to_display]
        sn_points_by_survey = [plt.scatter([sn['RA'] for sn in sns if sn['survey'] == survey],
                                           [sn['Dec'] for sn in sns if sn['survey'] == survey],
                                           c = survey_color_map[survey]) for survey in surveys_to_display]
        plt.legend(sn_points_by_survey, unique_surveys,   loc='lower left',
                   ncol=4)
        plt.title('SN Positions on sky')
        plt.xlabel('RA')
        plt.ylabel('Dec') 
        #CB_contour = plt.colorbar(shrink=0.8,extend='both',format='%.2e')

    if plot_type == 'residuals_v_z':
        f, axarr = plt.subplots(2, sharex = True)
        sn_points_by_survey = [axarr[0].scatter([sn['z'] for sn in sns if sn['survey'] == survey],
                                                [sn['muDiff'] for sn in sns if sn['survey'] == survey],
                                                c = survey_color_map[survey]) for survey in surveys_to_display]
        [axarr[0].errorbar([sn['z'] for sn in sns if sn['survey'] == survey],
                           [sn['muDiff'] for sn in sns if sn['survey'] == survey],
                           yerr = [sn['muErr'] for sn in sns if sn['survey'] == survey],
                           ecolor = survey_color_map[survey], fmt = None) for survey in surveys_to_display]
        #axarr[0].errorbar(sn_zs, mu_diffs, yerr = sn_muErrs, fmt = None)
        n_bins = 20.0 
        z_bin_centers, binned_Sn_data = binData([sn['z'] for sn in sns if sn['survey'] in surveys_to_display],
                                                [sn['muDiff'] for sn in sns if sn['survey'] in surveys_to_display],
                                                y_errs = [sn['muErr'] for sn in sns if sn['survey'] in surveys_to_display], n_bins = 20.0)
        #print z_bin_centers = '
        #print z_bin_centers
        #print 'binned_Sn_data = '
        #print binned_Sn_data
        axarr[1].scatter(z_bin_centers, binned_Sn_data[0]) 
        axarr[1].errorbar(z_bin_centers, binned_Sn_data[0], yerr = binned_Sn_data[1], fmt = None)
        #axarr[1].plot(dist_bins_centers, best_fit_curve,color = 'red') 
        axarr[0].set_xlabel('z')
        axarr[1].set_xlabel('z')
        axarr[0].set_ylabel('Individual Sn mu residual')
        axarr[1].set_ylabel('Binned Sn mu residual')
        axarr[0].set_title('Mu residuals vs redshift')

        axarr[0].legend(sn_points_by_survey, unique_surveys, loc = 'upper left' )

    if plot_type == 'on_sky_residual_contour':
        #print 'hi'
        min_n_sigma = np.min(n_sigma)
        
        #if min_n <= 0.0: min_n = 0.1
        max_n_sigma = np.max(n_sigma)
        #print 'min_n_sigma = ' + str(min_n_sigma)
        #print 'max_n_sigma = ' + str(max_n_sigma)
        #log_levels = logList(min_n,max_n,n_levels)
        lin_levels = np.arange(min_n_sigma, max_n_sigma + (max_n_sigma - min_n_sigma) / 2.0, (max_n_sigma - min_n_sigma) / 20.0)
        #print lin_levels 
        #plt.contour(RAmesh,Decmesh,sn_mean_diff_array, levels = log_levels, norm = LogNorm())
        plt.contour(RAmesh,Decmesh,n_sigma, levels = lin_levels)
        CB_contour = plt.colorbar(shrink=0.8,extend='both',format='%.2e')
        plt.title('n-sigma of expected vs measured mu') 
        plt.xlabel('RA')
        plt.ylabel('Dec')
        
    plt.show() 
        

    return 1.0 

                          

    

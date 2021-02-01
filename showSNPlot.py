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
from loadSN import loadSN
from scipy.interpolate import griddata
from FitStorer import FitStorer
from pylab import rcParams


#possible plot types: 'on_sky_residual_contour', 'sn_positions', 'residuals_v_z', 'residual_v_extinction_and_z', 'residual_v_extinction_and_z_contour'
def showSNPlot(data_set, bin_angular_size = None, RA_centers = None, Dec_centers=None, plot_types = ['sn_positions'], fit_information = {'funct': 'none'}, n_levels = 20, surveys_to_display = ['all'], legend = 1, n_bins = 20, bin_scheme = 'bin_size', z_range = None, resid_range = None, show_fit_label = 1, figsize = [12,8], show = 1, do_binning = 1, return_plot_object = 0):

    if RA_centers is None: RA_centers = np.arange(0.0,360.0, 1.0)
    if Dec_centers is None: Dec_centers = np.arange(-90.0, 90.0, 1.0)

    rcParams['figure.figsize'] = figsize[0],figsize[1]


    memberships = []
    
    sns = loadSN(data_set, surveys_to_display, pull_extinctions = 0)
    if surveys_to_display[0].lower() == 'all':
        surveys_to_display = np.unique([sn['survey'] for sn in sns]).tolist()
    if not(fit_information['funct'] is 'none' or fit_information['funct'] is None): 
        if 'SNAP' in surveys_to_display: surveys_to_display.remove('SNAP')
    fitters = [ FitStorer(fit_information) for survey in surveys_to_display ] 
    
    x_lims = [np.min(RA_centers), np.max(RA_centers)]
    y_lims = [np.min(Dec_centers), np.max(Dec_centers)]

    if 'sn_positions' in plot_types :
        print 'x_lims = '
        print x_lims
        print 'y_lims = '
        print y_lims 
        plt.xlim(x_lims)
        plt.ylim(y_lims)
        print 'surveys_to_display = '
        print surveys_to_display 
        #print [survey_color_map[survey] for survey in surveys_to_display]
        sn_points_by_survey = [plt.scatter([sn['RA'] for sn in sns if sn['survey'] == survey],
                                           [sn['Dec'] for sn in sns if sn['survey'] == survey],
                                           c = [sn['color'] for sn in sns if sn['survey'] == survey][0])
                               for survey in surveys_to_display ]
        if legend:
            plt.legend(sn_points_by_survey, surveys_to_display,
                       loc='lower left', ncol=4, prop = {'size':6})
        plt.title('SN Positions on sky')
        plt.xlabel('RA')
        plt.ylabel('Dec') 
        #CB_contour = plt.colorbar(shrink=0.8,extend='both',format='%.2e')

    if 'residuals_v_z' in plot_types:
        if do_binning:
            f, axarr = plt.subplots(2, sharex = True, squeeze = False)
        else:
            f, axarr = plt.subplots(1, squeeze = False)
        

        sn_points_by_survey = [axarr[0,0].scatter([sn['z'] for sn in sns if sn['survey'] == survey],
                                                [sn['muDiff'] for sn in sns if sn['survey'] == survey],
                                                c = [sn['color'] for sn in sns if sn['survey'] == survey][0])
                               for survey in surveys_to_display]
        
        [axarr[0,0].errorbar([sn['z'] for sn in sns if sn['survey'] == survey],
                           [sn['muDiff'] for sn in sns if sn['survey'] == survey],
                           yerr = [sn['muErr'] for sn in sns if sn['survey'] == survey],
                           ecolor = [sn['color'] for sn in sns if sn['survey'] == survey][0], fmt = None) for survey in surveys_to_display]
        if fitters[0].getDimension() == 1:
            [fitters[i].generateFit([sn['z'] for sn in sns if sn['survey'] == surveys_to_display[i]],
                                    [sn['muDiff'] for sn in sns if sn['survey'] == surveys_to_display[i]],
                                    [sn['muErr'] for sn in sns if sn['survey'] == surveys_to_display[i]],
                                    ) for i in range(len(surveys_to_display)) ]
            

            for i in range(len(surveys_to_display)):
                print 'For survey ' + surveys_to_display[i] + ', fit params are: '
                print (fitters[i].fit_params)
                #print [ [[ sn['z'] for sn in sns if sn['survey'] == surveys_to_display[i] ], fitters[i].getFitValues([ sn['z'] for sn in sns if sn['survey'] == surveys_to_display[i] ])] for i in range(len(surveys_to_display)) ]
            
            if z_range is None:
                extra_points_for_fit = []
            else:
                z_step = 0.001
                extra_points_for_fit = np.arange(z_range[0],z_range[1], z_step).tolist() 
            fitted_plots = [axarr[0,0].plot(sorted([ sn['z'] for sn in sns if sn['survey'] == surveys_to_display[i] ] + extra_points_for_fit ),
                                            fitters[i].getFitValues(sorted([ sn['z'] for sn in sns if sn['survey'] == surveys_to_display[i] ] + extra_points_for_fit)),
                                            c ='red'               )
                            for i in range(len(surveys_to_display))  ]
            print 'len(fitted_plots) = ' + str(len(fitted_plots))

            if show_fit_label:
                for fitter in fitters:
                    print 'Showing text: '
                    fit_string = fitter.fit_string
                    axarr[0,0].text(0.0,-0.5,fit_string, color = 'red')

            axarr[0,0].set_xlabel('z')
            axarr[0,0].set_ylabel('Individual Sn mu residual')
            axarr[0,0].set_title('Mu residuals vs redshift')

        if do_binning: 
            z_bin_centers, binned_Sn_data = binData([sn['z'] for sn in sns if sn['survey'] in surveys_to_display],
                                                    [sn['muDiff'] for sn in sns if sn['survey'] in surveys_to_display],
                                                    y_errs = [sn['muErr'] for sn in sns if sn['survey'] in surveys_to_display], n_bins = n_bins, bin_scheme = bin_scheme)
            axarr[1,0].scatter(z_bin_centers, binned_Sn_data[0]) 
            axarr[1,0].errorbar(z_bin_centers, binned_Sn_data[0], yerr = binned_Sn_data[1], fmt = None)
            axarr[1,0].set_xlabel('z')
            axarr[1,0].set_ylabel('Binned Sn mu residual')
        

        if not(z_range is None):
            axarr[0,0].set_xlim(z_range)
            if do_binning:
                axarr[1,0].set_xlim(z_range)
        if not(resid_range is None):
            axarr[0,0].set_ylim(resid_range)
            if do_binning: 
                axarr[1,0].set_ylim(resid_range)

        if legend:
            axarr[0,0].legend(sn_points_by_survey, surveys_to_display,
                            loc = 'upper left', prop={'size':8} )

        if return_plot_object: return [f, axarr]

    if 'residuals_v_extinction' in plot_types:
        if do_binning:
            f, axarr = plt.subplots(2, sharex = True, squeeze = False)
        else:
            f, axarr = plt.subplots(1, squeeze = False)
        
        sn_points_by_survey = [axarr[0,0].scatter([sn['extinction'] for sn in sns if sn['survey'] == survey],
                                                [sn['muDiff'] for sn in sns if sn['survey'] == survey],
                                                c = [sn['color'] for sn in sns if sn['survey'] == survey][0])
                               for survey in surveys_to_display]
        
        [axarr[0,0].errorbar([sn['extinction'] for sn in sns if sn['survey'] == survey],
                           [sn['muDiff'] for sn in sns if sn['survey'] == survey],
                           yerr = [sn['muErr'] for sn in sns if sn['survey'] == survey],
                           ecolor = [sn['color'] for sn in sns if sn['survey'] == survey][0], fmt = None) for survey in surveys_to_display]

        axarr[0,0].set_xlabel('extinction')
        axarr[0,0].set_ylabel('Individual Sn mu residual')
        axarr[0,0].set_title('Mu residuals vs extinction')
        if do_binning: 
            ext_bin_centers, binned_Sn_data = binData([sn['extinction'] for sn in sns if sn['survey'] in surveys_to_display],
                                                      [sn['muDiff'] for sn in sns if sn['survey'] in surveys_to_display],
                                                      y_errs = [sn['muErr'] for sn in sns if sn['survey'] in surveys_to_display], n_bins = n_bins, bin_scheme = bin_scheme)
            axarr[1,0].scatter(ext_bin_centers, binned_Sn_data[0]) 
            axarr[1,0].errorbar(ext_bin_centers, binned_Sn_data[0], yerr = binned_Sn_data[1], fmt = None)
        
            axarr[1,0].set_xlabel('extinction')
        
            axarr[1,0].set_ylabel('Binned Sn mu residual')
        

        if legend:
            axarr[0,0].legend(sn_points_by_survey, surveys_to_display,
                              loc = 'upper left' )

        if return_plot_object: return [f, axarr]

    if 'residual_v_extinction_and_z_contour' in plot_types:
        zs = [sn['z'] for sn in sns if (sn['survey'] in surveys_to_display) ]
        exts = [sn['extinction'] for sn in sns if (sn['survey'] in surveys_to_display) ]
        mu_diffs = [sn['muDiff'] for sn in sns if (sn['survey'] in surveys_to_display) ]
        mu_errs = [sn['muErr'] for sn in sns if (sn['survey'] in surveys_to_display) ]

        n_samplings_per_axis = 1000
        #z_mesh, ext_mesh = np.meshgrid( getListOfLength(0.0, max(zs), n_samplings_per_axis), getListOfLength(0.0, max(exts), n_samplings_per_axis) )
        z_mesh, ext_mesh = np.meshgrid(np.linspace(0.0, max(zs),n_samplings_per_axis), np.linspace(0.0, max(exts),n_samplings_per_axis) )
        
        val_to_plot = [ mu_diffs[i] / mu_errs[i] for i in range(len(mu_errs)) ] # pick your value to plot
        #val_to_plot = mu_diffs 
        #print 'val_to_plot = '
        #print val_to_plot
        val_to_plot_mesh = griddata([ [zs[i], exts[i]] for i in range(len(zs)) ], val_to_plot, (z_mesh, ext_mesh), method = 'linear' )

        #print val_to_plot_grid_vals
        print np.shape(val_to_plot_mesh)
        
        resid_contours = plt.contourf(z_mesh,ext_mesh,val_to_plot_mesh,50)
        #plt.imshow(val_to_plot_mesh, extent = [0.0,max(zs) * 1.1, 0.0, max(exts) * 1.1 ], origin = 'lower')
        plt.colorbar(resid_contours) 
        plt.xlabel('z')
        plt.ylabel('extinction')
        plt.title('(mu residuals) / (sigma) for SN in extinction, z space')

        if return_plot_object: return [f, axarr]
        

    if 'residual_v_extinction_and_z' in plot_types:
        sn_points_by_survey = [plt.scatter([sn['z'] for sn in sns if sn['survey'] == survey],
                                           [sn['extinction'] for sn in sns if sn['survey'] == survey],
                                           c = [sn['color'] for sn in sns if sn['survey'] == survey][0])
                               for survey in surveys_to_display ]
        if legend:
            plt.legend(sn_points_by_survey, surveys_to_display,
                       loc='upper right', ncol=4)
        x_lims = [0.0, max([sn['z'] for sn in sns]) * 1.1]
        y_lims = [0.0, max([sn['extinction'] for sn in sns]) * 1.1]
        plt.xlim(x_lims)
        plt.ylim(y_lims)
        
        #plt.title('SN Points in z, extinction space')
        plt.xlabel('z')
        plt.ylabel('extinction') 

    if 'on_sky_residual_contour' in plot_types:
        RAmesh,Decmesh = np.meshgrid(RA_centers,Dec_centers)
        angular_sep_meshes = [calculateAngularSeparation(RAmesh, Decmesh, sns[i]['RA'], sns[i]['Dec']) for i in range(len(sns)) ]
    
        sn_membership_array = [angular_sep_mesh <= bin_angular_size for angular_sep_mesh in angular_sep_meshes]
        n_sn_array = sum(sn_membership_array)
        sn_mean_diff_array = np.nan_to_num ( sum([ sn_membership_array[i] * sns[i]['muDiff'] for i in range(len(sns)) ]) / n_sn_array )  
        sn_mean_diff_err_array = np.nan_to_num ( np.sqrt(sum([ sn_membership_array[i] * (sns[i]['muErr']) ** 2.0 for i in range(len(sns)) ])) / n_sn_array )

        n_sigma = np.abs(sn_mean_diff_array) / sn_mean_diff_err_array
        n_sigma = np.nan_to_num(n_sigma) 
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

                          

    

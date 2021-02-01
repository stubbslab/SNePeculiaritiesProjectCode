import numpy as np
import math
from calculateAngularSeparation import calculateAngularSeparation
from SNDataArchive import SNDataArchive 
from RawSNDataStorer import RawSNDataStorer
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
from logList import logList
from computeMuForCosmology import computeMuForCosmology
import scipy.integrate as integrate
from binData import binData
from query3dDustMap import query
from PanStars1Archive import PanStars1Archive
from loadSN import loadSN
from getListOfLength import getListOfLength
from logList import logList
from matplotlib.colors import LogNorm

#possible plot types: 'on_sky_residual_contour', 'sn_positions', 'residuals_v_z'
def showResidualVsExtinctionOnSky(data_set, plot_type = 'sn_positions', n_contours = 20, surveys_to_display = ['PS1MD'],
                             legend = 1, meshsize = [100,100], min_display_n_sigma = -3.0, max_display_n_sigma = 3.0, n_colorbar_ticks = 10.0):

    sns = loadSN(data_set, surveys_to_display)
    
    if surveys_to_display[0].lower() == 'all':
        surveys_to_display = np.unique([sn['survey'] for sn in sns])

    print surveys_to_display

    star_archive = []

    if surveys_to_display == ['PS1MD']:
        star_archive = PanStars1Archive()

    sn_fields = star_archive.fields

    print 'sn_fields = '
    print sn_fields

    sn_by_field = {key: [sn for sn in sns
                         if (sn['RA'] >= sn_fields[key][0] and sn['RA'] <= sn_fields[key][1]
                             and sn['Dec'] >= sn_fields[key][2] and sn['Dec'] <= sn_fields[key][3])   ]
                   for key in sn_fields   }
    n_sigma_mu_diff = [sn['muDiff'] / sn['muErr'] for sn in sns ]
    print 'n_sigma_mu_diff = '
    print n_sigma_mu_diff
    max_n_sigma = max(n_sigma_mu_diff)
    min_n_sigma = min(n_sigma_mu_diff)
    print 'max_n_sigma = ' + str(max_n_sigma)
    print 'min_n_sigma = ' + str(min_n_sigma)
    
    #print [str(sn['RA']) + ',' + str(sn['Dec'])  for sn in sn_by_field[2]]

    n_fields = len(sn_fields)
    n_plots_per_field = 2 #number of distinct plots displayed per field 

    n_cols = 5
    n_rows = n_fields / n_cols
    n_ticks = 5 #number of ticks to display on each axis
    if n_rows * n_cols < n_fields: n_rows = n_rows + 1

    print 'n_cols = ' + str(n_cols)
    print 'n_rows = ' + str(n_rows)

    fig = plt.figure(figsize = (15, 10))
    outer = gridspec.GridSpec(n_rows, n_cols, wspace = 0.15, hspace = 0.1)

    RA_lim_arrays = [ [sn_fields[key][0],sn_fields[key][1]] for key in sn_fields ]
    Dec_lim_arrays = [ [sn_fields[key][2],sn_fields[key][3]] for key in sn_fields ]
    xtick_arrays = [ np.around(getListOfLength(RA_lim[0], RA_lim[1], n_ticks + 2 )[1:n_ticks + 2-1].tolist(),1) for RA_lim in RA_lim_arrays]
    ytick_arrays = [ np.around(getListOfLength(Dec_lim[0], Dec_lim[1], n_ticks + 2 )[1:n_ticks + 2-1].tolist(),1) for Dec_lim in Dec_lim_arrays]
    
    coord_mesh_arrays = [np.meshgrid( getListOfLength(RA_lim_arrays[i][0],RA_lim_arrays[i][1],meshsize[0]), getListOfLength(Dec_lim_arrays[i][0],Dec_lim_arrays[i][1],meshsize[1]) )
                          for i in range(len(sn_fields))]
    extinction_meshes_array = [np.array(query(coord_mesh[0].tolist(),coord_mesh[1].tolist(), coordsys = 'equ', mode = 'sfd')[u'EBV_SFD']) for coord_mesh in coord_mesh_arrays]
    
    min_ext = 0.0
    min_ext = min([np.min(extinction_mesh) for extinction_mesh in extinction_meshes_array])
    max_ext = max([np.max(extinction_mesh) for extinction_mesh in extinction_meshes_array])

    print 'min_ext = ' + str(min_ext)
    print 'max_ext = ' + str(max_ext)

    #ext_levels = getListOfLength(min_ext, max_ext, n_contours)
    ext_levels = logList(min_ext, max_ext, n_contours)
    print 'ext_levels = '
    print ext_levels

    

    for field_number in range(n_fields):
        field = sn_fields[field_number]
        
        inner =  gridspec.GridSpecFromSubplotSpec(n_plots_per_field, 1, subplot_spec = outer[field_number], wspace=0.0, hspace=0.0)
        
        ax_sn = plt.Subplot(fig,inner[1])
        ax_sn.set_xlim(RA_lim_arrays[field_number])
        ax_sn.set_xticks(xtick_arrays[field_number])
        ax_sn.set_ylim(Dec_lim_arrays[field_number])
        ax_sn.set_yticks(ytick_arrays[field_number])
        ax_sn.tick_params(labelsize = 8)
        scatter_for_colorbar = ax_sn.scatter([sn['RA'] for sn in sn_by_field[field_number]], [sn['Dec'] for sn in sn_by_field[field_number]],
                      c = [sn['muDiff'] / sn['muErr'] for sn in sn_by_field[field_number] ], vmin = min_display_n_sigma , vmax = max_display_n_sigma )
        
        ax_ext = plt.Subplot(fig,inner[0])
        ax_ext.set_xlim(RA_lim_arrays[field_number])
        ax_ext.set_xticks( xtick_arrays[field_number] )
        ax_ext.set_ylim(Dec_lim_arrays[field_number])
        ax_ext.set_yticks(ytick_arrays[field_number])
        ax_ext.set_xticks([])
        RAmesh,Decmesh = coord_mesh_arrays[field_number]
        ax_ext.tick_params(labelsize = 8)
        extinction_mesh = extinction_meshes_array[field_number]
        contourf_for_colorbar = ax_ext.contourf(RAmesh,Decmesh, extinction_mesh, levels = ext_levels, norm = LogNorm() )
        fig.add_subplot(ax_sn)
        fig.add_subplot(ax_ext) 

    fig.subplots_adjust(right = 0.9)
    cbar_ext_ax = fig.add_axes([0.91, 0.15, 0.05, 0.7])
    cbar_n_sigma_ax = fig.add_axes([0.015, 0.15, 0.05, 0.7])
    #cbar_ext_ax = fig.add_axes([0.09, 0.15, 0.05, 0.7])
    #print fig 
    cbar_ext = fig.colorbar (contourf_for_colorbar, cax = cbar_ext_ax, ticks =  ext_levels)
    cbar_n_sigma = fig.colorbar (scatter_for_colorbar, cax = cbar_n_sigma_ax,
                                 ticks =  np.around(getListOfLength(min_display_n_sigma, max_display_n_sigma, n_colorbar_ticks),2) )
    
    #cbar.set_label('This is the colorbar.' )
    print 'ext_levels = '
    print ext_levels
    #cbar.set_yticks(np.around(ext_levels,4))
    cbar_ext.ax.set_yticklabels(np.around(ext_levels,3))
    
    #plt.colorbar(contourf_for_colorbar)
    plt.suptitle('Extinction (upper) vs normalized mu residuals for sn in survey(s) '+','.join(surveys_to_display) + ' (lower). \n All frames are patches on sky (x-axis axis is RA, y-axis is Dec, both in degrees and standard equ. coords.) ', fontsize = 'large')
    #plt.text(1,2,'Here is some text') 
    fig.show() 
    
        

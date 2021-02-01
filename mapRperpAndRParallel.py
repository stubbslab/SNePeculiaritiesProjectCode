import numpy as np
import matplotlib.pyplot as plt 
import matplotlib

#Plot r_{\perp} / S_1 and r_{\parallel} / S_1 as a function of \Delta \theta and S_2 / S_1
# (Note S_2 / S_1 >= 1 since S_2 is always selected to be the larger of the two lengths
def calcComovingRperpAndRParallel(ang_seps, comoving_dist_rats, 
                                  show_fig = 1, save_fig = 1, figsize = (10, 10)): 

    save_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/SNTwoPointCorrelationsProject/plots/'
    save_fig_name = 'rperpAndRParallelBaseOnSkyCoords.png' 
    comoving_center_funct = lambda ang_sep, comoving_dist_rat: 0.5 * (1 ** 2.0 + 2.0 * 1 * comoving_dist_rat * np.cos(ang_sep) + comoving_dist_rat ** 2.0) ** 0.5

    comoving_perps = lambda ang_sep, comoving_dist_rat: np.abs(np.sin(ang_sep)) * 1.0 * comoving_dist_rat / comoving_center_funct(ang_sep, comoving_dist_rat) 
    comoving_parallels = lambda ang_sep, comoving_dist_rat: 0.5 * abs(comoving_dist_rat ** 2.0 - 1.0 ** 2.0) / comoving_center_funct(ang_sep, comoving_dist_rat)

    ang_sep_mesh, dist_rat_mesh = np.meshgrid(ang_seps, comoving_dist_rats)
    comoving_perp_mesh = comoving_perps(ang_sep_mesh, dist_rat_mesh)
    comoving_parallel_mesh = comoving_parallels(ang_sep_mesh, dist_rat_mesh)
    
    if show_fig or save_fig: 
        matplotlib.rc('xtick', labelsize=14) 
        matplotlib.rc('ytick', labelsize=14) 
        f, axarr = plt.subplots(2,1, figsize = figsize)
        im0 = axarr[0].contourf(ang_sep_mesh, dist_rat_mesh, comoving_perp_mesh, levels = 21)
        im1 = axarr[1].contourf(ang_sep_mesh, dist_rat_mesh, comoving_parallel_mesh, levels = 21)
        axarr[0].set_xlabel(r'$\Delta \theta$ (rad)', fontsize = 18)
        axarr[1].set_xlabel(r'$\Delta \theta$ (rad)', fontsize = 18)
        axarr[0].set_ylabel(r'$S_2 / S_1 \ (\geq 1$ by def)', fontsize = 18)
        axarr[1].set_ylabel(r'$S_2 / S_1 \ (\geq 1$ by def)', fontsize = 18)
        axarr[0].set_title(r'$r_{\perp} / S_1$', fontsize = 18)
        axarr[1].set_title(r'$r_{\parallel} / S_1$', fontsize = 18)
        plt.tight_layout() 
        f.subplots_adjust(right = 0.8)
        cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
        f.colorbar(im1, cax = cbar_ax) 

    if show_fig:
        plt.show()

    if save_fig:
        f.savefig(save_dir + save_fig_name) 

    plt.close('all') 
    

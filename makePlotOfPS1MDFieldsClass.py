import numpy as np
import scipy
import cantrips as cant
import math
import numpy as np
import time
from PanStars1Archive import PanStars1Archive
from loadSN import loadSN
from binData import binData
import matplotlib.pyplot as plt
from DirectoryArchive import DirectoryArchive
#from FitStorer import FitStorer
from SDSSArchive import SDSSArchive
from AstronomicalParameterArchive import AstronomicalParameterArchive
from matplotlib.patches import Rectangle
import loadSN as lsn
import GeneralMCMC as mcmc
import scipy.optimize as optimize
import scipy
import randomSortData as rsd
import pickle
import sys

# import makePlotOfPS1MDFieldsClass as mpfc
# import cantrips as cant
# import numpy as np
# import matplotlib.pyplot as plt
# sdssdir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/SNIsotropyProject/SDSSGalaxies/'
# fulldata_fastRead = cant.readInColumnsToList('SDSS_fullCoverage_SDSSGals_pzAll.csv', sdssdir, delimiter = ',', n_ignore = 2, all_np_readable = 1)
# field_plotter = mpfc.PanStarsFieldManager(1, full_sdss_gal_data_file = 'SDSS_PSSuperField3_SDSSGals_Allpz.csv', preloaded_sdss_gals = fulldata_fastRead)

#Now we want to read in generated data.  Do these in python
"""
import makePlotOfPS1MDFieldsClass as mpfc
import cantrips as cant
import numpy as np
import matplotlib.pyplot as plt
# I may do this by field.  Just know that
data_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/SNIsotropyProject/SNeFieldFits/'
fields = [0,1,2,3,4,5,6,7,8,9]
chi_sqr_index = 6
overall_prob_index = 12
n_bins = 20
param_bounds = [(0.00, 3.0), (-2.0, 2.0), (-5.0, 5.0), (-5.0, 5.0)]
fit_param_indeces = [1,2,3,4]

rand_files_all_without_gal = ['RandFieldPS1MDFieldFitter_field0_1_2_3_4_5_6_7_8_9' + '_GalWeight0p0_N' + str(i) + '.csv' for i in list(range(0, 299))]
rand_files_all_with_gal = [['RandFieldPS1MDFieldFitter_field' + str(field) + '_GalWeight0p5_N' + str(i) + '.csv' for i in list(range(0, 101))] for field in fields]
loaded_data_all_without_gal = [cant.readInColumnsToList(data_dir + rand_file, delimiter = ',') for rand_file in rand_files_all_without_gal]
loaded_data_all_with_gal = [[cant.readInColumnsToList(data_dir + rand_file, delimiter = ',') for rand_file in rand_files] for rand_files in rand_files_all_with_gal]
good_rand_indeces_with_gal = [[np.all([ param_bounds[j][0] <= float(single_fit[fit_param_indeces[j]][fields[i]+1]) <= param_bounds[j][1] for j in range(len(param_bounds)) ]) for single_fit in loaded_data_all_with_gal[i]] for i in range(len(fields))]
loaded_rand_r_chi_sqrs_with_gal = [[float(single_fit[chi_sqr_index][fields[i]+1]) for single_fit in loaded_data_all_with_gal[i]] for i in range(len(fields))]
loaded_rand_overall_probs_with_gal = [[float(single_fit[overall_prob_index][fields[i]+1]) for single_fit in loaded_data_all_with_gal[i]] for i in range(len(fields))]
good_rand_indeces_without_gal = [[np.all([ param_bounds[j][0] <= float(single_fit[fit_param_indeces[j]][fields[i]+1]) <= param_bounds[j][1] for j in range(len(param_bounds)) ]) for single_fit in loaded_data_all_without_gal] for i in range(len(fields))]
loaded_rand_r_chi_sqrs_without_gal = [[float(single_fit[chi_sqr_index][fields[i]+1]) for single_fit in loaded_data_all_without_gal] for i in range(len(fields))]
loaded_rand_overall_probs_without_gal = [[float(single_fit[overall_prob_index][fields[i]+1]) for single_fit in loaded_data_all_without_gal] for i in range(len(fields))]
true_file_without_gal = 'TruePS1MDFieldFitter_field0_1_2_3_4_5_6_7_8_9' + '_GalWeight0p0_N' + str(0) + '.csv'
true_file_with_gal = 'TruePS1MDFieldFitter_field0_1_2_3_4_5_6_7_8_9' + '_GalWeight0p5_N' + str(0) + '.csv'
true_data_without_gal = cant.readInColumnsToList(data_dir + true_file_without_gal, delimiter = ',')
true_data_with_gal = cant.readInColumnsToList(data_dir + true_file_with_gal, delimiter = ',')
true_r_chi_sqrs_without_gal = [float(true_data_without_gal[chi_sqr_index][fields[i]+1]) for i in range(len(fields))]
true_r_chi_sqrs_with_gal = [float(true_data_with_gal[chi_sqr_index][fields[i]+1]) for i in range(len(fields))]
true_overall_probs_without_gal = [float(true_data_without_gal[overall_prob_index][fields[i]+1]) for i in range(len(fields))]
true_overall_probs_with_gal = [float(true_data_with_gal[overall_prob_index][fields[i]+1]) for i in range(len(fields))]
print ('good_rand_indeces_with_gal ')

exclude_out_of_field_fits = 1
show_without_gal = 1
if show_without_gal:
    plot_r_chi_sqr = 1
    loaded_rand_overall_probs = loaded_rand_overall_probs_without_gal
    loaded_rand_r_chi_sqrs = loaded_rand_r_chi_sqrs_without_gal
    true_overall_probs = true_overall_probs_without_gal
    true_r_chi_sqrs = true_r_chi_sqrs_without_gal
    good_indeces = good_rand_indeces_without_gal
else:
    plot_r_chi_sqr = 0
    loaded_rand_overall_probs = loaded_rand_overall_probs_with_gal
    loaded_rand_r_chi_sqrs = loaded_rand_r_chi_sqrs_with_gal
    true_overall_probs = true_overall_probs_with_gal
    true_r_chi_sqrs = true_r_chi_sqrs_with_gal
    good_indeces = good_rand_indeces_with_gal


n_hists_per_row = 5
fig_unit_len = 3


f, axarr = plt.subplots(int(np.ceil(len(fields) / n_hists_per_row)), n_hists_per_row  , figsize = (n_hists_per_row * fig_unit_len,   int(np.ceil(len(fields) / n_hists_per_row)) * fig_unit_len), squeeze = 'False')
for field_num in range(len(fields)):
    if len(fields) == 0:
        ax = axarr
    elif len(fields) <= n_hists_per_row:
        ax = axarr[field_num]
    else:
        ax = axarr[field_num // n_hists_per_row, field_num%  n_hists_per_row]
    field = fields[field_num]
    rand_r_chi_sqrs = loaded_rand_r_chi_sqrs[field_num]
    rand_overall_probs = loaded_rand_overall_probs[field_num]
    good_indeces_for_field = good_indeces[field_num]
    #print ('rand_res = ' + str(rand_res))
    true_r_chi_sqr = true_r_chi_sqrs[field_num]
    true_overall_prob = true_overall_probs[field_num]
    if plot_r_chi_sqr:
        if exclude_out_of_field_fits:
            rand_res = [ rand_r_chi_sqrs[i] for i in range(len(rand_r_chi_sqrs)) if good_indeces_for_field[i] ]
        else:
            rand_res = rand_r_chi_sqrs
        true_res = true_r_chi_sqr
        xlabel = r'$\chi_{\nu}^2$'
    else:
        if exclude_out_of_field_fits:
            rand_res = [rand_overall_probs[i] for i in range(len(rand_overall_probs)) if good_indeces_for_field[i]]
        else:
            rand_res = rand_overall_probs
        true_res = true_overall_prob
        xlabel = r'$P_{overall}$'
    hist_res = ax.hist(rand_res, bins = n_bins, fill = None)
    ax.axvline(x = true_res, c = 'cyan', ls = '-')
    better_rand = [ rand_res_elem for rand_res_elem in rand_res if rand_res_elem <= true_res]
    worse_rand = [ rand_res_elem for rand_res_elem in rand_res if rand_res_elem > true_res]
    #print('For field ' + str(field) + ', worse_chi_sqrs = ' + str(worse_chi_sqrs))
    n_better = len( better_rand )
    ax.text(true_res, max(hist_res[0]) / 2, str(n_better) + '/' + str(len(rand_res)) + ' rands <= true' , c = 'r', rotation = 90, horizontalalignment='right', verticalalignment='center')
    ax.text((hist_res[1][-1] - hist_res[1][0]) * 0.02 + hist_res[1][0], max(hist_res[0]) * 0.9, 'f' + str(field) + ': ' + xlabel + r', $_{T}$ = ' + str(cant.round_to_n(true_res, 3)), c = 'r', horizontalalignment='left', verticalalignment='center')
    ax.set_xlabel(xlabel )
    ax.set_ylabel(r'$\Delta \mu$ (mags)')

plt.tight_layout()
plt.show()
"""


class PanStarsFieldManager:

    #Basically, we need to invert the redshift to mu funct and then the vel to z funct
    def computePecVFromMuResids(self, full_zs, deltaMus):
        zs_from_mus = self.getPeculiarZsFromMuResids(full_zs, deltaMus)
        vs_from_zs = self.getLOSVsFromZs(zs_from_mus)
        #plt.scatter(deltaMus, zs_from_mus)
        #plt.xlabel(r'$\Delta \mu$')
        #plt.ylabel(r'$z(\mu)$')
        #plt.show()
        return vs_from_zs

    def getPeculiarZsFromMuResids(self, ref_zs, deltaMus):
        r_of_z_interp = self.r_of_z_interp
        ref_zs = np.array(ref_zs)
        deltaMus = np.array(deltaMus)
        print('[ref_zs, deltaMus] = ' + str([ref_zs, deltaMus]) )
        if r_of_z_interp is None:
            print ('Have not implemented this yet')
            zs_from_mus = [0.0 for mu in deltaMus]
        else:
            #delta_mu = 5 * np.log10( (1 + zf) * r_of_z_interp((zf - zg)/(zg + 1)) ) - 5 * np.log10( (1 + zf) * r_of_z_interp(zf) )
            #         = 5 log10(r_of_z_interp((zf - zg)/(zg + 1)) / r_of_z_interp(zf))
            intermed_term = (10.0 ** (deltaMus / 5.0) - 1.0) * np.array([scipy.integrate.quad(lambda z_int: 1/self.H_of_z(z_int), 0, ref_z)[0] for ref_z in ref_zs]) * (-self.H_of_z(ref_zs) / (1.0 + ref_zs))
            print ('intermed_term = ' + str(intermed_term))
            zs_from_mus = intermed_term / (1.0- intermed_term)

        return zs_from_mus


    def getLOSVsFromZs(self, zs):
        one_p_zs = 1.0 + np.array(zs)
        vs = (one_p_zs ** 2.0 - 1.0) / (1.0 + one_p_zs ** 2.0)
        return vs

    def getErroneousRedshiftPlot(self, z_extra_space, z_obs_space,
                                 zHD = 1, z_include_range = [0.0, np.inf], surveys_of_interest = ['all'], surveys_to_excise = [''],
                                 ):

        r_of_z_interp = self.r_of_z_interp
        c = self.astro_arch.getc() # in km/s

        #Note H0 is given in km/s/Mpc.  c is in km/s.  So c/H0 gives this ratio in Mpc.
        if r_of_z_interp is None:
            dLInMpc = lambda zf, zg: c / H0 * (1 + zf) * scipy.integrate.quad(lambda z_int: 1/np.sqrt((1 + z_int) ** 3 * OmM + (1 + z_int) ** 0 * OmL + (1 + z_int) ** 4 * OmR), 0, (zf - zg)/(zg + 1))[0]
        else:
            #dLInMpc = lambda zf, zg: c / H0 * (1 + zf) * r_of_z_interp((zf - zg)/(zg + 1))
            dLInMpc = lambda zf, zg: (1 + zf) * r_of_z_interp((zf - zg)/(zg + 1))
        mu =  lambda zf, zg: 5 * np.log10(dLInMpc(zf, zg)) + 25
        #vals = [mu(z, zg_of_z_funct(z)) - mu(z, 0.0) for z in z_space]
        vals = [mu(z_obs_space[z_index], z_extra_space[z_index]) - mu(z_obs_space[z_index], 0.0) for z_index in range(len(z_obs_space))]
        vals = np.where(np.isnan(vals), np.inf, vals)
        #for z_index in range(len(z_obs_space)):
        #    if dLInMpc(z_obs_space[z_index], z_extra_space[z_index]) == 0.0:
        #        print ('[z_obs_space[z_index], z_extra_space[z_index], dLInMpc(z_obs_space[z_index], z_extra_space[z_index]), mu(z_obs_space[z_index], z_extra_space[z_index]), vals[z_index] ] = ' + str([z_obs_space[z_index], z_extra_space[z_index], dLInMpc(z_obs_space[z_index], z_extra_space[z_index]), mu(z_obs_space[z_index], z_extra_space[z_index]), vals[z_index] ]))

        return vals

    def plot_mwd(self, RA, Dec, ax, org=0,title='Aitoff Sky', projection='aitoff'):
        ''' RA, Dec are arrays of the same length.
        RA takes values in [0,360), Dec in [-90,90],
        which represent angles in degrees.
        org is the origin of the plot, 0 or a multiple of 30 degrees in [0,360).
        title is the title of the figure.
        projection is the kind of projection: 'mollweide', 'aitoff', 'hammer', 'lambert'
        '''
        x = np.remainder(RA+360-org,360) # shift RA values
        ind = x>180
        x[ind] -=360    # scale conversion to [-180, 180]
        x=-x    # reverse the scale: East to the left
        tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
        tick_labels = np.remainder(tick_labels+360+org,360)
        #fig = plt.figure(figsize=(10, 5))
        #ax = fig.add_subplot(111, projection=projection, axisbg ='LightCyan')
        ax.scatter(np.radians(x),np.radians(Dec))  # convert degrees to radians
        ax.set_xticklabels(tick_labels)     # we add the scale on the x axis
        ax.set_title(title)
        ax.title.set_fontsize(15)
        ax.set_xlabel("RA")
        ax.xaxis.label.set_fontsize(12)
        ax.set_ylabel("Dec")
        ax.yaxis.label.set_fontsize(12)
        ax.grid(True)

    def plotPantheonSNOnSky(self, data_set,
                            projection = 'aitoff', save = 1, show = 0, surveys_to_display = 'all', z_range = [-0.1, 3.0],
                            pull_extinctions = 0, figsize = [10.0, 8.0], zHD = 1):

        dir_archive = DirectoryArchive()
        plot_dir = dir_archive.getPlotDirectory()
        deg_to_rad = self.astro_arch.getDegToRad()

        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(2,1)
        sky_plot = fig.add_subplot(gs[0,0], projection="aitoff")
        sky_plot.grid(True)
        sky_plot.set_xlabel('R.A.')
        sky_plot.set_ylabel('Decl.')
        tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
        tick_labels = np.remainder(tick_labels+360+0,360)

        #fig = plt.figure(figsize=(10, 5))
        #ax = fig.add_subplot(111, projection=projection, axisbg ='LightCyan')
        sky_plot.set_xticklabels(tick_labels)     # we add the scale on the x axis
        all_sns = loadSN(1, ['all'], pull_extinctions = pull_extinctions, zHD = zHD)
        all_sns = [sn for sn in all_sns if (sn['z'] >= z_range[0] and sn['z'] <= z_range[1]) ]
        if surveys_to_display in ['all','All','ALL']:
            surveys_to_display = [sn['survey'] for sn in all_sns]
            surveys_to_display = np.unique(surveys_to_display).tolist()
        sn_by_survey = [[sn for sn in all_sns if sn['survey'] == survey] for survey in surveys_to_display if len([sn for sn in all_sns if sn['survey'] == survey]) > 0]
        #print ('sn_by_survey = '  +str(sn_by_survey))
        colors = [sn_set[0]['color'] for sn_set in sn_by_survey]
        plots_for_legend = [0 for  sn_in_survey in sn_by_survey]
        for i in range(len(plots_for_legend)):
            sn_in_survey = sn_by_survey[i]
            RA = np.array([(sn['RA']) for sn in sn_in_survey])
            Dec = [(sn['Dec']) for sn in sn_in_survey]
            x = np.remainder(RA+360-0,360) # shift RA values
            ind = x>180
            x[ind] -=360    # scale conversion to [-180, 180]
            x=-x    # reverse the scale: East to the left
            plots_for_legend[i] = sky_plot.scatter(np.radians(x), np.radians(Dec) , color = [sn['color']  for sn in sn_in_survey], marker = 'o', s = 2.0 )
        sky_plot.legend(plots_for_legend, surveys_to_display, ncol = 3, fontsize = 8, loc='lower right', bbox_to_anchor=(1.1, -0.05))
        redshift_plot = fig.add_subplot(gs[1,0])
        for i in range(len(sn_by_survey)):
             sn_in_survey = sn_by_survey[i]
             redshift_plot.scatter([sn['z'] for sn in sn_in_survey], [sn['muDiff'] for sn in sn_in_survey] , c = [sn['color'] for sn in sn_in_survey])
        if save:
            plt.savefig(plot_dir + 'SN_on_sky_' +str(z_range[0]) + '<' + str(z_range[1]) + '.pdf')
        if show:
            plt.show()

    def printingFitFunct(self, fit_funct, zs, muDiffs, muErrs, params, print_params = 0):
        fit_vals = fit_funct(zs, *params)
        diffs_from_fit = (np.array(fit_vals) - np.array(muDiffs))
        weighted_mean = np.sum(diffs_from_fit * np.array(muErrs) ** (-2.0)) / np.sum( np.array(muErrs) ** (-2.0))
        chi_sqr = np.sum((diffs_from_fit - weighted_mean) ** 2.0 / np.array(muErrs) ** 2.0)
        if print_params: print ('np.array(params).tolist() + [chi_sqr] = ' +str(np.array(params).tolist() + [chi_sqr]) )
        return fit_vals, chi_sqr, weighted_mean

    def getFittedMCMCObject(self, fit_funct, all_zs, all_resids, all_errs, dof, mcmc_start_params, mcmc_step_sizes, bounds, ):
        MCMC_fit_funct = lambda params: np.sum(((fit_funct(all_zs, *params) - all_resids) / all_errs) ** 2.0 )

        new_MCMC_fitter = mcmc.MCMCFitter(MCMC_fit_funct, mcmc_start_params, mcmc_step_sizes, self.n_mcmc_steps, bounds = bounds, likelihood_from_chisqr  = 1)
        return new_MCMC_fitter

    def doMinimization(self, fit_funct, zs, resids, errs, dof, fit_params_A, fit_params_B, fit_params_C, fit_params_D = None, fit_params_E = None, bounds= None):
        print('[fit_params_A, fit_params_B, fit_params_C, fit_params_D, fit_params_E] = ' + str([fit_params_A, fit_params_B, fit_params_C, fit_params_D, fit_params_E] ))
        brute_minimization = []
        brute_chisqr = np.inf
        start = time.time()
        if fit_params_D is None:
            fit_params_D = [None]
        if fit_params_E is None:
            fit_params_E = [None]
        print ('[fit_params_A, fit_params_B, fit_params_C, dof, fit_params_D, fit_params_E] = ' + str([fit_params_A, fit_params_B, fit_params_C, dof, fit_params_D, fit_params_E] ))
        for fit_param_A in fit_params_A:
            print ('working on fit_param_A = ' + str(fit_param_A))
            for fit_param_B in fit_params_B:
                for fit_param_C in fit_params_C:
                    for fit_param_D in fit_params_D:
                        #this_fit = [fit_z, fit_beta, fit_rwidth, fit_power]
                        for fit_param_E in fit_params_E:
                            this_fit = [fit_param_A, fit_param_B, fit_param_C]
                            if not(fit_param_D is None) :
                                this_fit = this_fit + [fit_param_D]
                            if not(fit_param_E is None) :
                                this_fit = this_fit + [fit_param_E]
                            #print ('[fit_param_A, fit_param_B, fit_param_C, dof, fit_param_D, fit_param_E] = ' + str([fit_param_A, fit_param_B, fit_param_C, dof, fit_param_D, fit_param_E] ))
                            #this_fit_resids = fit_funct(zs, *this_fit)
                            #this_fit_mean = np.sum( ((this_fit_resids - resids) * errs ** (-2.0))) / np.sum(errs ** (-2.0))
                            #this_fit_mean = c.weighted_mean(this_fit_resids - resids, errs)
                            #print ('this_fit_mean = ' + str(this_fit_mean ))
                            #print ('this_fit = ' + str(this_fit))
                            this_fit_resids, this_fit_chisqr, this_fit_weighted_mean = self.printingFitFunct(fit_funct, zs, resids, errs, np.array(this_fit))
                            #print ('this_fit + [this_fit_chisqr] = ' + str(this_fit + [this_fit_chisqr]))
                            #this_fit_rchisqr_argmin = np.argmin(this_fit_rchisqrs)
                            if this_fit_chisqr < brute_chisqr:
                                brute_chisqr = this_fit_chisqr
                                brute_minimization = this_fit[:]
        end = time.time()
        print ('Brute minimization took ' + str(end - start) + 's')
        print ('[brute_minimization, brute_chisqr] = ' + str([brute_minimization, brute_chisqr]))
        print ('Now trying to refine using curve_fit...')
        true_fit_rchisqr = brute_chisqr
        funct_fit = brute_minimization
        try:
            print ('brute_minimization = ' + str(brute_minimization))
            #final_minimization = scipy.optimize.curve_fit(fit_funct, all_zs, all_resids, p0 = brute_minimization, sigma=all_errs, maxfev = maxfev, bounds = bounds)
            print ('[bounds, brute_minimization] = ' + str([bounds, brute_minimization]))
            final_minimization = scipy.optimize.minimize(lambda params: self.printingFitFunct(fit_funct, zs, resids, errs, params, print_params = 0)[1], brute_minimization, bounds , method='L-BFGS-B')
            #fit_res = optimize.curve_fit(lambda zs, A, mu: fit_funct(zs, A, mu, sig, 0.2, 0.0), all_zs, all_resids, p0 = [0.0, max(all_zs) / 2.0], sigma=all_errs, maxfev = maxfev)
            this_funct_fit = final_minimization['x']
        except RuntimeError:
            print("Curve_fit failed!.  Plotting initial guess. ")
            this_funct_fit = np.array(brute_minimization)
        final_chisqr = self.printingFitFunct(fit_funct, zs, resids, errs, np.array(this_funct_fit))[1]
        print ('[this_funct_fit, final_chisqr] = ' + str([this_funct_fit, final_chisqr]))
        final_rchisqr = final_chisqr / dof

        return [this_funct_fit, final_chisqr]


    def redshift_of_v_funct(self, beta, coses_of_los_angle, max_beta = 0.999):
        beta = np.where(np.abs(beta) < 1.0, beta, max_beta * np.sign(beta))
        redshift = ((1 + beta * coses_of_los_angle)/np.sqrt(1-abs(beta) ** 2.0) - 1)
        return redshift

    def computeVelocityFieldFromNFWHalo(self, zs_to_calc, central_z, critical_mass, concentration, impact_param,
                                        r_of_z_interp, ):

        #extra_z_of_zMeas_funct = lambda single_z, central_z, v, r_width: redshift_of_v_funct(v) if abs( r_of_z_interp(central_z) - r_of_z_interp(single_z) ) < r_width / 2.0 else 0.0
        #D_growth = lambda single_z: H_of_z(single_z) * integrate.quad(lambda z_int: (1 + z_int) / np.sqrt((1 + z_int) ** 3 * OmM + (1 + z_int) ** 0 * OmL + (1 + z_int) ** 4 * OmR) ** 3.0, single_z, np.inf)[0]
        #d_z_D_growth = lambda single_z (1.0 + )
        d_z_D_growth_over_D_term1 = lambda single_z: d_H_d_z (single_z) / H_of_z(single_z)
        d_z_D_growth_over_D_term2 = lambda single_z: -(1.0 + single_z) / H_of_z(single_z) ** 3.0 / scipy.integrate.quad(lambda z_int: (1 + z_int) / H_of_z(z_int) ** 3.0, single_z, np.inf)[0]
        #d_z_D_growth_over_D = lambda single_z: (1.0 + single_z) * H_of_z(single_z) * (d_z_D_growth_over_D_term1(single_z) + d_z_D_growth_over_D_term2(single_z))
        d_z_D_growth_over_D = lambda single_z: (d_z_D_growth_over_D_term1(single_z) + d_z_D_growth_over_D_term2(single_z))
        #The geometry of how our line of sight intersects the halo
        r_incidence_leg = lambda single_r, central_r, r_incidence: (np.sqrt(central_r ** 2.0 - r_incidence ** 2.0 ) - single_r)
        r_well_from_geometry = lambda single_r, central_r, r_incidence: np.sqrt(r_incidence ** 2.0 + r_incidence_leg(single_r, central_r, r_incidence) ** 2.0)
        sin_r_angle_from_geometry = lambda single_r, central_r, r_incidence:  r_incidence_leg(single_r, central_r, r_incidence) / r_well_from_geometry(single_r, central_r, r_incidence)

        delta_vir_of_z = lambda single_z: (18.0 * np.pi + 82.0 * ((Om0 *(1.0 + single_z) ** 3.0) / H_of_z(single_z) ** 2.0) - 39.0 * ((Om0 *(1.0 + single_z) ** 3.0) / H_of_z(single_z) ** 2.0) ** 2.0)
        #for NFW; r_s should be given in Mpc
        nfw_unitless_int = lambda s: np.log(s + 1) - s / (s+1)
        nfw_mass_enclosed = lambda single_r, concentration, critical_mass, critical_r: nfw_unitless_int(concentration * single_r / critical_r) / nfw_unitless_int(concentration) * critical_mass
        int_of_rsqr_delta_nfw = lambda scaled_r, c: (nfw_unitless_int(scaled_r) - nfw_unitless_int(0.0) - 1.0 / 3.0 * (scaled_r) ** 3.0 / (c * (1 + c) ** 2.0))
        #delta_portion_funct = lambda single_z, central_z, delta_s, rs, c, b_impact: rs ** 1.0 * (rs / r_well_from_geometry(r_of_z_interp(single_z), r_of_z_interp(central_z), b_impact * rs)) ** 2.0 * delta_s * (int_of_rsqr_delta_nfw( min(c, r_well_from_geometry(r_of_z_interp(single_z), r_of_z_interp(central_z), b_impact * rs) / rs), c) ) if (single_z != central_z or abs(b_impact) > 0.0) else 0.0
        delta_portion_funct = lambda single_z, central_z, delta_s, rs, c, b_impact: rs ** 1.0 * (rs / r_well_from_geometry(r_of_z_interp(single_z), r_of_z_interp(central_z), b_impact * r_of_z_interp(central_z))) ** 2.0 * delta_s * (int_of_rsqr_delta_nfw( min(c, r_well_from_geometry(r_of_z_interp(single_z), r_of_z_interp(central_z), b_impact * rs) / rs), c) ) if (single_z != central_z or abs(b_impact) > 0.0) else 0.0



    def redshiftDueToGravPotential (self, zs, central_z, delta, R):
        speedol = self.astro_arch.getc()
        well_depth = -delta * self.OmM * self.H0 ** 2.0 * R ** 2.0 / 2.0
        central_r = self.r_of_z_interp(central_z)
        phi_of_delta = [well_depth * (2 - (abs(self.r_of_z_interp(single_z) - central_r) / R) ** 2.0 if abs(self.r_of_z_interp(single_z) - central_r) < R else R / abs(self.r_of_z_interp(single_z) - central_r)) for single_z in zs]
        extra_z_of_zMeas_funct_G = (phi_of_delta(single_z, float(self.r_of_z_interp(central_z)), delta, R) - phi_of_delta(0.0, float(self.r_of_z_interp(central_z)), delta, R)) / (speedol) ** 2.0
        fitted_G = np.array(self.getErroneousRedshiftPlot(lambda single_z: extra_z_of_zMeas_funct_G(single_z, central_z, delta, R), zs, zHD = self.zHD, astro_arch = self.astro_arch, r_of_z_interp = self.r_of_z_interp ))
        return fitted_G


    def getFittingFunctG(self):
        n_G_fit_params = 3
        G_bounds = [(0.0, 4.0), (-1.0, 10.0), (50.0, 5000.0)]
        fit_central_zs = np.linspace(0.0, 1.0, 15)
        fit_deltas = np.linspace(-1.0, 1.0, 6)
        fit_Rs = np.linspace(100.0, 1000.0, 6) #in Mpc
        fit_funct = self.redshiftDueToGravPotential

        return [fit_funct, n_G_fit_params, G_bounds, [fit_central_zs, fit_deltas, fit_Rs]]

    def d_z_D_growth_over_D(self, single_z):
        d_z_D_growth_over_D_term1 = self.d_H_d_z (single_z) / self.H_of_z(single_z)
        d_z_D_growth_over_D_term2 = -(1.0 + single_z) / self.H_of_z(single_z) ** 3.0 / scipy.integrate.quad(lambda z_int: (1 + z_int) / self.H_of_z(z_int) ** 3.0, single_z, np.inf)[0]
        d_z_D_growth_over_D = d_z_D_growth_over_D_term1 + d_z_D_growth_over_D_term2
        return d_z_D_growth_over_D


    #The geometry of how our line of sight intersects the halo
    def r_incidence_leg (self, single_r, central_r, r_incidence):
        leg =  (np.sqrt(central_r ** 2.0 - r_incidence ** 2.0 ) - single_r)
        if np.isnan(leg).any():
            print ('[single_r, central_r, r_incidence, leg] = ' + str([single_r, central_r, r_incidence, leg]  ))
        return leg

    def r_well_from_geometry(self, single_r, central_r, r_incidence ):
        r_well =  np.sqrt(r_incidence ** 2.0 + self.r_incidence_leg(single_r, central_r, r_incidence) ** 2.0)
        return r_well

    #This seems wrong?
    def sin_of_infall_angle_relative_to_los(self, single_r, central_r, r_incidence):
        sin_r = self.r_incidence_leg(single_r, central_r, r_incidence) / self.r_well_from_geometry(single_r, central_r, r_incidence)
        return sin_r

    def cos_of_infall_angle_relative_to_los(self, single_r, central_r, r_incidence):
        cos_r = self.r_incidence_leg(single_r, central_r, r_incidence) / self.r_well_from_geometry(single_r, central_r, r_incidence)
        return cos_r

    def overdensty_vir_of_z (self, single_z):
        overdensity = (18.0 * np.pi + 82.0 * ((self.Om0 *(1.0 + single_z) ** 3.0) / H_of_z(single_z) ** 2.0) - 39.0 * ((self.Om0 *(1.0 + single_z) ** 3.0) / H_of_z(single_z) ** 2.0) ** 2.0)
        return overdensity

    def nfw_unitless_over_int(self, radius_per_rs):
        nfw_int = np.log(radius_per_rs + 1) - radius_per_rs / (radius_per_rs+1)
        return nfw_int

    def nfw_mass_enclosed(self, r, mass_overdensity, concentration, critical_density, overdensity):
        r_overdensity = self.r_overdensity(mass_overdensity, critical_density, overdensity)
        M_enc = self.nfw_unitless_over_int(concentration * r / r_overdensity) / self.nfw_unitless_over_int(concentration) * mass_overdensity

        return M_enc

    def nfw_mass_enclosed_ours(self, halo_density, radius_in_scale_radii):
        #r_overdensity = self.r_overdensity(mass_overdensity, critical_density, overdensity)
        #M_enc = self.nfw_unitless_over_int(concentration * r / r_overdensity) / self.nfw_unitless_over_int(concentration) * mass_overdensity
        M_enc = halo_density * (np.log(radius_in_scale_radii + 1.0) + 1.0 / (1.0 + radius_in_scale_radii))
        return M_enc

    def point_mass_enclosed(self, halo_mass):
        #r_overdensity = self.r_overdensity(mass_overdensity, critical_density, overdensity)
        #M_enc = self.nfw_unitless_over_int(concentration * r / r_overdensity) / self.nfw_unitless_over_int(concentration) * mass_overdensity
        return halo_mass

    def r_overdensity(self, M_overdensity, critical_density, overdensity):
        r = ((M_overdensity * 3.0) / (4.0 * np.pi * overdensity * critical_density )) ** (1.0 / 3.0)
        return r

    # returns the critical mass with units of T M_{sun} Mpc ^ -3
    def crit_density(self, z):
        G = self.astro_arch.getGravitationalConstant()
        crit_density_unitless_part = 3.0 * self.H_of_z(z) ** 2.0 / (8.0 * np.pi)
        pc_to_m = self.astro_arch.getParsecToM()
        Mpc_to_km = pc_to_m  * 10.0 ** 3.0
        Msun_to_kg = self.astro_arch.getSolarMass()
        crit_density_unitfull_part = (self.H0 / Mpc_to_km) ** 2.0 / G / (10.0 ** 17.0 * Msun_to_kg) * (pc_to_m * 10.0 ** 6.0) ** 3.0
        crit_density = crit_density_unitless_part * crit_density_unitfull_part
        #print ('[crit_density_unitless_part, crit_density_unitfull_part] = ' + str([crit_density_unitless_part, crit_density_unitfull_part]))
        return crit_density

    #Mass should be in M_sun X 10 ^(12) ; returns r ** 2.0 * int(delta r*2 4 pi dr) in units of Mpc
    def delta_portion_funct(self, single_z, central_z, M_overdensity, concentration, b_impact, overdensity, ):
        central_r = self.r_of_z_interp(central_z) #units of Mpc
        critical_density = self.crit_density(single_z) #units of GMsun / Mpc ^ 3
        r_overdensity = self.r_overdensity(M_overdensity, critical_density, overdensity)
        rs = r_overdensity / concentration
        r_well = self.r_well_from_geometry(self.r_of_z_interp(single_z), self.r_of_z_interp(central_z), b_impact * central_r)

        M_enc = self.nfw_mass_enclosed(r_well, M_overdensity, concentration, critical_density, overdensity)

        delta_over_rsqr_int_over_rsqr = 1.0 / (4.0 * np.pi * r_well ** 2.0) * M_enc / critical_density

        return delta_over_rsqr_int_over_rsqr

    #Mass should be in M_sun X 10 ^(12) ; returns r ** 2.0 * int(delta r*2 4 pi dr) in units of Mpc
    # delta_portion_funct_ours(single_z, central_z, halo_density, r_scale, b_impact, overdensity)
    def delta_portion_funct_ours(self, single_z, central_z, halo_mass, r_scale, r_incidence  ):
        single_r = self.r_of_z_interp(single_z)
        central_r = self.r_of_z_interp(central_z) #units of Mpc
        critical_density = self.crit_density(single_z) #units of GMsun / Mpc ^ 3
        #r_overdensity = self.r_overdensity(M_overdensity, critical_density, overdensity)
        #rs = r_overdensity / concentration
        r_well = self.r_well_from_geometry(single_r, central_r, r_incidence)

        #M_enc = self.nfw_mass_enclosed_ours(halo_density, r_well / r_scale )
        #M_enc_scaled = M_enc * r_scale ** 3.0

        M_enc = self.point_mass_enclosed(halo_mass)
        M_enc_scaled = M_enc / (self.OmM * critical_density)

        delta_over_rsqr_int_over_rsqr = - 1.0 / (4.0 * np.pi * r_well ** 2.0) * M_enc_scaled

        return delta_over_rsqr_int_over_rsqr

    def computeVirialOverdensity(self, single_z):
        local_OmM = self.OmM * (1.0 + single_z) ** 3.0 / (self.OmL + self.OmM * (1.0 + single_z) ** 3.0)
        local_x = local_OmM - 1
        delta_vir = 18.0 * np.pi ** 2.0 + 82.0 * local_x - 39.0 * local_x ** 2.0
        return delta_vir

    def getLOSVelocityFieldAlongLOSTraditional(self, zs, central_z, M_overdensity, concentration, b_impact, overdensity = 'virial'):
        zs = calc_zs + ref_zs
        if overdensity == 'virial':
            overdensity = np.array([self.computeVirialOverdensity(single_z) for single_z in zs])
        Hs_of_z = [self.H_of_z(single_z) for single_z in zs]
        d_zs_D_over_D = [self.d_z_D_growth_over_D(single_z) for single_z in zs]


        single_rs = [self.r_of_z_interp(single_z) for single_z in zs]
        central_r = self.r_of_z_interp(central_z)
        #critical_densities = np.array([self.crit_density(single_z) for single_z in zs])
        critical_density = self.crit_density(central_z)

        r_overdensity = self.r_overdensity(M_overdensity, critical_density, overdensity)
        r_scale = r_overdensity / concentration
        sins_of_los_angle = self.sin_of_infall_angle_relative_to_los(single_rs, central_r, central_r * b_impact)

        delta_portions = [self.delta_portion_funct(single_z, central_z, M_overdensity, concentration, b_impact, overdensity) for single_z in zs]

        speedol = self.astro_arch.getc()
        los_velocity_field = self.H0 / speedol * np.array(Hs_of_z) * np.array(d_zs_D_over_D) * np.array(sins_of_los_angle) * delta_portions

        return los_velocity_field


    #parameters are the mass in one scale radius, the scale radius, the central z, and angles describing the impact parameter
    #def getLOSVelocityFieldAlongLOSOurs(self, zs, RAs, Decs, field_center, central_z, halo_scaling_power, comoving_scale_radius_power, central_sky_phi, central_sky_theta, ):
    def getLOSVelocityFieldAlongLOSOurs(self, zs, RAs, Decs, field_center, central_z, halo_mass_power, central_sky_RA_offset, central_sky_Dec_offset, ):
        #For now, we are going to see what happens when we limit our model to just having the halo_scaling_power adjust the halo nature -- basically making it a point mass
        comoving_scale_radius_power = self.fixed_comoving_scale_radius_power
        deg_to_rad = self.astro_arch.getDegToRad()
        field_RA, field_Dec = field_center
        #central_RA, central_Dec = [field_RA + np.sin(central_sky_phi * deg_to_rad) * central_sky_theta * np.cos(field_Dec * deg_to_rad), field_Dec + np.cos(central_sky_phi * deg_to_rad) * central_sky_theta]
        central_RA, central_Dec = [field_RA + central_sky_RA_offset, field_Dec + central_sky_Dec_offset]
        #print ('[central_RA, central_Dec, central_sky_phi, central_sky_theta] = ' + str([central_RA, central_Dec, central_sky_phi, central_sky_theta]))
        #print ('field_center = ' + str(field_center))
        Hs_of_z = [self.H_of_z(single_z) for single_z in zs]
        d_zs_D_over_D = [self.d_z_D_growth_over_D(single_z) for single_z in zs]

        single_rs = np.array([self.r_of_z_interp(single_z) for single_z in zs])
        central_r = self.r_of_z_interp(central_z)
        #critical_densities = np.array([self.crit_density(single_z) for single_z in zs])
        critical_density = self.crit_density(central_z)
        r_scale = 10.0 ** comoving_scale_radius_power #scale radius in Mpc
        halo_mass = np.sign(halo_mass_power) * (10.0 ** abs(halo_mass_power) - 1.0) #Allow the sampling to go from -whatever to + whatever, passing through 0 when halo_mass_power = 0
        #vel_at_scale_radius = np.sign(vel_at_scale_radius_power) * (10.0 ** (-abs(vel_at_scale_radius_power) - 1.0)
        #rolled_v_scale_power = np.where(vel_at_scale_radius_power > self.min_v_scale_power, vel_at_scale_radius_power, -vel_at_scale_radius_power + 2.0 * self.min_v_scale_power) #Rolls over to correspond to negative velocities at the minimium allowed power
        #vel_at_scale_radius = np.where(vel_at_scale_radius_power > self.min_v_scale_power, 1.0, -1.0) * (10.0 ** rolled_v_scale_power  - 10.0 ** self.min_v_scale_power) #self.max_v_scale_power is the given power at which the velocity at z = 0 is c
        #speedol = self.astro_arch.getc()
        #halo_scaling  = 4.0 * np.pi * vel_at_scale_radius / (r_scale * self.H0 / speedol * self.nfw_mass_enclosed_ours(1.0, 1.0))
        #halo_scaling = mass_scaling / (self.nfw_mass_enclosed_ours(1.0, 1.0) * r_scale ** 3.0)
        #print ('[vel_at_scale_radius_power, vel_at_scale_radius, halo_scaling, r_scale] = ' + str([vel_at_scale_radius_power, vel_at_scale_radius, halo_scaling, r_scale]))
        angular_offsets = np.array([cant.measureAngularSeparationOnSky([RAs[i], Decs[i]], [central_RA, central_Dec], return_radian = 1) for i in range(len(zs))])

        #r_offsets = np.sqrt(single_rs ** 2.0 + central_r ** 2.0 - 2.0 * single_rs * central_r * np.cos(angular_offsets * deg_to_rad))
        r_incidences = np.sin(angular_offsets) * central_r
        #print ('angular_offsets.tolist() = ' + str(angular_offsets.tolist()))
        #print ('r_incidences.tolist() = ' + str(r_incidences.tolist()))
        coses_of_los_angle = self.cos_of_infall_angle_relative_to_los(single_rs, central_r, r_incidences)
        #print ('[zs, central_z, halo_scaling, r_scale, r_incidences] = ' + str([zs, central_z, halo_scaling, r_scale, r_incidences] ))
        delta_portions = np.array([self.delta_portion_funct_ours(zs[i], central_z, halo_mass, r_scale, r_incidences[i]) for i in range(len(zs))]) #in units of Mpc

        speedol = self.astro_arch.getc()
        #self.H0 / speedol gives H0 / c in units of 1/MPC
        #print ('[np.isnan(Hs_of_z).any(), np.isnan(d_zs_D_over_D).any(), np.isnan(sins_of_los_angle).any(), np.isnan(delta_portions).any()] = ' + str([np.isnan(Hs_of_z).any(), np.isnan(d_zs_D_over_D).any(), np.isnan(sins_of_los_angle).any(), np.isnan(delta_portions).any()]))
        #print ('[np.isinf(Hs_of_z).any(), np.isinf(d_zs_D_over_D).any(), np.isinf(sins_of_los_angle).any(), np.isinf(delta_portions).any()] = ' + str([np.isinf(Hs_of_z).any(), np.isinf(d_zs_D_over_D).any(), np.isinf(sins_of_los_angle).any(), np.isinf(delta_portions).any()]))

        #los_velocity_field = self.H0 / speedol * np.array(Hs_of_z) * np.array(d_zs_D_over_D) * np.array(sins_of_los_angle) * delta_portions
        los_velocity_field = self.H0 / speedol * np.array(Hs_of_z) * np.array(d_zs_D_over_D) * delta_portions
        #print ('los_velocity_field = ' + str(los_velocity_field))

        return los_velocity_field, coses_of_los_angle

    #Here
    def getMuDiffOfVelocityField(self, calc_zs, calc_RAs, calc_Decs, ref_zs, ref_RAs, ref_Decs, ref_muDiffs, ref_muErrs, field_center, vel_field_params, vel_field_funct):
        all_zs = np.array(calc_zs).tolist() + np.array(ref_zs).tolist()
        all_RAs = np.array(calc_RAs).tolist() + np.array(ref_RAs).tolist()
        all_Decs = np.array(calc_Decs).tolist() + np.array(ref_Decs).tolist()

        scaled_velocities, coses_of_los_angles = vel_field_funct(all_zs, all_RAs, all_Decs, field_center, *vel_field_params)
        vel_redshifts = self.redshift_of_v_funct(scaled_velocities, coses_of_los_angles)

        #for z_index in range(len(all_zs)):
        #    if vel_redshifts[z_index] > all_zs[z_index]:
        #        print ('[all_zs[z_index], all_RAs[z_index], all_Decs[z_index], scaled_velocities[z_index], coses_of_los_angles[z_index], vel_redshifts[z_index]] = ' + str([all_zs[z_index], all_RAs[z_index], all_Decs[z_index], scaled_velocities[z_index], coses_of_los_angles[z_index], vel_redshifts[z_index]]  ))
        #print('[vel_redshifts, all_zs] = ' + str([vel_redshifts, all_zs]))
        erroneous_mus = self.getErroneousRedshiftPlot(vel_redshifts, all_zs,)
        #print ('erroneous_mus = ' + str(erroneous_mus))
        #print ('[all_zs, scaled_velocities, vel_redshifts, erroneous_mus] = ' + str([all_zs, scaled_velocities, vel_redshifts, erroneous_mus]))
        #f1, axarr = plt.subplots(4,1)
        #axarr[0].scatter(scaled_velocities, vel_redshifts)
        #axarr[1].scatter(all_zs, scaled_velocities, marker = '.')
        #axarr[2].scatter(all_zs, vel_redshifts, marker = '.')
        #axarr[3].scatter(all_zs, erroneous_mus, marker = '.')
        #plt.show()
        ref_calced_muDiffs = erroneous_mus[-len(ref_zs):]
        weights = np.array(ref_muErrs) ** (-2.0)
        #print ('weights = ' + str(weights))
        #print ('np.isinf(weights).any() = ' + str(np.isinf(weights).any()))
        #print ('ref_muDiffs = ' + str(ref_muDiffs))
        #print ('np.isinf(ref_muDiffs).any() = ' + str(np.isinf(ref_muDiffs).any()))
        #print ('ref_calced_muDiffs = ' + str(ref_calced_muDiffs))
        #print ('np.isinf(ref_calced_muDiffs).any() = ' + str(np.isinf(ref_calced_muDiffs).any()))
        weighted_diffs = (np.array(ref_calced_muDiffs) - np.array(ref_muDiffs) ) * weights
        bad_diffs = np.where(np.isinf(weighted_diffs), 0.0, 1.0)
        weighted_mean_diff = np.sum(np.where(bad_diffs, weighted_diffs, 0.0)) / np.sum(weights * bad_diffs)
        corrected_muDiffs = np.array(erroneous_mus[0:-len(ref_zs)]) - weighted_mean_diff
        return corrected_muDiffs


    def getFittingFunctVelNFW(self):
        n_vel_fit_params = 4
        vel_bounds = [(0.0, 2.0), (1, 1000.0), (1, 20), (0, 5)]
        fit_central_zs = np.linspace(0.01, 1.0, 21)
        fit_M_overdensities = np.linspace(1, 1000.0, 5) #In T M_{sun}
        fit_concentrations = np.linspace(1, 20, 3)
        fit_impact_params = np.linspace(0.1, 2.0, 4)
        fit_overdensities = [200]
        fit_funct = self.getMuDiffOfVelocityField

        return [fit_funct, n_vel_fit_params, vel_bounds, [fit_central_zs, fit_central_zs, fit_M_overdensities, fit_concentrations]]

    #def getFittingFunctVelNFWOurs(self, n_seeds_in_z = 50, n_seeds_around_sky = 4):
    def getFittingFunctVelNFWOurs(self, seed_z_step = 0.01, n_seeds_around_sky = 4):
        n_vel_fit_params = 4
        #Note: for scaling (second) and scale radius (third) params, the numbers are logarithmic
        #vel_bounds = [(0.0, 1.0), (0.0, 10.0), (0.5, 4.0), (0.0, 360.0), (0.01, 10.0)]
        #vel_bounds = [(0.0, 1.0), (0.0, 5.0), (0.0, 360.0), (0.01, 10.0)]
        vel_bounds = [[0.01, self.z_range[1]], (0.0, 2.0), (-5.0, 5.0), (-5.0, 5.0)]
        #vel_bounds = [[0.28, 0.32]]
        #fit_central_zs = np.linspace(0.0, 1.0, 21)
        #fit_halo_scaling_powers = np.linspace(0.1, 4.0, 9) #In T M_{sun}
        #fit_r_scale_powers = np.linspace(0.1, 3, 10)
        #fit_impact_params = np.linspace(0.1, 3.0, 6)
        #fit_overdensities = [200]
        #mcmc_step_sizes = [0.005, 0.02, 0.015, 2.0, 0.05]
        mcmc_step_sizes = [0.005, 0.1, 0.05, 0.05]
        sub_funct = self.getLOSVelocityFieldAlongLOSOurs
        muDiff_of_z_funct = self.getMuDiffOfVelocityField
        MCMC_chain_starts  = []
        seed_zs = np.arange(*vel_bounds[0], seed_z_step)
        for seed_z in seed_zs:
            for j in range(n_seeds_around_sky):
                MCMC_chain_starts = MCMC_chain_starts  + [[seed_z, 0.0 * np.log10(2.0), np.cos(2.0 * np.pi / n_seeds_around_sky * (0.5 + j)) * 1.0, np.sin(2.0 * np.pi / n_seeds_around_sky * (0.5 + j)) * 1.0]]
        #MCMC_chain_starts = MCMC_chain_starts  + [[vel_bounds[0][1] / (n_seeds_in_z - 1) * i + vel_bounds[0][0], np.log10(2.0), 0.05, 0.05] for i in range(n_seeds_in_z)]
        #print ('MCMC_chain_starts = ' + str(MCMC_chain_starts))
        return [sub_funct, muDiff_of_z_funct, n_vel_fit_params, vel_bounds, MCMC_chain_starts, mcmc_step_sizes]

    def computeRChiSqr(self, all_zs, all_RAs, all_Decs, all_resids, all_errs, field_center, n_model_params, fit_params, print_params = 0):
        if print_params: print ('fit_params = ' + str(fit_params))
        fitted_resids = self.muDiff_of_z_funct(all_zs, all_RAs, all_Decs, all_zs, all_RAs, all_Decs, all_resids, all_errs, field_center, fit_params, self.fit_funct)
        #weights = np.array(all_errs) ** (-2.0)
        #weighted_mean_diff = np.sum(np.array(fitted_resids - all_resids)  * weights) / np.sum(weights)
        #weighted_mean_diff = 0.0
        #print ('weighted_mean_diff = ' + str(weighted_mean_diff))
        #print('fitted_resids = ' + str(fitted_resids))
        total_dof = len(all_zs) - n_model_params
        chi_sqr = np.sum(((np.array(fitted_resids)) - all_resids) ** 2.0 / (np.array(all_errs) ** 2.0))
        r_chi_sqr = chi_sqr / total_dof
        if print_params: print ('np.array(fit_params).tolist() + [r_chi_sqr] = ' + str(np.array(fit_params).tolist()+ [r_chi_sqr]))
        return r_chi_sqr

    def getSDSSGalDensity(self, fit_params, field_center, comoving_sample_rad, bg_comoving_annulus_inner_rad, bg_comoving_annulus_outer_rad):
        center_redshift = fit_params[0]
        center_RA, center_Dec = [fit_params[2] + field_center[0], fit_params[3] + field_center[1]]
        if center_Dec > 180.0:
            center_Dec = 180.0 - (center_Dec - 180.0)
            center_RA = center_RA + 180.0
        if center_Dec < 180.0:
            center_Dec =  - (center_Dec)
            center_RA = center_RA + 180.0
        center_in_footprint = self.archive.checkIfSkyCoordsInFootprint([center_RA], [center_Dec])
        if center_in_footprint:
            relative_gal_n = self.archive.measureRelativeGalDensityAroundCoord(center_redshift, center_RA, center_Dec, comoving_sample_rad, bg_comoving_annulus_inner_rad, bg_comoving_annulus_outer_rad, verbose = 0)
        else:
            relative_gal_n = 1

        return relative_gal_n


    def getSNeDataForSurveysInField(self, field_num):
        surveys_to_include = self.surveys_to_include
        included_surveys = self.included_surveys
        sn_in_survey_by_field = self.sn_in_survey_by_field
        all_zs = [[] for i in range(len(surveys_to_include))]
        all_resids = [[] for i in range(len(surveys_to_include))]
        all_errs = [[] for i in range(len(surveys_to_include))]
        all_RAs = [[] for i in range(len(surveys_to_include))]
        all_Decs = [[] for i in range(len(surveys_to_include))]
        print ('field_num = ' + str(field_num))
        for i in range(len(included_surveys)):
            #print ('[len(sn_in_survey_by_field), len(sn_in_survey_by_field[i]), i, field_num] = ' + str([len(sn_in_survey_by_field), len(sn_in_survey_by_field[i]), i, field_num]))
            sns_in_field_in_survey = sn_in_survey_by_field[i][field_num]
            zs = [sn['z'] for sn in sns_in_field_in_survey]
            dls = [sn['dl'] for sn in sns_in_field_in_survey]
            rs = [dls[j] / (1.0 + zs[j]) for j in range(len(zs))]
            all_zs[i] =  zs
            #print ('all_zs = ' + str(all_zs))
            muDiffs = [sn['muDiff'] - sn['muDiffWMean'] for sn in sns_in_field_in_survey]
            all_resids[i] = muDiffs
            muErrs = [sn['muErr'] for sn in sn_in_survey_by_field[i][field_num]]
            all_errs[i] = muErrs
            RAs = [sn['RA'] for sn in sn_in_survey_by_field[i][field_num]]
            Decs = [sn['Dec'] for sn in sn_in_survey_by_field[i][field_num]]
            all_RAs[i] = RAs
            all_Decs[i] = Decs
        return [all_zs, all_resids, all_errs, all_RAs, all_Decs]


    def makePlotOfPS1MDFields(self, fields_to_plot = 'all', plot_fit = 0,
                              plot_delta_vs = 0,
                              save = 1, show = 0, fit_information = {'funct':'none'}, show_fit_label = 1, n_plot_zs = 1001, z_range_to_plot = [-0.1, 3.0],
                              mu_plot_lims = [-0.6, 0.6], vel_plot_lims = [-0.07, 0.07], z_plot_lims = [0.0, 1.0], archive_to_use = 'PS1MD' , figsize = [16.0, 8.0], save_name = 'SN_residuals_v_z_PS1MD_fields_.pdf',
                              xlabel = r'$z$', ylabel = r'$\Delta \mu$ (mags)',
                              n_z_bins = 10, super_plot_title = r'$\Delta \mu$ of Pan Stars 1 medium-deep fields', plot_param_round_to = 3, plot_chi_sqr_round_to = 3, min_fig_side = 3 ):

        master_start = time.time()
        dir_archive = DirectoryArchive()
        plot_dir = dir_archive.getPlotDirectory()
        deg_to_rad = self.astro_arch.getDegToRad()
        pc_to_m = self.astro_arch.getParsecToM()
        Mpc_to_km = pc_to_m * 10.0 ** 6.0 * 10.0 ** -3.0
        speedol = self.astro_arch.getc()
        zs_to_plot = np.linspace(self.min_z, self.max_z, n_plot_zs)


        if fields_to_plot == 'all':
            fields = self.fields
            field_nums = list(fields.keys())
        else:
            fields = {field_key:self.fields[field_key] for field_key in self.fields.keys() if field_key in fields_to_plot}
            field_nums = fields_to_plot[:]

        chi_squares_by_field = {field:[-1, -1] for field in fields}

        colors = [sn_set[0]['color'] for sn_set in self.sn_by_survey]

        sn_in_survey_by_field = self.sn_in_survey_by_field
        sn_by_survey = self.sn_by_survey

        n_fields = len(fields)
        fig_side = 0
        while fig_side * 4 - 4 < n_fields:
            fig_side = fig_side + 1
        fig_side = max(fig_side, min_fig_side )

        fig = plt.figure(constrained_layout=True, figsize = figsize)

        gs = fig.add_gridspec(fig_side, fig_side)
        plot_indeces = [[fig_side - 1, i] for i in range(fig_side)] + [[0, i] for i in range(fig_side)] + cant.flattenListOfLists( [[[i, 0], [i, fig_side - 1]] for i in range(1, fig_side-1)] )
        for plot_index in plot_indeces:
            fig.add_subplot(gs[plot_index[0], plot_index[1]])
        #Add the central sky plot
        sky_plot = fig.add_subplot(gs[1:fig_side-1, 1:fig_side-1], projection="aitoff")
        sky_plot.grid(True)
        sky_plot.set_xlabel('R.A.')
        sky_plot.set_ylabel('Decl.')
        tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
        tick_labels = np.remainder(tick_labels+360+0,360)
        #fig = plt.figure(figsize=(10, 5))
        #ax = fig.add_subplot(111, projection=projection, axisbg ='LightCyan')
        sky_plot.set_xticklabels(tick_labels)     # we add the scale on the x axis
        plots_for_legend = [0 for  sn_in_survey in sn_by_survey]
        for i in range(len(plots_for_legend)):
            sn_in_survey = sn_by_survey[i]
            RA = np.array([(sn['RA']) for sn in sn_in_survey])
            Dec = [(sn['Dec']) for sn in sn_in_survey]
            x = np.remainder(RA+360-0,360) # shift RA values
            ind = x>180
            x[ind] -=360    # scale conversion to [-180, 180]
            x=-x    # reverse the scale: East to the left
            plots_for_legend[i] = sky_plot.scatter(np.radians(x), np.radians(Dec) , color = [sn['color']  for sn in sn_in_survey], marker = 'o', s = 2.0 )
        legend_plot_index = plot_indeces[-1]
        ax = fig.add_subplot(gs[legend_plot_index[0], legend_plot_index[1]])
        ax.legend(plots_for_legend, self.included_surveys, ncol = 3, fontsize = 8)

        for j in range(len(field_nums)):
            field_num = field_nums[j]
        #for field_num in range(1):
            field_str = 'f' + str(field_num)
            field = fields[field_num]
            print ('field = ' + str(field))
            field_center = [(field[0] + field[1]) / 2.0, (field[2] + field[3]) / 2.0]
            plot_index = plot_indeces[j]
            ax = fig.add_subplot(gs[plot_index[0], plot_index[1]])
            rect_coords1 =[ -(field[1] if field[1] < 180.0 else field[1] - 360), field[2]]
            rect_coords2 =[ -(field[0] if field[0] < 180.0 else field[0] - 360) - rect_coords1[0] , field[3] - field[2]]
            field_rect = Rectangle([rect_coord * deg_to_rad for rect_coord in rect_coords1], *[rect_coord * deg_to_rad for rect_coord in rect_coords2], edgecolor='black', facecolor='black', alpha = 0.25)
            sky_plot.add_patch(field_rect)
            sky_plot.text((rect_coords1[0] + rect_coords2[0]) * deg_to_rad, (rect_coords1[1] + rect_coords2[1]) * deg_to_rad, field_str)

            print ('Displaying field number ' + str(field_num))
            survey_plots = []
            fitted_plots = []
            fit_position_index = 0
            all_zs = []
            all_RAs = []
            all_Decs = []
            all_dls = []
            all_resids = []
            all_errs = []
            #fitters = [ [FitStorer(fit_information) for survey in self.surveys_to_include ] for field in fields]
            zs_by_survey, resids_by_survey, mu_errs_by_survey, ras_by_survey, decs_by_survey = self.getSNeDataForSurveysInField(field_num)
            for i in range(len(self.included_surveys)):
                color = colors[i]
                if len(zs_by_survey[i]) > 0:
                    zs_data = zs_by_survey[i]
                    resids_data = resids_by_survey[i]
                    mu_errs_data = mu_errs_by_survey[i]
                    if plot_delta_vs:
                        vs_data, v_errs_data = [self.computePecVFromMuResids(zs_data, resids_data), self.computePecVFromMuResids(zs_data, mu_errs_data)]
                        xs_to_plot, ys_to_plot, y_errs_to_plot = [zs_data, vs_data, v_errs_data]
                        #print ('[xs_to_plot, ys_to_plot, y_errs_to_plot] = ' + str([xs_to_plot, ys_to_plot, y_errs_to_plot]))
                    else:
                        xs_to_plot, ys_to_plot, y_errs_to_plot = [zs_data, resids_data, mu_errs_data]
                    survey_plots = survey_plots + [ax.scatter(xs_to_plot, ys_to_plot, c = color) ]
                    ax.errorbar(xs_to_plot, ys_to_plot, yerr = y_errs_to_plot, ecolor = color, fmt = 'none')
                else:
                    survey_plots = survey_plots + [ax.scatter([],[])]

            all_zs = cant.flattenListOfLists(zs_by_survey)
            all_RAs = cant.flattenListOfLists(ras_by_survey)
            all_Decs = cant.flattenListOfLists(decs_by_survey)
            all_resids  = cant.flattenListOfLists(resids_by_survey)
            all_errs  = cant.flattenListOfLists(mu_errs_by_survey)

            label_str = field_str
            if plot_fit and len(self.fit_params_by_field[field_num]) > 0:
                fit_params = self.fit_params_by_field[field_num]
                null_params = self.null_params
                RAs_to_plot =  [field_center[0] for z in zs_to_plot]
                Decs_to_plot = [field_center[1] for z in zs_to_plot]
                print ('fit_params = ' + str(fit_params))
                muDiffs_to_plot = self.muDiff_of_z_funct(zs_to_plot, RAs_to_plot, Decs_to_plot, all_zs, all_RAs, all_Decs, all_resids, all_errs, field_center, fit_params, self.fit_funct)
                all_muDiffs = self.muDiff_of_z_funct(all_zs, all_RAs, all_Decs, all_zs, all_RAs, all_Decs, all_resids, all_errs, field_center, fit_params, self.fit_funct)
                if plot_delta_vs:
                    #xs_to_plot, ys_to_plot = [zs_to_plot, muDiffs_to_plot]
                    xs_to_plot, ys_to_plot = [zs_to_plot, self.computePecVFromMuResids(zs_to_plot, muDiffs_to_plot)]
                    #all_xs, all_ys = [all_zs, all_muDiffs]
                    all_xs, all_ys = [all_zs, self.computePecVFromMuResids(all_zs, all_muDiffs)]
                else:
                    xs_to_plot, ys_to_plot = [zs_to_plot, muDiffs_to_plot]
                    all_xs, all_ys = [all_zs, all_muDiffs]
                ax.plot(xs_to_plot, ys_to_plot, c = 'k')
                ax.scatter(all_xs,  all_ys, c = 'k', marker = 'x')
                chi_sqr = self.computeRChiSqr(all_zs, all_RAs, all_Decs, all_resids, all_errs, field_center, self.n_fit_params + 1, fit_params)
                null_chi_sqr = self.computeRChiSqr(all_zs, all_RAs, all_Decs, all_resids, all_errs, field_center, 1, null_params)
                print ('chi_sqr = ' + str(chi_sqr))
                label_str = label_str + 'params ' + str( [cant.round_to_n(fit_param, plot_param_round_to) for fit_param in fit_params]) + r' with $\chi^2_{\nu}$ ' + str(cant.round_to_n(chi_sqr, plot_chi_sqr_round_to)) + r'; $\chi^2_{\nu,0}$ ' + str(cant.round_to_n(null_chi_sqr, plot_chi_sqr_round_to))
                chi_squares_by_field[field_num] = [chi_sqr, null_chi_sqr ]

            print ('label_str = ' + str(label_str))
            ax.set_xlim(z_plot_lims)
            if plot_delta_vs:
                ax.set_ylim(vel_plot_lims)
            else:
                ax.set_ylim(mu_plot_lims)

            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.text(0.5, 0.8, label_str,  horizontalalignment='center', transform=ax.transAxes, fontsize = 7)


        plt.suptitle(super_plot_title)

        if save: plt.savefig(plot_dir + save_name )
        if show: plt.show()
        plt.close('all')

        return chi_squares_by_field

    #Note, we return H normalized by H_0, usually denoted E(z)
    def H_of_z(self, z_for_H):
        H = np.sqrt((1 + z_for_H) ** 3 * self.OmM + (1 + z_for_H) ** 0 * self.OmL + (1 + z_for_H) ** 4 * self.OmR)
        return H
    #Note, we return H normalized by H_0, usually denoted E(z)
    def d_H_d_z(self, z_for_H):
        d_H = 0.5 * 1.0 / self.H_of_z(z_for_H) * (3.0 * (1 + z_for_H) ** 2.0 * self.OmM + 4.0 * (1 + z_for_H) ** 3.0 * self.OmR)
        return d_H

    def r_of_z(self, z_meas):
        speedol = self.astro_arch.getc()
        r = (speedol) / self.H0 * scipy.integrate.quad( lambda z_int: 1.0 / self.H_of_z(z_int), 0, z_meas )[0]
        return r

    def initialize_r_of_z_interp(self, interp_z_params):
        self.interp_zs = np.linspace(interp_z_params[0], interp_z_params[1], interp_z_params[2]).tolist() + [1000.0]
        self.r_of_z_interp = scipy.interpolate.interp1d([-1000.0] + self.interp_zs, [0.0] + [self.r_of_z(z) for z in self.interp_zs], fill_value = (0.0, self.r_of_z(1000)), bounds_error = False)
        return 1

    #Randomly pair the SNe residuals and residual uncertainties with zs,
    def randomizeSNBySurvey(self):
        all_surveys = list(set([sn['survey'] for sn in self.all_sns]))
        randomized_sn = []
        sn_by_survey = {survey:[] for survey in all_surveys}
        for sn in self.all_sns:
            sn_by_survey[sn['survey']] = sn_by_survey[sn['survey']] + [sn]
        for survey in all_surveys:
            sn_in_survey = sn_by_survey[survey]
            all_muDiffs = [sn['muDiff'] for sn in sn_in_survey]
            all_muErrs = [sn['muErr'] for sn in sn_in_survey]
            rand_muDiffs, rand_muErrs = rsd.randomShuffleListOfLists([all_muDiffs, all_muErrs])
            for j in range(len(sn_in_survey)):
                sn_in_survey[j]['muDiffTrue'] = sn_in_survey[j]['muDiff']
                sn_in_survey[j]['muDiff'] = rand_muDiffs[j]
                sn_in_survey[j]['muErrTrue'] = sn_in_survey[j]['muErr']
                sn_in_survey[j]['muErr'] = rand_muErrs[j]
            randomized_sn = randomized_sn + sn_in_survey
        self.all_sn = randomized_sn

    def initializeSN(self, z_range, surveys_to_include, surveys_to_ignore):

        all_sns = loadSN(self.data_set, surveys_to_include, pull_extinctions = self.pull_extinctions, zHD = self.zHD, OmM = self.OmM, OmL = self.OmL, Om0 = self.Om0, OmR = self.OmR, H0 = self.H0)
        if surveys_to_include[0] == 'all':
            surveys_to_include = [sn['survey'] for sn in all_sns]
            surveys_to_include = list(set(surveys_to_include))
        self.surveys_to_include = surveys_to_include
        self.surveys_to_ignore = surveys_to_ignore
        all_sns = [sn for sn in all_sns if not(sn['survey'] in self.surveys_to_ignore)]
        all_sns = [sn for sn in all_sns if (sn['z'] >= z_range[0] and sn['z'] <= z_range[1]) ]
        all_zs = [sn['z'] for sn in all_sns]
        self.all_sns = all_sns
        self.all_zs = all_zs
        self.min_z, self.max_z = [np.min(self.all_zs), np.max(self.all_zs)]
        return 1

    def randomlySortMusInFields(self):
        included_surveys = self.included_surveys
        sn_in_survey_by_field = self.sn_in_survey_by_field
        #print ('self.sn_in_survey_by_field = ' + str(sn_in_survey_by_field))
        #print ('self.sn_in_survey_by_field[1][0][0] = ' + str(sn_in_survey_by_field[1][0][0]))
        for field_key in self.fields:
            for i in range(len(included_surveys)):
                sns_in_field_in_survey = sn_in_survey_by_field[i][field_key]
                true_muDiffs = [sn['muDiff'] for sn in sns_in_field_in_survey]
                true_muErrs = [sn['muErr'] for sn in sns_in_field_in_survey]
                #print ('[true_muDiffs, true_muErrs] = ' + str([true_muDiffs, true_muErrs]))
                if len(sns_in_field_in_survey) > 0:
                    rand_muDiffs, rand_muErrs = rsd.randomShuffleListOfLists([true_muDiffs, true_muErrs])

                for j in range(len(sns_in_field_in_survey)):
                    sns_in_field_in_survey[j]['muErr'] = rand_muErrs[j]
                    sns_in_field_in_survey[j]['muErrTrue'] = true_muErrs[j]
                    sns_in_field_in_survey[j]['muDiff'] = rand_muDiffs[j] + sns_in_field_in_survey[j]['muDiffWMean']
                    sns_in_field_in_survey[j]['muDiffTrue'] = true_muDiffs[j] + sns_in_field_in_survey[j]['muDiffWMean']
                self.sn_in_survey_by_field[i][field_key] = sns_in_field_in_survey
        #print ('self.sn_in_survey_by_field[1][0][0] = ' + str(sn_in_survey_by_field[1][0][0]))

        return 1

    def initializeSNByField(self, archive_to_use, randomize = 0 ):

        PSArch = PanStars1Archive(OmM = self.OmM, OmL = self.OmL, OmR = self.OmR, H0 = self.H0, n_healpix_sides = self.archive_healpix_sides, full_sdss_gal_data_file = self.full_sdss_gal_data_file, preloaded_sdss_gals = self.preloaded_sdss_gals)
        SDSSArch = SDSSArchive()

        if self.archive_to_use.lower() == 'sdss':
            self.archive = SDSSArch
        else:
            self.archive = PSArch
        self.fields = self.archive.fields
        self.field_centers = self.archive.field_centers
        print('self.fields = ' + str(self.fields))

        print ('self.surveys_to_include = ' + str(self.surveys_to_include ) )
        self.sn_by_survey = [[sn for sn in self.all_sns if sn['survey'] == survey] for survey in self.surveys_to_include if len([sn for sn in self.all_sns if sn['survey'] == survey]) > 0]
        self.included_surveys = [sn_by_single_survey[0]['survey'] for sn_by_single_survey in self.sn_by_survey]
        print ('self.included_surveys = ' + str(self.included_surveys ))
        #print ('sn_by_survey = '  +str(sn_by_survey))
        colors = [sn_set[0]['color'] for sn_set in self.sn_by_survey]

        sn_in_survey_by_field = []

        for survey_index in range(len(self.sn_by_survey)):
            sn_set = self.sn_by_survey[survey_index]
            print ('Working on survey ' + self.included_surveys[survey_index])
            sn_by_field = {}
            for key in self.fields:
                field = self.fields[key]
                #sn_by_field[key] = [sn for sn in sn_set if (sn['RA']>field[0] and sn['RA'] < field[1] and sn['Dec'] > field[2] and sn['Dec'] < field[3] )]
                sn_by_field[key] = [sn for sn in sn_set if (sn['RA']>field[0] and sn['RA'] < field[1] and sn['Dec'] > field[2] and sn['Dec'] < field[3] )]
                #print ('For [field, key] ' + str([field, key]) + ', len(sn_by_field[key]) = ' + str(len(sn_by_field[key])))
            sn_in_survey_by_field = sn_in_survey_by_field + [sn_by_field]
        self.sn_in_survey_by_field = sn_in_survey_by_field

        print('randomize = ' + str(randomize))
        if randomize:
            self.randomlySortMusInFields()

        return 1


    def computeMCMCFits(self, field_nums, n_mcmc_steps = None, mcmc_step_sizes = None, n_chains = 'all'):
        if n_mcmc_steps == None:
            n_mcmc_steps = self.n_mcmc_steps
        if mcmc_step_sizes == None:
            mcmc_step_sizes = self.mcmc_step_sizes
        if n_chains == 'all':
            mcmc_start_params = self.mcmc_start_params[:]
        else:
            mcmc_start_params = self.mcmc_start_params[0:n_chains]
        MCMC_fitters_by_field = [0 for field_num in field_nums]

        for i in range(len(field_nums)):
            field_num = field_nums[i]
            all_zs_in_field_by_survey, all_muResids_in_field_by_survey, all_muErrs_in_field_by_survey, all_RAs_in_field_by_survey, all_Decs_in_field_by_survey = self.getSNeDataForSurveysInField(field_num)
            all_zs = cant.flattenListOfLists(all_zs_in_field_by_survey)
            all_muResids = cant.flattenListOfLists(all_muResids_in_field_by_survey)
            all_muErrs = cant.flattenListOfLists(all_muErrs_in_field_by_survey)
            MCMC_fit_funct = lambda params: self.computeRChiSqr(all_zs, all_muResids, all_muErrs, params)
            new_MCMC_fitter = mcmc.MCMCFitter(MCMC_fit_funct, mcmc_start_params, mcmc_step_sizes, n_mcmc_steps, bounds = self.bounds, likelihood_from_chisqr  = 1)
            null_chi_sqr = MCMC_fit_funct(self.null_params)
            MCMC_fitters_by_field [field_num] = {'fit':new_MCMC_fitter, 'null_chi_square': null_chi_sqr}
        self.MCMC_fits = MCMC_fitters_by_field

        return 1

    #This is the cumulative density function
    def computeChiSqrProb(self, all_zs, all_RAs, all_Decs, all_muResids, all_muErrs, field_center, n_fit_params, params, print_params = 0):
        dof = len(all_zs) - n_fit_params
        chi_sqr = self.computeRChiSqr(all_zs, all_RAs, all_Decs, all_muResids, all_muErrs, field_center, n_fit_params, params, print_params = print_params) * dof
        chi_sqr_prob = scipy.special.gammainc(dof / 2.0, chi_sqr / 2.0) # chi_sqr ** ((dof - 2) / 2) * np.exp(-chi_sqr / 2.0)
        #print ('[chi_sqr, chi_sqr_prob] = ' + str([chi_sqr, chi_sqr_prob]))
        return chi_sqr_prob

    #This is the cumulative density function
    #def computeGalDensProb(self, targ_z, targ_RA, targ_Dec, inner_comoving_rad, annulus_inner_comoving_rad = None, annulus_outer_comoving_rad = None):
    def computeGalDensProb(self, targ_z, targ_RA, targ_Dec, halo_excess_mass, annulus_inner_comoving_rad = None, annulus_outer_comoving_rad = None):
        #if inner_comoving_rad_frac is None:
        #    inner_comoving_rad_frac = self.inner_comoving_rad_frac
        #inner_comoving_rad = inner_comoving_rad_frac * self.r_of_z_interp(targ_z)
        inner_comoving_rad = self.r_of_z_interp(targ_z) * self.start_spherical_comoving_rad_frac
        #print ('inner_comoving_rad = ' + str(inner_comoving_rad))
        if annulus_inner_comoving_rad is None:
            annulus_inner_comoving_rad = self.annulus_inner_comoving_rad_scaling * inner_comoving_rad
        if annulus_outer_comoving_rad is None:
            annulus_outer_comoving_rad = self.annulus_outer_comoving_rad_scaling * inner_comoving_rad
        obj_center = [targ_RA, targ_Dec]
        n_gal_density, n_gal_density_err = self.archive.measureRelativeGalDensityAroundCoord(targ_z, *obj_center, inner_comoving_rad, annulus_inner_comoving_rad, annulus_outer_comoving_rad, verbose = 1)
        #print ('[n_gal_density, n_gal_density_err] = ' + str([n_gal_density, n_gal_density_err]))
        if halo_excess_mass > 0.0:
            normal_gal_dens_prob = 1.0 - 0.5 * (1.0 + scipy.special.erf((n_gal_density - 1.0) / (n_gal_density_err * np.sqrt(2.0))))
        else:
            normal_gal_dens_prob = 0.5 * (1.0 + scipy.special.erf((n_gal_density - 1.0) / (n_gal_density_err * np.sqrt(2.0))))

        return n_gal_density, n_gal_density_err, normal_gal_dens_prob

    def overall_prob_funct(self, all_zs, all_RAs, all_Decs, all_muResids, all_muErrs, field_center, ext_params, print_params = 1, gal_dens_weight = 0, n_model_params = None) :
        #print('ext_params = ' + str(ext_params))
        if n_model_params is None:
            n_model_params = self.n_fit_params
        sne_fit = self.computeChiSqrProb(all_zs, all_RAs, all_Decs, all_muResids, all_muErrs, field_center, n_model_params + 1, ext_params[0:4], print_params = 0) * (1.0 - gal_dens_weight)
        #print ('sne_prob = ' + str(sne_prob))
        #gal_dens_prob = self.computeGalDensProb(ext_params[0], ext_params[2] + field_center[0], ext_params[3] + field_center[1], ext_params[4])
        sphere_center = [(ext_params[2] + field_center[0]), (ext_params[3] + field_center[1])]
        if sphere_center[1] > 90.0:
            sphere_center[1] = sphere_center[1] - 90.0
            n_loops = (sphere_center[1] ) // 180.0 + 1
            sphere_center[0] = sphere_center[0] + 180.0 * n_loops
            sphere_center[1] = sphere_center[1] % 360.0
            if sphere_center[1] > 180.0: sphere_center[1] = 360.0 - sphere_center[1]
            sphere_center[1] = 90 - sphere_center[1]
        elif sphere_center[1] <  -90.0:
            sphere_center[1] = sphere_center[1] + 90.0
            n_loops = (sphere_center[1]) // -180.0 + 1
            sphere_center[0] = sphere_center[0] + 180.0 * n_loops
            sphere_center[1] = sphere_center[1] % -360.0
            if sphere_center[1] < 180.0: sphere_center[1]: sphere_center[1] = -360.0 - sphere_center[1]
            sphere_center[1] = -90 - sphere_center[1]

        if gal_dens_weight > 0.0:
            gal_dens_fit = self.computeGalDensProb(ext_params[0], sphere_center[0], sphere_center[1], ext_params[1])[-1] * gal_dens_weight
        else:
            gal_dens_fit = 0.0
        #print ('gal_dens_prob = ' + str(gal_dens_prob))
        overall_fit = sne_fit + gal_dens_fit
        #overall_prob = sne_probf
        if print_params:
            print ('[ext_params, sne_fit, gal_dens_fit, overall_fit] = ' + str([ext_params, sne_fit, gal_dens_fit, overall_fit]))
        return overall_fit

    ##method =  'trust_constr'
    def computeByHandFits(self, field_nums, n_mcmc_steps = None, mcmc_step_sizes = None, n_chains = 'all', method = 'Nelder-Mead', add_to_master_min_dict = 1, update_field_fits_as_best_fit = 1, update_null_chi_by_field = 1, start_comoving_rad_frac = None,  frac_seeds_to_search_for_gals = 0.1):
        if start_comoving_rad_frac is None:
            start_comoving_rad_frac = self.start_spherical_comoving_rad_frac
        master_start = time.time()
        if n_mcmc_steps == None:
            n_mcmc_steps = self.n_mcmc_steps
        if mcmc_step_sizes == None:
            mcmc_step_sizes = self.mcmc_step_sizes
        if n_chains == 'all':
            mcmc_start_params = self.mcmc_start_params[:]
        else:
            mcmc_start_params = self.mcmc_start_params[0:n_chains]

        min_res_by_field = {}
        r_chi_square_by_field = {}
        null_chi_square_by_field = {}
        fit_funct_no_gal_by_field = {}
        fit_funct_with_gal_by_field = {}
        n_sig_gal_by_field = {}
        prob_of_gal_by_field = {}
        overall_prob_of_gal_by_field = {}
        #For now, I want to fix rs to be its starting value
        for i in range(len(field_nums)):
            field_num = field_nums[i]
            all_zs_in_field_by_survey, all_muResids_in_field_by_survey, all_muErrs_in_field_by_survey, all_RAs_in_field_by_survey, all_Decs_in_field_by_survey = self.getSNeDataForSurveysInField(field_num)
            all_zs = cant.flattenListOfLists(all_zs_in_field_by_survey)
            all_muResids = cant.flattenListOfLists(all_muResids_in_field_by_survey)
            all_muErrs = cant.flattenListOfLists(all_muErrs_in_field_by_survey)
            all_RAs = cant.flattenListOfLists(all_RAs_in_field_by_survey)
            all_Decs = cant.flattenListOfLists(all_Decs_in_field_by_survey)
            field = self.fields[field_num]
            field_center = self.field_centers[field_num]
            print ('!!!! Starting working on [field_num, field, field_center] = ' + str([field_num, field, field_center]))
            fit_funct_no_gal_dens = lambda ext_params: self.overall_prob_funct(all_zs, all_RAs, all_Decs, all_muResids, all_muErrs, field_center, ext_params, print_params = 0, gal_dens_weight = 0.0)
            null_fit_funct_no_gal_dens = lambda ext_params: self.overall_prob_funct(all_zs, all_RAs, all_Decs, all_muResids, all_muErrs, field_center, ext_params, print_params = 0, gal_dens_weight = 0.0, n_model_params = 0)
            fit_funct_with_gal_dens = lambda ext_params: self.overall_prob_funct(all_zs, all_RAs, all_Decs, all_muResids, all_muErrs, field_center, ext_params, print_params = 0, gal_dens_weight = self.gal_dens_weighting )
            no_gal_min_res = [0 for seed in range(len(mcmc_start_params)) ]
            for min_seed_index in range(len(mcmc_start_params)):
                seed_start = time.time()
                start_min_seed = mcmc_start_params[min_seed_index]
                #start_comoving_rad = start_comoving_rad_frac * self.r_of_z_interp(max(min_seed[0], 0.001))
                #min_seed = min_seed + [start_comoving_rad]
                if (start_min_seed[0] > min(all_zs) and start_min_seed[0] < max(all_zs)):
                    print ('Starting minimization with min_seed: ' + str(start_min_seed))
                    no_gal_new_min_res = optimize.minimize(fit_funct_no_gal_dens, start_min_seed, method = method, tol = 10.0 ** -5.0, bounds = self.bounds)
                    #print ('new_min_res = ' + str(new_min_res))
                    print ('Found minimizing params without including galaxy density: ' + str(no_gal_new_min_res['x']) + ' with chi_sqr, chi_sqr likelihood: ' + str(self.computeRChiSqr(all_zs, all_RAs, all_Decs, all_muResids, all_muErrs, field_center, self.n_fit_params + 1, no_gal_new_min_res['x'], print_params = 0)) + ', ' + str(fit_funct_no_gal_dens(no_gal_new_min_res['x']))  )
                    no_gal_min_res[min_seed_index] = [no_gal_new_min_res['x'], no_gal_new_min_res['fun']]
                    #if self.gal_dens_weighting > 0.0:
                    #    print ('Starting minimization with that min_seed, now including galaxy density: ' + str(no_gal_new_min_res['x']))
                    #    new_min_res = optimize.minimize(fit_funct_with_gal_dens, no_gal_new_min_res['x'], method = method, tol = 10.0 ** -4.0, bounds = self.bounds)
                    #    print ('Found minimizing params with including galaxy density: ' + str(new_min_res['x']) + ' with chi_sqr_likelihood, gal_n_dens_likelihood, and overall likelihoods: ' + str([fit_funct_no_gal_dens(new_min_res['x']) * (1.0 - self.gal_dens_weighting), fit_funct_with_gal_dens(new_min_res['x']) - fit_funct_no_gal_dens(new_min_res['x']) * (1.0 - self.gal_dens_weighting) , fit_funct_with_gal_dens(new_min_res['x'])]))
                    #else:
                    #    new_min_res = no_gal_new_min_res
                    #min_res[min_seed_index] = [new_min_res['x'], new_min_res['fun']]
                else:
                    print ('Not going to do fit for seed ' + str(start_min_seed) + ' because seed central z is outside the z range in this field. ')
                    no_gal_min_res[min_seed_index] = [start_min_seed, fit_funct_no_gal_dens(start_min_seed)]
                seed_end = time.time()
                print ('Single seed no-galaxy  fit took ' + str(seed_end - seed_start) + 's')


            if self.gal_dens_weighting > 0.0:
                min_res = [0 for seed_num in range(int(len(no_gal_min_res) * frac_seeds_to_search_for_gals))]
                no_gal_min_values = [single_seed_res[1] for single_seed_res in no_gal_min_res]
                indeces_sorted_by_min_val, no_gal_sorted_min_values = cant.safeSortOneListByAnother(no_gal_min_values, [list(range(len(no_gal_min_values))), no_gal_min_values])
                seeds_to_reexplore =  [no_gal_min_res[good_fit_index][0] for good_fit_index in indeces_sorted_by_min_val[0:int(len(no_gal_min_values) * frac_seeds_to_search_for_gals)]]
                print ('We will reexplore the following seeds for this field: ' + str(seeds_to_reexplore))
                for i in range(len(seeds_to_reexplore)):
                    seed_start = time.time()
                    new_seed_to_explore = seeds_to_reexplore[i]
                    print ('Re checking seed ' + str(new_seed_to_explore) + ', now including galaxy density. ')
                    new_min_res = optimize.minimize(fit_funct_with_gal_dens, new_seed_to_explore, method = method, tol = 10.0 ** -4.0, bounds = self.bounds)
                    min_res[i] = [new_min_res['x'], new_min_res['fun']]
                    print ('Found minimizing params with including galaxy density: ' + str(new_min_res['x']) + ' with chi_sqr_likelihood, gal_n_dens_likelihood, and overall likelihoods: ' + str([fit_funct_no_gal_dens(new_min_res['x']) * (1.0 - self.gal_dens_weighting), fit_funct_with_gal_dens(new_min_res['x']) - fit_funct_no_gal_dens(new_min_res['x']) * (1.0 - self.gal_dens_weighting) , fit_funct_with_gal_dens(new_min_res['x'])]))
                    seed_end = time.time()
                    print ('Single seed with-galaxy fit took ' + str(seed_end - seed_start) + 's')
            else:
                min_res = no_gal_min_res
            #print ('no_gal_min_res = ' + str(no_gal_min_res))
            print ('min_res = ' + str(min_res))

            best_fit_min_res = min_res[np.argmin( [min_res[j][1] for j in range(len(min_res)) ]) ]
            best_fit_params = best_fit_min_res[0]
            best_fit_val = best_fit_min_res[1]
            r_chi_square = self.computeRChiSqr(all_zs, all_RAs, all_Decs, all_muResids, all_muErrs, field_center, self.n_fit_params + 1, best_fit_params, print_params = 1)
            null_chi_square = self.computeRChiSqr(all_zs, all_RAs, all_Decs, all_muResids, all_muErrs, field_center, 0 + 1, self.null_params)


            min_res_by_field[field_num] = min_res
            r_chi_square_by_field[field_num] = r_chi_square
            null_chi_square_by_field[field_num] = null_chi_square
            fit_funct_no_gal_by_field[field_num] = fit_funct_no_gal_dens(best_fit_params)
            n_sig_gal_by_field[field_num] = list(self.computeGalDensProb(best_fit_params[0], field_center[0] + best_fit_params[2], field_center[1] + best_fit_params[3], best_fit_params[1])[0:2])
            prob_of_gal_by_field[field_num] = (best_fit_val - fit_funct_no_gal_by_field[field_num] * (1.0 - self.gal_dens_weighting)) / self.gal_dens_weighting
            overall_prob_of_gal_by_field[field_num] = best_fit_val


        if add_to_master_min_dict:
            for field_num in min_res_by_field.keys():
                self.min_res[field_num] = min_res_by_field[field_num]
        if update_field_fits_as_best_fit:
            #print ('self.min_res = ' + str(self.min_res))
            print ('self.fit_params_by_field = ' + str(self.fit_params_by_field))
            for key in self.min_res.keys():
                print('key = ' + str(key))
                print ('self.min_res[key] = ' + str(self.min_res[key]))
                best_fit_min_res = self.min_res[key][np.argmin( [self.min_res[key][j][1] for j in range(len(self.min_res[key])) ]) ]
                best_fit_params = best_fit_min_res[0]
                best_fit_val = best_fit_min_res[1]
                self.fit_params_by_field[key] = list(best_fit_params)
                self.prob_fits_by_field[key] = {'sne_chi_sqr':r_chi_square_by_field[key], 'sne_null_chi_sqr':null_chi_square_by_field[key], 'sne_chi_sqr_prob':fit_funct_no_gal_by_field[key], 'n_sig_gal_dens':n_sig_gal_by_field[key], 'n_sig_gal_prob':prob_of_gal_by_field[key], 'sne_overall_prob':overall_prob_of_gal_by_field[key]}
                print ('Set self.fit_params_by_field[' + str(key) + '] to ' + str(self.fit_params_by_field[key]))
        if update_null_chi_by_field:
            self.null_chi_square_by_field = null_chi_square_by_field
        master_end = time.time()
        print ('self.prob_fits_by_field = ' + str(self.prob_fits_by_field))
        print ('Minimizing all given fields took ' + str(master_end - master_start) + 's')

        return min_res_by_field

    def assignBestFitCurves(self):
        #Sasha's by eye best fit params:
        self.fit_params_by_field[0] = [0.248, -0.959, -3.52, 0.95]
        self.fit_params_by_field[1] = [0.266, 1.064, -0.42, -2.99]
        self.fit_params_by_field[2] = [0.420, 1.132, -1.50, -1.13]
        self.fit_params_by_field[3] = [0.349, 0.079, 0.03, -0.41]
        self.fit_params_by_field[4] = [0.438, -0.023, 0.504, 0.52]
        self.fit_params_by_field[5] = [0.296, 0.850, 2.41, 0.98]
        self.fit_params_by_field[6] = [0.371, 0.143, 0.28, 0.07]
        self.fit_params_by_field[7] = [0.298, 0.676, 2.18, 0.08]
        self.fit_params_by_field[8] = [0.137, 0.010, -0.61, -0.14]
        self.fit_params_by_field[9] = [0.298, 0.234, -0.61, -0.70]
        return 1

    #Saves the best fit information for each field to a file
    def saveFitInformation(self, save_file_name, save_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/SNIsotropyProject/SNeFieldFits/'):

        header = ['Field', 'z_c', 'log10(M/(Msun10^-17))', 'RAOff(deg)', 'DecOff(deg)', 'WeightGal', 'rChiSqr', 'rChiSqr0', 'P(rChiSqr)', 'deltaGal', 'sigDeltaGal', 'P(deltaGal/sigDeltaGal)', 'FullProb']
        rows_to_save = []
        for field_key in  self.prob_fits_by_field.keys():
            row_part1 = [field_key]+ self.fit_params_by_field[field_key]+ [self.gal_dens_weighting]
            row_part2 = [self.prob_fits_by_field[field_key]['sne_chi_sqr'], self.prob_fits_by_field[field_key]['sne_null_chi_sqr'], self.prob_fits_by_field[field_key]['sne_chi_sqr_prob']]
            row_part3 = self.prob_fits_by_field[field_key]['n_sig_gal_dens']
            row_part4 = [self.prob_fits_by_field[field_key]['n_sig_gal_prob'], self.prob_fits_by_field[field_key]['sne_overall_prob']]
            print ('[row_part1, row_part2, row_part3, row_part4] = ' + str([row_part1, row_part2, row_part3, row_part4]))
            row_to_save =  row_part1 + row_part2 + row_part3 + row_part4
            rows_to_save = rows_to_save + [row_to_save]

        cant.saveListToFile(rows_to_save, save_file_name, save_dir = save_dir, sep = ', ', append = False, header = header)

        return 1


    def __init__(self, data_set,
                      randomize_each_survey = 0, randomize_each_field = 0, n_z_bins = 10, zHD = 0, OmM = 0.3, OmL = 0.7, Om0 = 1.0, OmR = 0.0, H0 = 70.0,
                      archive_healpix_sides = 32, annulus_inner_comoving_rad_scaling = 1.0, annulus_outer_comoving_rad_scaling = 2.0, gal_dens_weighting = 0.5,
                      full_sdss_gal_data_file = 'SDSS_PSSuperField3_SDSSGals_Allpz.csv', preloaded_sdss_gals = None, start_spherical_comoving_rad_frac = 0.03,
                       interp_z_params = [0.0, 100.0, 1001], z_range = [-0.1, 3.0], pull_extinctions = 0, surveys_to_include = ['all'], surveys_to_ignore = [],
                       n_mcmc_steps = 1000, archive_to_use = 'ps1md', resid_fitting_funct = 'redshift_from_vel', min_v_scale_power = -5, fixed_comoving_scale_radius_power = 1.0,
                        ):

        self.data_set = data_set
        self.n_z_bins = n_z_bins
        self.pull_extinctions = pull_extinctions
        self.zHD = zHD
        self.OmM = OmM
        self.OmL = OmL
        self.Om0 = Om0
        self.OmR = OmR
        self.H0 = H0
        self.interp_z_params = interp_z_params
        self.astro_arch = AstronomicalParameterArchive()
        self.initialize_r_of_z_interp(interp_z_params)
        self.z_range = z_range
        self.n_mcmc_steps = n_mcmc_steps
        self.surveys_to_ignore = surveys_to_ignore
        self.surveys_to_include = surveys_to_include
        #self.min_v_scale_power = min_v_scale_power #10^(%this number) = smallest velocity, as a fraction of c, that we consider # if number goes below this value, we flip it around and look at negative velocities
        self.fixed_comoving_scale_radius_power = fixed_comoving_scale_radius_power
        self.archive_healpix_sides =  archive_healpix_sides
        self.annulus_inner_comoving_rad_scaling = annulus_inner_comoving_rad_scaling
        self.annulus_outer_comoving_rad_scaling = annulus_outer_comoving_rad_scaling
        self.full_sdss_gal_data_file = full_sdss_gal_data_file
        self.preloaded_sdss_gals = preloaded_sdss_gals
        self.start_spherical_comoving_rad_frac = start_spherical_comoving_rad_frac
        self.gal_dens_weighting = gal_dens_weighting


        self.initializeSN(self.z_range, surveys_to_include, surveys_to_ignore)
        if randomize_each_survey:
            self.randomizeSNBySurvey()

        self.archive_to_use = archive_to_use
        self.initializeSNByField(archive_to_use, randomize = randomize_each_field)

        if resid_fitting_funct == 'redshift_from_grav':
            #For a redshift due to a gravitational well
            self.fit_funct, self.muDiff_of_z_funct, self.n_fit_params, self.bounds, self.mcmc_start_params, self.mcmc_step_sizes  = self.getFittingFunctG()
        if resid_fitting_funct == 'redshift_from_vel':
            #for a redshift due to a velocity profile (v given in km/s)
            self.fit_funct, self.muDiff_of_z_funct, self.n_fit_params, self.bounds, self.mcmc_start_params, self.mcmc_step_sizes = self.getFittingFunctVelNFWOurs()
        self.fit_params_by_field = {field_key:[] for field_key in self.fields.keys()}
        self.prob_fits_by_field = {field_key:{'sne_chi_sqr':0.0, 'sne_null_chi_sqr':0.0, 'sne_chi_sqr_prob':0.0, 'n_sig_gal_dens':[0.0, 0.0], 'n_sig_gal_prob':0.0, 'sne_overall_prob':0.0} for field_key in self.fields.keys()}
        self.assignBestFitCurves()
        self.min_res = { }

        self.null_params = [0.5, 0.0, 0.1, 0.1]


def loadPickledPlotter(file_to_load):
    loaded_fielder = pickle.load(open(file_to_load, 'rb'))

    return loaded_fielder


if __name__ == "__main__":
    line_args = sys.argv[1:]
    fitter_id = line_args[0]
    field_ids = line_args[1]
    print ('field_ids = ' + str(field_ids))
    field_ids[1:-1].split(',')
    field_ids = [int(id) for id in field_ids[1:-1].split(',')]
    print('field_ids = ' + str(field_ids))
    do_randomization_by_field = 0
    do_randomization_by_survey = 0
    gal_dens_weighting = 0.5
    save_plot = 1
    fit_params_file =  ('RandField' if do_randomization_by_field else 'RandSurvey' if do_randomization_by_survey else 'True') + 'PS1MDFieldFitter_field' + '_'.join([str(elem) for elem in field_ids])  + '_GalWeight0p' + str(int(10 * gal_dens_weighting)) + '_N' + str(fitter_id) + '.csv'
    plot_file_name = ('RandField' if do_randomization_by_field else 'RandSurvey' if do_randomization_by_survey else 'True') + 'PS1MDFieldFitter_field' + '_'.join([str(elem) for elem in field_ids]) + '_GalWeight0p' + str(int(10 * gal_dens_weighting)) + '_N' + str(fitter_id) + '.pdf'
    z_range = [0.0, 1.0]
    sdssdir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/SNIsotropyProject/SDSSGalaxies/'
    fulldata_fastRead = cant.readInColumnsToList('SDSS_fullCoverage_SDSSGals_pzAll.csv', sdssdir, delimiter = ',', n_ignore = 2, all_np_readable = 1)
    print('fitter_id = ' + str(fitter_id))
    print('fit_params_file =  ' + str(fit_params_file))
    field_plotter = PanStarsFieldManager(1, full_sdss_gal_data_file = 'SDSS_fullCoverage_SDSSGals_pzAll.csv', preloaded_sdss_gals = fulldata_fastRead, randomize_each_field = do_randomization_by_field, randomize_each_survey = do_randomization_by_survey, gal_dens_weighting = gal_dens_weighting, z_range = z_range)# , surveys_to_include = ['PS1MD' ,'SDSS', 'SNLS'])
    #field_fitter =  PanStarsFieldManager(1, randomize_each_field = do_randomization_by_field, randomize_each_survey = do_randomization_by_survey, surveys_to_include = ['PS1MD' ,'SDSS', 'SNLS'])


    fit_res = field_plotter.computeByHandFits(field_ids)

    #fit_res = field_plotter.computeByHandFits([1])

    #fit_res = field_plotter.computeByHandFits([0,1,2])
    if save_plot:
        field_plotter.makePlotOfPS1MDFields(show = 0, save = 1, save_name = plot_file_name,  plot_fit = 1,)

    #field_plotter.preloaded_sdss_gals = None

    print ('field_plotter.prob_fits_by_field = ' + str(field_plotter.prob_fits_by_field))
    #pickle.dump(field_plotter, open(fit_file, "wb"))
    field_plotter.saveFitInformation(fit_params_file , save_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/SNIsotropyProject/SNeFieldFits/')

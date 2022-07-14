import numpy as np
import cantrips as can
import math
import numpy as np
import time
from PanStars1Archive import PanStars1Archive
from binData import binData
import matplotlib.pyplot as plt
from DirectoryArchive import DirectoryArchive
#from FitStorer import FitStorer
from SDSSArchive import SDSSArchive
from AstronomicalParameterArchive import AstronomicalParameterArchive
from matplotlib.patches import Rectangle
import loadSN as lsn
import GeneralMCMC as mcmc
from scipy import special
from scipy import integrate
from scipy import interpolate
from scipy import optimize
import randomSortData as rsd
import pickle
import sys
import os
import emcee
import corner
import CosmicDensityProfilesForSNePerturbationsClass as cdp

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
loaded_data_all_without_gal = [can.readInColumnsToList(data_dir + rand_file, delimiter = ',') for rand_file in rand_files_all_without_gal]
loaded_data_all_with_gal = [[can.readInColumnsToList(data_dir + rand_file, delimiter = ',') for rand_file in rand_files] for rand_files in rand_files_all_with_gal]
good_rand_indeces_with_gal = [[np.all([ param_bounds[j][0] <= float(single_fit[fit_param_indeces[j]][fields[i]+1]) <= param_bounds[j][1] for j in range(len(param_bounds)) ]) for single_fit in loaded_data_all_with_gal[i]] for i in range(len(fields))]
loaded_rand_r_chi_sqrs_with_gal = [[float(single_fit[chi_sqr_index][fields[i]+1]) for single_fit in loaded_data_all_with_gal[i]] for i in range(len(fields))]
loaded_rand_overall_probs_with_gal = [[float(single_fit[overall_prob_index][fields[i]+1]) for single_fit in loaded_data_all_with_gal[i]] for i in range(len(fields))]
good_rand_indeces_without_gal = [[np.all([ param_bounds[j][0] <= float(single_fit[fit_param_indeces[j]][fields[i]+1]) <= param_bounds[j][1] for j in range(len(param_bounds)) ]) for single_fit in loaded_data_all_without_gal] for i in range(len(fields))]
loaded_rand_r_chi_sqrs_without_gal = [[float(single_fit[chi_sqr_index][fields[i]+1]) for single_fit in loaded_data_all_without_gal] for i in range(len(fields))]
loaded_rand_overall_probs_without_gal = [[float(single_fit[overall_prob_index][fields[i]+1]) for single_fit in loaded_data_all_without_gal] for i in range(len(fields))]
true_file_without_gal = 'TruePS1MDFieldFitter_field0_1_2_3_4_5_6_7_8_9' + '_GalWeight0p0_N' + str(0) + '.csv'
true_file_with_gal = 'TruePS1MDFieldFitter_field0_1_2_3_4_5_6_7_8_9' + '_GalWeight0p5_N' + str(0) + '.csv'
true_data_without_gal = can.readInColumnsToList(data_dir + true_file_without_gal, delimiter = ',')
true_data_with_gal = can.readInColumnsToList(data_dir + true_file_with_gal, delimiter = ',')
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
    ax.text((hist_res[1][-1] - hist_res[1][0]) * 0.02 + hist_res[1][0], max(hist_res[0]) * 0.9, 'f' + str(field) + ': ' + xlabel + r', $_{T}$ = ' + str(can.round_to_n(true_res, 3)), c = 'r', horizontalalignment='left', verticalalignment='center')
    ax.set_xlabel(xlabel )
    ax.set_ylabel(r'$\Delta \mu$ (mags)')

plt.tight_layout()
plt.show()
"""


class PanStarsFieldManager:

    def getComovingCrossSectionOfAngularScaleAtRedshift(self, z, angular_scale_in_deg):
        comoving_cross_section = self.r_of_z_interp(z) * np.deg2rad(angular_scale_in_deg)
        return comoving_cross_section

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
            intermed_term = (10.0 ** (deltaMus / 5.0) - 1.0) * np.array([integrate.quad(lambda z_int: 1/self.H_of_z(z_int), 0, ref_z)[0] for ref_z in ref_zs]) * (-self.H_of_z(ref_zs) / (1.0 + ref_zs))
            print ('intermed_term = ' + str(intermed_term))
            zs_from_mus = intermed_term / (1.0- intermed_term)

        return zs_from_mus

    def getLOSVsFromZs(self, zs):
        one_p_zs = 1.0 + np.array(zs)
        vs = (one_p_zs ** 2.0 - 1.0) / (1.0 + one_p_zs ** 2.0)
        return vs

    def getErroneousRedshiftPlot(self, z_extra_space, z_obs_space,
                                 z_obs_space_errs = None,
                                 zHD = 1, z_include_range = [0.0, np.inf], surveys_of_interest = ['all'], surveys_to_excise = [''],
                                 ):
        if z_obs_space_errs == None:
            z_obs_space_errs = [0.0 for z in z_obs_space]
        r_of_z_interp = self.r_of_z_interp
        c = self.astro_arch.getc() # in km/s
        z_eff = lambda zf, zg: (zf - zg)/(zg + 1)

        #Note H0 is given in km/s/Mpc.  c is in km/s.  So c/H0 gives this ratio in Mpc.
        if r_of_z_interp is None:
            dLInMpc = lambda zf, zg: c / H0 * (1 + zf) * integrate.quad(lambda z_int: 1/self.H_of_z(z_int), 0, z_eff(zf, zg) )[0]
        else:
            #dLInMpc = lambda zf, zg: c / H0 * (1 + zf) * r_of_z_interp((zf - zg)/(zg + 1))
            dLInMpc = lambda zf, zg: (1 + zf) * r_of_z_interp(z_eff(zf, zg))
        mu =  lambda zf, zg: (5 * np.log10(dLInMpc(zf, zg)) + 25 if zf > zg else -np.inf)
        #vals = [mu(z, zg_of_z_funct(z)) - mu(z, 0.0) for z in z_space]
        vals = [mu(z_obs_space[z_index], z_extra_space[z_index]) - mu(z_obs_space[z_index], 0.0) for z_index in range(len(z_obs_space))]
        dLs_extraEffect, dLs_noEffect = [ np.array([dLInMpc(z_obs_space[z_index], z_extra_space[z_index]) for z_index in range(len(z_obs_space))]), np.array([dLInMpc(z_obs_space[z_index], 0) for z_index in range(len(z_obs_space))]) ]

        #if 0.0 in dLs_extraEffect:
        #    print ('Zero proper motion corrected dL detected in array: ' + str(dLs_extraEffect.tolist()))
        #    print ( '[z_eff(z_obs_space[z_index], z_extra_space[z_index]) for z_index in range(len(z_obs_space))] = ' + str([z_eff(z_obs_space[z_index], z_extra_space[z_index]) for z_index in range(len(z_obs_space))]))
        #if (0.0 in dLs_noEffect):
        #    print ('Zero straight up dL detected in array: ' + str(dLs_extraEffect))

        mu_diffs_with_extra_a = np.where(dLs_extraEffect > 0.0, 1 / (dLs_extraEffect), 0.0)
        mu_diffs_with_extra_b =  np.array( [self.H_of_z(z_eff(z_obs_space[z_index], z_extra_space[z_index])) * (1.0 + z_extra_space[z_index]) for z_index in range(len(z_obs_space))] )
        mu_diffs_with_extra = mu_diffs_with_extra_a * mu_diffs_with_extra_b
        #if 0.0 in dLs_extraEffect:
        #    print ('mu_diff_with_extra = ' + str(mu_diffs_with_extra))
        mu_diff_no_extra = np.array([1 / (dLs_noEffect * self.H_of_z(z_obs_space[z_index])) for z_index in range(len(z_obs_space))])
        partial_mu_diffs_in_z = 5.0 / np.log(10) * (mu_diffs_with_extra - mu_diff_no_extra)
        mu_diff_errs_from_z_errs = np.abs(np.array(partial_mu_diffs_in_z) * np.array(z_obs_space_errs))
        vals = np.where(np.isnan(vals), np.inf, vals)
        mu_diff_errs_from_z_errs = np.where(np.isnan(vals), 0, mu_diff_errs_from_z_errs)

        return vals, mu_diff_errs_from_z_errs

    def plotPantheonSNOnSky(self, data_set,
                            projection = 'aitoff', save = 1, show = 0, surveys_to_display = 'all', z_range = [-0.1, 3.0],
                            pull_extinctions = 0, figsize = [10.0, 8.0], zHD = 1):

        #dir_archive = DirectoryArchive()
        plot_dir = self.dir_archive.getPlotDirectory()

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
        all_sns = lsn.loadSN(1, ['all'], pull_extinctions = pull_extinctions, zHD = zHD)
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

    def redshift_of_v_funct(self, beta, coses_of_los_angle, max_beta = 0.999):
        beta = np.where(np.abs(beta) < 1.0, beta, max_beta * np.sign(beta))
        redshift = ((1 + beta * coses_of_los_angle)/np.sqrt(1-abs(beta) ** 2.0) - 1)
        return redshift

    def redshiftDueToGravPotential (self, zs, central_z, delta, R):
        speedol = self.astro_arch.getc()
        well_depth = -delta * self.OmM * self.H0 ** 2.0 * R ** 2.0 / 2.0
        central_r = self.r_of_z_interp(central_z)
        phi_of_delta = [well_depth * (2 - (abs(self.r_of_z_interp(single_z) - central_r) / R) ** 2.0 if abs(self.r_of_z_interp(single_z) - central_r) < R else R / abs(self.r_of_z_interp(single_z) - central_r)) for single_z in zs]
        extra_z_of_zMeas_funct_G = (phi_of_delta(single_z, float(self.r_of_z_interp(central_z)), delta, R) - phi_of_delta(0.0, float(self.r_of_z_interp(central_z)), delta, R)) / (speedol) ** 2.0
        fitted_G, fitted_G_errs = np.array(self.getErroneousRedshiftPlot(lambda single_z: extra_z_of_zMeas_funct_G(single_z, central_z, delta, R), zs, zHD = self.zHD, astro_arch = self.astro_arch, r_of_z_interp = self.r_of_z_interp ))
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
        d_z_D_growth_over_D_term2 = -(1.0 + single_z) / self.H_of_z(single_z) ** 3.0 / integrate.quad(lambda z_int: (1 + z_int) / self.H_of_z(z_int) ** 3.0, single_z, np.inf)[0]
        d_z_D_growth_over_D = d_z_D_growth_over_D_term1 + d_z_D_growth_over_D_term2
        return d_z_D_growth_over_D

    #The geometry of how our line of sight intersects the halo
    def r_incidence_leg (self, single_r, central_r, r_incidence):
        leg =  (np.sqrt(central_r ** 2.0 - r_incidence ** 2.0 ) - single_r)
        if np.isnan(leg).any():
            print ('[single_r, central_r, r_incidence, leg] = ' + str([single_r, central_r, r_incidence, leg]  ))
        return leg

    #This is what I was doing previously which I think incorrectly assumes right angle
    #def r_well_from_geometry(self, single_r, central_r, r_incidence ):
    #    r_well =  np.sqrt(r_incidence ** 2.0 + self.r_incidence_leg(single_r, central_r, r_incidence) ** 2.0)
    #    return r_well

    #This is the correct formula, I think.
    def r_well_from_geometry(self, single_rs, central_r, angular_offsets ):
        #delta_dists = single_rs - central_r
        r_wells =  np.sqrt(single_rs ** 2.0 + central_r ** 2.0 - 2.0 * single_rs * central_r * np.cos(angular_offsets))
        return r_wells

    #This seems wrong?
    def sin_of_infall_angle_relative_to_los(self, single_r, central_r, r_incidence):
        sin_r = self.r_incidence_leg(single_r, central_r, r_incidence) / self.r_well_from_geometry(single_r, central_r, r_incidence)
        return sin_r

    #I believe this is the correct formula, though there is the outstanding question of handling this in an expanding Universe.
    def cos_of_infall_angle_relative_to_los(self, single_rs, central_r, angular_offsets):
        r_seps = self.r_well_from_geometry(single_rs, central_r, angular_offsets)
        cos_rs = (single_rs ** 2.0 + r_seps ** 2.0 - central_r ** 2.0) / (2.0 * single_rs * r_seps)
        return cos_rs

    #This is the formula I was using previously...
    #def cos_of_infall_angle_relative_to_los(self, single_r, central_r, r_incidence):
    #    cos_r = 1 - np.cos()
    #    return cos_r

    """
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
        M_enc = halo_density * (np.log(radius_in_scale_radii + 1.0) - radius_in_scale_radii / (1.0 + radius_in_scale_radii))
        return M_enc

    def point_mass_enclosed(self, halo_mass):
        #r_overdensity = self.r_overdensity(mass_overdensity, critical_density, overdensity)
        #M_enc = self.nfw_unitless_over_int(concentration * r / r_overdensity) / self.nfw_unitless_over_int(concentration) * mass_overdensity
        return halo_mass
    """

    """
    def r_overdensity(self, M_overdensity, critical_density, overdensity):
        r = ((M_overdensity * 3.0) / (4.0 * np.pi * overdensity * critical_density )) ** (1.0 / 3.0)
        return r
    """

    # returns the critical mass with units of T M_{sun} Mpc ^ -3 - remember T = tera = 10 ^ 12
    def crit_density(self, z):
        G = self.astro_arch.getGravitationalConstant()
        crit_density_unitless_part = 3.0 * self.H_of_z(z) ** 2.0 / (8.0 * np.pi)
        pc_to_m = self.astro_arch.getParsecToM()
        Mpc_to_m = pc_to_m * 10.0 ** 6.0
        km_to_m = 1000.0
        Msun_to_kg = self.astro_arch.getSolarMass()
        kg_to_TMsun = 1.0 / (Msun_to_kg * 10.0 ** 12)
        crit_density_unitfull_part = self.H0 ** 2.0 / G * (km_to_m) ** 2.0 * (Mpc_to_m) * (kg_to_TMsun)
        #(1 / Mpc_to_km) ** 2.0 / (10.0 ** 17.0 * Msun_to_kg) * (pc_to_m * 10.0 ** 6.0) ** 3.0
        crit_density = crit_density_unitless_part * crit_density_unitfull_part
        #print ('[crit_density_unitless_part, crit_density_unitfull_part] = ' + str([crit_density_unitless_part, crit_density_unitfull_part]))
        return crit_density



    #Mass should be in M_sun X 10 ^(12) ; returns r ** 2.0 * int(delta r*2 4 pi dr) in units of Mpc
    def delta_portion_funct(self, single_rs, central_z, angular_offsets, model_params):
        central_r = self.r_of_z_interp(central_z) #units of Mpc
        critical_density = self.crit_density(central_z) #units of GMsun / Mpc ^ 3
        OmM_of_z = self.OmM * (1.0 + central_z) ** 3.0
        r_wells = self.r_well_from_geometry(single_rs, central_r, angular_offsets)

        #Those models that are based on unitless overdensities rather than massess
        #    need the overdensity to be passed as an extra parameter.
        if self.resid_profile_funct in ['exp_void']:
            extra_params = [critical_density, OmM_of_z]
        else:
            extra_params = None
        M_encs = self.radial_mass_funct(r_wells, model_params, extra_params = extra_params)
        #print ('angular_offsets = ' + str(angular_offsets))
        #print ('r_wells = ' + str(r_wells))
        #print ('model_params = ' + str(model_params))
        M_encs_scaled = M_encs / (4.0 * np.pi * r_wells ** 2.0 * OmM_of_z * critical_density)
        deltas_over_rsqr_int_over_rsqr = M_encs_scaled

        return deltas_over_rsqr_int_over_rsqr


    '''
    #Mass should be in M_sun X 10 ^(12) ; returns r ** 2.0 * int(delta r*2 4 pi dr) in units of Mpc
    # delta_portion_funct_ours(single_z, central_z, halo_density, r_scale, b_impact, overdensity)
    def delta_portion_funct_NFW_ours(self, single_zs, central_z, r_cutoff, r_scale, angular_offsets, overdensity_param = None  ):

        if overdensity_param == None:
            overdensity_param = self.overdensity_param
        single_rs = self.r_of_z_interp(single_zs)
        central_r = self.r_of_z_interp(central_z) #units of Mpc
        critical_density = self.crit_density(central_z) #units of GMsun / Mpc ^ 3
        OmM_of_z = self.OmM * (1.0 + central_z) ** 3.0
        #r_overdensity = self.r_overdensity(M_overdensity, critical_density, overdensity)
        #rs = r_overdensity / concentration
        r_wells = self.r_well_from_geometry(single_rs, central_r, angular_offsets)
        #concentration = self.computeConcentrationForNFW(OmM_of_z, critical_density, halo_mass, overdensity_param)
        halo_mass = r_cutoff ** 3.0 * (4.0 * np.pi * overdensity_param * OmM_of_z * critical_density) / (3.0)

        M_encs = halo_mass * self.nfw_unitless_over_int(r_wells / r_scale )
        print ('[r_cutoff, r_scale, halo_mass] = ' + str([r_cutoff, r_scale, halo_mass] ))
        M_encs_scaled = M_encs / (4.0 * np.pi * r_wells ** 2.0 * self.nfw_unitless_over_int(r_cutoff / r_scale) * OmM_of_z * critical_density)
        print ('M_encs_scaled = ' + str(M_encs_scaled))

        deltas_over_rsqr_int_over_rsqr = - M_encs_scaled

        return deltas_over_rsqr_int_over_rsqr

    #Unitless scaling should be unitless scalar (can be positive or negative) ; returns 4 pi / 3 delta_0 r^3  e^(-(r/r_0) ^ 2) * rho_c in units of Mpc^3
    def delta_portion_funct_exponential(self, single_zs, central_z, unitless_overdensity, r_scale, angular_offsets):

        single_rs = self.r_of_z_interp(single_zs)
        central_r = self.r_of_z_interp(central_z) #units of Mpc
        critical_density = self.crit_density(central_z) #units of GMsun / Mpc ^ 3
        OmM_of_z = self.OmM * (1.0 + central_z) ** 3.0
        #r_overdensity = self.r_overdensity(M_overdensity, critical_density, overdensity)
        #rs = r_overdensity / concentration
        r_wells = self.r_well_from_geometry(single_rs, central_r, angular_offsets)
        #concentration = self.computeConcentrationForNFW(OmM_of_z, critical_density, halo_mass, overdensity_param)
        M_encs = np.pi * 4.0 / 3.0 * critical_density * OmM_of_z * unitless_overdensity * r_wells ** 3.0 * np.exp(-(r_wells / r_scale) ** 2.0)
        M_encs_scaled = M_encs / (4.0 * np.pi * r_wells ** 2.0 * OmM_of_z * critical_density)
        #print ('M_encs_scaled = ' + str(M_encs_scaled))

        deltas_over_rsqr_int_over_rsqr = - M_encs_scaled

        return -deltas_over_rsqr_int_over_rsqr

    #Mass should be in M_sun X 10 ^(12) ; returns r ** 2.0 * int(delta r*2 4 pi dr) in units of Mpc
    # delta_portion_funct_ours(single_z, central_z, halo_density, r_scale, b_impact, overdensity)
    def delta_portion_funct_point_mass(self, single_zs, central_z, angular_offsets, halo_mass ):
        single_rs = self.r_of_z_interp(single_zs) #Units in Mpc
        central_r = self.r_of_z_interp(central_z)
        critical_density = self.crit_density(central_z) #units of GMsun / Mpc ^ 3
        OmM_of_z = self.OmM * (1.0 + central_z) ** 3.0
        r_wells = self.r_well_from_geometry(single_rs, central_r, angular_offsets)
        M_enc = halo_mass
        M_enc_scaled = M_enc / (OmM_of_z * critical_density)
        delta_over_rsqr_int_over_rsqr =  1.0 / (4.0 * np.pi * r_wells ** 2.0) * M_enc_scaled

        return delta_over_rsqr_int_over_rsqr

    #Mass should be in M_sun X 10 ^(12) ; returns r ** 2.0 * int(delta r*2 4 pi dr) in units of Mpc
    # delta_portion_funct_ours(single_z, central_z, halo_density, r_scale, b_impact, overdensity)
    def delta_portion_funct_uniform_mass(self, single_zs, central_z, angular_offsets, halo_mass, cutoff_radius ):
        single_rs = self.r_of_z_interp(single_zs) #Units in Mpc
        central_r = self.r_of_z_interp(central_z)
        critical_density = self.crit_density(central_z) #units of GMsun / Mpc ^ 3
        OmM_of_z = self.OmM * (1.0 + central_z) ** 3.0
        r_wells = self.r_well_from_geometry(single_rs, central_r, angular_offsets)
        #M_enc = self.uniform_mass_enclosed_ours(halo_mass, radius_in_scale_radii)
        #M_enc_scaled = M_enc / (OmM_of_z * critical_density)
        delta_over_rsqr_int_over_rsqr = np.where(r_wells < cutoff_radius, halo_mass * (r_wells / cutoff_radius) ** 3.0, halo_mass ) / (OmM_of_z * critical_density * 4.0 * np.pi * r_wells ** 2.0)
        return delta_over_rsqr_int_over_rsqr

    def getHaloMassFromVariedParam(self, halo_mass_power, linear_scaling_term = 10 ** 3.0):
        """
        We vary a stand-in parameter for the mass.  Here is where
           we convert back to the physical mass, in 10^12 Mpc.
        """
        halo_mass = np.sign(halo_mass_power) * (10.0 ** abs(halo_mass_power) - 1.0) #Allow the sampling to go from -whatever to + whatever, passing through 0 when halo_mass_power = 0
        halo_mass = halo_mass * linear_scaling_term
        return halo_mass

    def getRadiusFromVariedParam(self, radius_power, linear_scaling_term = 10 ** 0.0):
        """
        We vary a stand-in parameter for the radii parameters.  Here is where
           we convert back to the physical mass, in 10^12 Mpc.
        """
        radius = 10.0 ** radius_power
        radius = radius * linear_scaling_term
        return radius

    def getLOSVelocityFieldAlongLOSExponentialVoid(self, zs, RAs, Decs, field_center, central_z, unitless_overdensity, comoving_scale_radius_power, central_sky_RA_offset, central_sky_Dec_offset, ):
        """
        For a set of sky coordinates in 3 space (redshift, angles on sky), determine
           their line of sight (los) peculiar velocities induced by an exponential DM
           halo over/under density (negative masses are underdensities).
        We rotate the spherical coordinates to be centered on some central
            coordinates, central_sky_RA_offset, central_sky_Dec_offset.
        The exponential halo mass and scale radius give a characteristic halo mass
            density, which is the paramter usually used to specify halo density.
            delta(r) = -delta_0 (1 - 2 r ^ 2 / (3 r_0 ^ 3)) e ^ (-r^2/r_0^2).  Note
            the TOTAL mass in this void is 0.
        """
        #comoving_scale_radius_power = self.fixed_comoving_scale_radius_power
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
        #total_mass = 10.0 ** total_mass_power #mass, in M_{sun} * 10 ^ 12
        r_scale = 10.0 ** comoving_scale_radius_power #scale radius in Mpc
        angular_offsets = np.array([can.measureAngularSeparationOnSky([RAs[i], Decs[i]], [central_RA, central_Dec], return_radian = 1) for i in range(len(zs))])
        coses_of_los_angle = self.cos_of_infall_angle_relative_to_los(single_rs, central_r, angular_offsets)
        delta_portions = self.delta_portion_funct_exponential(zs, central_z, unitless_overdensity, r_scale, angular_offsets) #in units of Mpc^2
        speedol = self.astro_arch.getc()
        los_velocity_field = self.H0 / speedol * np.array(Hs_of_z) * np.array(d_zs_D_over_D) * delta_portions
        return los_velocity_field, coses_of_los_angle

    #parameters are the mass in one scale radius, the scale radius, the central z, and angles describing the impact parameter
    #def getLOSVelocityFieldAlongLOSOurs(self, zs, RAs, Decs, field_center, central_z, halo_scaling_power, comoving_scale_radius_power, central_sky_phi, central_sky_theta, ):
    def getLOSVelocityFieldAlongLOSNFWOurs(self, zs, RAs, Decs, field_center, central_z, comoving_scale_radius_power, cutoff_radius, overdensity_param, central_sky_RA_offset, central_sky_Dec_offset, ):
        """
        For a set of sky coordinates in 3 space (redshift, angles on sky), determine
           their line of sight (los) peculiar velocities induced by an NFW DM halo
           over/under density (negative masses are underdensities).
        We rotate the spherical coordinates to be centered on some central
            coordinates, central_sky_RA_offset, central_sky_Dec_offset.
        The NFW halo mass is specified for some concentration parameter, Delta, and
            a specified cutoff radius (the distance at which the mass density falls)
            to 0.
        """
        #comoving_scale_radius_power = self.fixed_comoving_scale_radius_power
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
        #halo_mass = self.getHaloMassFromVariedParam(halo_mass_power)
        #vel_at_scale_radius = np.sign(vel_at_scale_radius_power) * (10.0 ** (-abs(vel_at_scale_radius_power) - 1.0)
        #rolled_v_scale_power = np.where(vel_at_scale_radius_power > self.min_v_scale_power, vel_at_scale_radius_power, -vel_at_scale_radius_power + 2.0 * self.min_v_scale_power) #Rolls over to correspond to negative velocities at the minimium allowed power
        #vel_at_scale_radius = np.where(vel_at_scale_radius_power > self.min_v_scale_power, 1.0, -1.0) * (10.0 ** rolled_v_scale_power  - 10.0 ** self.min_v_scale_power) #self.max_v_scale_power is the given power at which the velocity at z = 0 is c
        #speedol = self.astro_arch.getc()
        #halo_scaling  = 4.0 * np.pi * vel_at_scale_radius / (r_scale * self.H0 / speedol * self.nfw_mass_enclosed_ours(1.0, 1.0))
        #halo_scaling = mass_scaling / (self.nfw_mass_enclosed_ours(1.0, 1.0) * r_scale ** 3.0)
        #print ('[vel_at_scale_radius_power, vel_at_scale_radius, halo_scaling, r_scale] = ' + str([vel_at_scale_radius_power, vel_at_scale_radius, halo_scaling, r_scale]))
        angular_offsets = np.array([can.measureAngularSeparationOnSky([RAs[i], Decs[i]], [central_RA, central_Dec], return_radian = 1) for i in range(len(zs))])

        #r_offsets = np.sqrt(single_rs ** 2.0 + central_r ** 2.0 - 2.0 * single_rs * central_r * np.cos(angular_offsets * deg_to_rad))
        #r_incidences = np.sin(angular_offsets / 2.0) * central_r * 2.0
        #print ('angular_offsets.tolist() = ' + str(angular_offsets.tolist()))
        #print ('r_incidences.tolist() = ' + str(r_incidences.tolist()))
        coses_of_los_angle = self.cos_of_infall_angle_relative_to_los(single_rs, central_r, angular_offsets)
        #print ('[zs, central_z, halo_scaling, r_scale, r_incidences] = ' + str([zs, central_z, halo_scaling, r_scale, r_incidences] ))
        delta_portions = self.delta_portion_funct_NFW_ours(zs, central_z, cutoff_radius, r_scale, angular_offsets, overdensity_param = overdensity_param) #in units of Mpc
        #print ('delta_portions = ' + str(delta_portions))
        speedol = self.astro_arch.getc()
        #self.H0 / speedol gives H0 / c in units of 1/MPC
        #print ('[np.isnan(Hs_of_z).any(), np.isnan(d_zs_D_over_D).any(), np.isnan(sins_of_los_angle).any(), np.isnan(delta_portions).any()] = ' + str([np.isnan(Hs_of_z).any(), np.isnan(d_zs_D_over_D).any(), np.isnan(sins_of_los_angle).any(), np.isnan(delta_portions).any()]))
        #print ('[np.isinf(Hs_of_z).any(), np.isinf(d_zs_D_over_D).any(), np.isinf(sins_of_los_angle).any(), np.isinf(delta_portions).any()] = ' + str([np.isinf(Hs_of_z).any(), np.isinf(d_zs_D_over_D).any(), np.isinf(sins_of_los_angle).any(), np.isinf(delta_portions).any()]))

        #los_velocity_field = self.H0 / speedol * np.array(Hs_of_z) * np.array(d_zs_D_over_D) * np.array(sins_of_los_angle) * delta_portions
        print ('r_scale = ' + str(r_scale))
        los_velocity_field = self.H0 / speedol * np.array(Hs_of_z) * np.array(d_zs_D_over_D) * delta_portions
        print ('los_velocity_field = ' + str(los_velocity_field))
        #print ('los_velocity_field = ' + str(los_velocity_field))
        #print ('coses_of_los_angle = ' + str(coses_of_los_angle))
        return los_velocity_field, coses_of_los_angle

    #parameters are the point mass in 10^XX M_sun, the central z, and angles describing the impact parameter
    def getLOSVelocityFieldAlongLOSUniformMass(self, zs, RAs, Decs, field_center, central_z, halo_mass_power, cutoff_radius_power, central_sky_RA_offset, central_sky_Dec_offset, ):
        """
        For a set of sky coordinates in 3 space (redshift, angles on sky), determine
           their line of sight (los) peculiar velocities induced by a spherical,
           uniform over/under density (negative masses are underdensities).
        We rotate the spherical coordinates to be centered on some central
            coordinates, central_sky_RA_offset, central_sky_Dec_offset.
        """
        #print ('[zs, RAs, Decs, field_center, central_z, halo_mass_power, cutoff_radius_power, central_sky_RA_offset, central_sky_Dec_offset] = ' + str([zs, RAs, Decs, field_center, central_z, halo_mass_power, cutoff_radius_power, central_sky_RA_offset, central_sky_Dec_offset]))
        deg_to_rad = self.astro_arch.getDegToRad()
        field_RA, field_Dec = field_center
        #central_RA, central_Dec = [field_RA + np.sin(central_sky_phi * deg_to_rad) * central_sky_theta * np.cos(field_Dec * deg_to_rad), field_Dec + np.cos(central_sky_phi * deg_to_rad) * central_sky_theta]
        central_RA, central_Dec = [field_RA + central_sky_RA_offset, field_Dec + central_sky_Dec_offset]
        #print ('[central_RA, central_Dec] = ' + str([central_RA, central_Dec]))
        Hs_of_z = [self.H_of_z(single_z) for single_z in zs]
        d_zs_D_over_D = [self.d_z_D_growth_over_D(single_z) for single_z in zs]
        #Relationship between comoving radial distance and redshift is determined by cosmology.
        single_rs = np.array([self.r_of_z_interp(single_z) for single_z in zs])
        central_r = self.r_of_z_interp(central_z)
        critical_density = self.crit_density(central_z)
        halo_mass = self.getHaloMassFromVariedParam(halo_mass_power)
        cutoff_radius = self.getRadiusFromVariedParam(cutoff_radius_power)
        #print ('[halo_mass_power, halo_mass, cutoff_radius_power, cutoff_radius] = ' + str([halo_mass_power, halo_mass, cutoff_radius_power, cutoff_radius]))
        angular_offsets = np.array([can.measureAngularSeparationOnSky([RAs[i], Decs[i]], [central_RA, central_Dec], return_radian = 1) for i in range(len(zs))])
        #r_incidences = np.sin(angular_offsets / 2.0) * central_r * 2.0
        #r_seps = self.r_well_from_geometry(single_rs, central_r, angular_offsets)
        coses_of_los_angle = self.cos_of_infall_angle_relative_to_los(single_rs, central_r, angular_offsets)
        delta_portions = self.delta_portion_funct_uniform_mass(zs, central_z, angular_offsets, halo_mass, cutoff_radius) #in units of Mpc

        speedol = self.astro_arch.getc()
        #self.H0 / speedol gives H0 / c in units of 1/MPC
        los_velocity_field = self.H0 / speedol * np.array(Hs_of_z) * np.array(d_zs_D_over_D) * delta_portions
        return los_velocity_field, coses_of_los_angle

    '''

    def getVelocityField(self, zs, RAs, Decs, field_center, central_sky_RA_offset, central_sky_Dec_offset, profile_params, ):
         """
         For a set of sky coordinates in 3 space (redshift, angles on sky), determine
            their line of sight (los) peculiar velocities induced by a the spherical
            overdensity model of the data.
         We rotate the spherical coordinates to be centered on some central
             coordinates, central_sky_RA_offset, central_sky_Dec_offset.
         """
         #print ('field_center = ' + str(field_center))
         #print ('profile_params = ' + str(profile_params))
         deg_to_rad = self.astro_arch.getDegToRad()
         central_z, field_RA, field_Dec = field_center
         #central_RA, central_Dec = [field_RA + np.sin(central_sky_phi * deg_to_rad) * central_sky_theta * np.cos(field_Dec * deg_to_rad), field_Dec + np.cos(central_sky_phi * deg_to_rad) * central_sky_theta]
         central_RA, central_Dec = [field_RA + central_sky_RA_offset, field_Dec + central_sky_Dec_offset]
         #print ('[central_RA, central_Dec] = ' + str([central_RA, central_Dec] ))
         #print ('RAs = ' + str(RAs))
         #print ('Decs = ' + str(Decs))
         Hs_of_z = [self.H_of_z(single_z) for single_z in zs]
         d_zs_D_over_D = [self.d_z_D_growth_over_D(single_z) for single_z in zs]
         #Relationship between comoving radial distance and redshift is determined by cosmology.
         single_rs = np.array([self.r_of_z_interp(single_z) for single_z in zs])
         central_r = self.r_of_z_interp(central_z)
         critical_density = self.crit_density(central_z)
         angular_offsets = np.array([can.measureAngularSeparationOnSky([RAs[i], Decs[i]], [central_RA, central_Dec], return_radian = 1) for i in range(len(zs))])
         #print ('angular_offsets = ' + str(angular_offsets))
         coses_of_los_angle = self.cos_of_infall_angle_relative_to_los(single_rs, central_r, angular_offsets)

         delta_portions = self.delta_portion_funct(single_rs, central_z, angular_offsets, profile_params) #in units of Mpc
         speedol = self.astro_arch.getc()
         #print ('[Hs_of_z, d_zs_D_over_D, delta_portions] = ' + str([Hs_of_z, d_zs_D_over_D, delta_portions]))
         los_velocity_field = self.H0 / speedol * np.array(Hs_of_z) * np.array(d_zs_D_over_D) * delta_portions
         return los_velocity_field, coses_of_los_angle

    '''
    #parameters are the point mass in 10^XX M_sun, the central z, and angles describing the impact parameter
    def getLOSVelocityFieldAlongLOSPointMass(self, zs, RAs, Decs, field_center, halo_mass_power, central_sky_RA_offset, central_sky_Dec_offset, ):
        """
        For a set of sky coordinates in 3 space (redshift, angles on sky), determine
           their line of sight (los) peculiar velocities induced by a point mass
           over/under density (negative masses are underdensities).
        We rotate the spherical coordinates to be centered on some central
            coordinates, central_sky_RA_offset, central_sky_Dec_offset.
        """
        central_z, field_RA, field_Dec = field_center
        comoving_scale_radius_power = self.fixed_comoving_scale_radius_power
        deg_to_rad = self.astro_arch.getDegToRad()
        #central_RA, central_Dec = [field_RA + np.sin(central_sky_phi * deg_to_rad) * central_sky_theta * np.cos(field_Dec * deg_to_rad), field_Dec + np.cos(central_sky_phi * deg_to_rad) * central_sky_theta]
        central_RA, central_Dec = [field_RA + central_sky_RA_offset, field_Dec + central_sky_Dec_offset]
        Hs_of_z = [self.H_of_z(single_z) for single_z in zs]
        d_zs_D_over_D = [self.d_z_D_growth_over_D(single_z) for single_z in zs]
        #Relationship between comoving radial distance and redshift is determined by cosmology.
        single_rs = np.array([self.r_of_z_interp(single_z) for single_z in zs])
        central_r = self.r_of_z_interp(central_z)
        critical_density = self.crit_density(central_z)
        #r_scale = 10.0 ** comoving_scale_radius_power #scale radius in Mpc
        halo_mass = self.getHaloMassFromVariedParam(halo_mass_power)
        angular_offsets = np.array([can.measureAngularSeparationOnSky([RAs[i], Decs[i]], [central_RA, central_Dec], return_radian = 1) for i in range(len(zs))])
        #r_incidences = np.sin(angular_offsets / 2.0) * central_r * 2.0
        #r_seps = self.r_well_from_geometry(single_rs, central_r, angular_offsets)
        coses_of_los_angle = self.cos_of_infall_angle_relative_to_los(single_rs, central_r, angular_offsets)

        delta_portions = self.delta_portion_funct(zs, central_z, angular_offsets, halo_mass) #in units of Mpc

        speedol = self.astro_arch.getc()
        #self.H0 / speedol gives H0 / c in units of 1/MPC
        los_velocity_field = self.H0 / speedol * np.array(Hs_of_z) * np.array(d_zs_D_over_D) * delta_portions
        return los_velocity_field, coses_of_los_angle

    '''

    def getMuDiffOfVelocityField(self, calc_zs, calc_RAs, calc_Decs, ref_zs, ref_RAs, ref_Decs, ref_muDiffs, ref_muErrs, field_center_on_sky, vel_field_params, print_steps = 0):
        all_zs = np.array(calc_zs).tolist() + np.array(ref_zs).tolist()
        all_RAs = np.array(calc_RAs).tolist() + np.array(ref_RAs).tolist()
        all_Decs = np.array(calc_Decs).tolist() + np.array(ref_Decs).tolist()

        central_z, RA_offset, Dec_offset = vel_field_params[0:3]
        RA_offset_to_center = (180.0 - RA_offset)
        #print ('field_center_on_sky = ' + str(field_center_on_sky[0]))
        #print ('all_RAs = ' + str(all_RAs))
        all_RAs = self.centerRAs(all_RAs,  RA_offset_to_center )
        #print ('all_RAs = ' + str(all_RAs))
        central_RA, central_Dec = [ self.centerRAs([field_center_on_sky[0]],  RA_offset_to_center )[0], field_center_on_sky[1] ]
        #print ('central_RA = ' + str(central_RA))

        halo_params = vel_field_params[3:]
        field_center = [central_z] + [central_RA, central_Dec]
        scaled_velocities, coses_of_los_angles = self.getVelocityField(all_zs, all_RAs, all_Decs, field_center, RA_offset, Dec_offset, halo_params)
        vel_redshifts = self.redshift_of_v_funct(scaled_velocities, coses_of_los_angles)
        #print ('vel_redshifts = ' + str(vel_redshifts))
        #for z_index in range(len(all_zs)):
        #    if vel_redshifts[z_index] > all_zs[z_index]:
        #        print ('[all_zs[z_index], all_RAs[z_index], all_Decs[z_index], scaled_velocities[z_index], coses_of_los_angles[z_index], vel_redshifts[z_index]] = ' + str([all_zs[z_index], all_RAs[z_index], all_Decs[z_index], scaled_velocities[z_index], coses_of_los_angles[z_index], vel_redshifts[z_index]]  ))
        #print('[vel_redshifts, all_zs] = ' + str([vel_redshifts, all_zs]))
        erroneous_mus, erroneous_mu_errs = self.getErroneousRedshiftPlot(vel_redshifts, all_zs,)
        ref_calced_muDiffs = erroneous_mus[-len(ref_zs):]
        weights = np.array(ref_muErrs) ** (-2.0)
        weighted_diffs = (np.array(ref_calced_muDiffs) - np.array(ref_muDiffs) ) * weights
        #bad_diffs = np.where(np.isinf(weighted_diffs), 0.0, 1.0)
        #weighted_mean_diff = np.sum(np.where(bad_diffs, weighted_diffs, 0.0)) / np.sum(weights * bad_diffs)
        weighted_mean_diff = np.sum(weighted_diffs) / np.sum(weights)
        #Do we want to subtract away the mean residual?
        weighted_mean_diff = 0.0
        if 0:
            print ('erroneous_mus.tolist() = ' + str(erroneous_mus.tolist()))
            print ('weighted_mean_diff = ' + str(weighted_mean_diff))
            print ('That is for vel_field_params = '  + str(vel_field_params))
            print ('weights.tolist() = ' + str(weights.tolist()))
            print ('ref_calced_muDiffs.tolist() = ' + str(ref_calced_muDiffs.tolist()))
            print ('vel_redshifts.tolist() = ' + str(vel_redshifts.tolist()) )
            print ('scaled_velocities.tolist() = ' + str(scaled_velocities.tolist()))
            print ('coses_of_los_angles.tolist() = ' + str(coses_of_los_angles.tolist()))
            print ('all_zs = ' + str(all_zs))
        if np.isinf(weighted_mean_diff):
            corrected_muDiffs = [weighted_mean_diff for elem in np.array(erroneous_mus[0:-len(ref_zs)])]
        else:
            corrected_muDiffs = np.array(erroneous_mus[0:-len(ref_zs)]) - weighted_mean_diff
        return corrected_muDiffs

    '''
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

    def getFittingFunctVelPointMass(self, seed_z_step = 0.01, n_seeds_around_sky = 4):
        n_vel_fit_params = 4
        #Note: for mass scaling (second) parameter, the number is logarithmic
        #BOUNDS ARE: offsets from centered z, logarithm of mass, relative to 10 ^ 15 solar masses,
        #            offsets from object position on sky, in degrees
        vel_bounds = [[0.322, 0.333], (0.0, 2.0), (-5.0, 5.0), (-5.0, 5.0)]
        mcmc_step_sizes = [0.005, 0.1, 0.05, 0.05]
        sub_funct = self.getLOSVelocityFieldAlongLOSPointMass
        muDiff_of_z_funct = self.getMuDiffOfVelocityField
        MCMC_chain_starts  = []
        seed_zs = np.arange(vel_bounds[0][0], vel_bounds[0][1], seed_z_step)
        for seed_z in seed_zs:
            MCMC_chain_starts = MCMC_chain_starts + [[seed_z, 0.0, 0.0, 0.0]] #Put one seed at the center of the field
            for j in range(n_seeds_around_sky):
                MCMC_chain_starts = MCMC_chain_starts  + [[seed_z, 0.0 * np.log10(2.0), np.cos(2.0 * np.pi / n_seeds_around_sky * (0.5 + j)) * 1.0, np.sin(2.0 * np.pi / n_seeds_around_sky * (0.5 + j)) * 1.0]]
        return [sub_funct, muDiff_of_z_funct, n_vel_fit_params, vel_bounds, MCMC_chain_starts, mcmc_step_sizes]

    #def getFittingFunctVelNFWOurs(self, n_seeds_in_z = 50, n_seeds_around_sky = 4):
    def getFittingFunctVelNFWOurs(self, seed_z_step = 0.01, n_seeds_around_sky = 4):
        n_vel_fit_params = 4
        #Note: for scaling (second) and scale radius (third) params, the numbers are logarithmic
        #vel_bounds = [(0.0, 1.0), (0.0, 10.0), (0.5, 4.0), (0.0, 360.0), (0.01, 10.0)]
        #vel_bounds = [(0.0, 1.0), (0.0, 5.0), (0.0, 360.0), (0.01, 10.0)]
        #BOUNDS ARE: offsets from centered z, logarithm of mass, relative to 10 ^ 15 solar masses,
        #            scale radius, relative to XXXMpc, offsets from object position on sky, in degrees
        #vel_bounds = [(-0.10, 0.10), (-3.0, 3.0), (-20.0, 20.0), (-20.0, 20.0)]
        #vel_bounds = [[0.01, self.z_range[1]], (0.0, 2.0), (-5.0, 5.0), (-5.0, 5.0)]
        vel_bounds = [[0.322, 0.333], (0.0, 2.0), (-5.0, 5.0), (-5.0, 5.0)]
        #fit_central_zs = np.linspace(0.0, 1.0, 21)
        #fit_halo_scaling_powers = np.linspace(0.1, 4.0, 9) #In T M_{sun}
        #fit_r_scale_powers = np.linspace(0.1, 3, 10)
        #fit_impact_params = np.linspace(0.1, 3.0, 6)
        #fit_overdensities = [200]
        #mcmc_step_sizes = [0.005, 0.02, 0.015, 2.0, 0.05]
        mcmc_step_sizes = [0.005, 0.1, 0.05, 0.05]
        loaded_los_velocity_funct = lambda zs, RAs, Decs, field_center, central_z, comoving_scale_radius_power, central_sky_RA_offset, central_sky_Dec_offset: self.getLOSVelocityFieldAlongLOSNFWOurs(zs, RAs, Decs, field_center, central_z, mass_scale_power, comoving_scale_radius_power, central_sky_RA_offset, central_sky_Dec_offset)
        sub_funct = loaded_los_velocity_funct
        muDiff_of_z_funct = self.getMuDiffOfVelocityField
        MCMC_chain_starts  = []
        seed_zs = np.arange(vel_bounds[0][0], vel_bounds[0][1], seed_z_step)
        for seed_z in seed_zs:
            MCMC_chain_starts = MCMC_chain_starts + [[seed_z, 0.0, 0.0, 0.0]] #Put one seed at the center of the field
            for j in range(n_seeds_around_sky):
                MCMC_chain_starts = MCMC_chain_starts  + [[seed_z, 0.0 * np.log10(2.0), np.cos(2.0 * np.pi / n_seeds_around_sky * (0.5 + j)) * 1.0, np.sin(2.0 * np.pi / n_seeds_around_sky * (0.5 + j)) * 1.0]]
        #MCMC_chain_starts = MCMC_chain_starts  + [[vel_bounds[0][1] / (n_seeds_in_z - 1) * i + vel_bounds[0][0], np.log10(2.0), 0.05, 0.05] for i in range(n_seeds_in_z)]
        #print ('MCMC_chain_starts = ' + str(MCMC_chain_starts))
        return [sub_funct, muDiff_of_z_funct, n_vel_fit_params, vel_bounds, MCMC_chain_starts, mcmc_step_sizes]

    def getFittingFunctVelExponentialVoid(self, seed_z_step = 0.01, n_seeds_around_sky = 4):
        n_vel_fit_params = 4
        #Note: for density (first) and scale radius (second) params, the numbers are logarithmic
        #vel_bounds = [(0.0, 1.0), (0.0, 10.0), (0.5, 4.0), (0.0, 360.0), (0.01, 10.0)]
        #vel_bounds = [(0.0, 1.0), (0.0, 5.0), (0.0, 360.0), (0.01, 10.0)]
        #BOUNDS ARE: offsets from centered z, logarithm of mass, relative to 10 ^ 15 solar masses,
        #            scale radius, relative to XXXMpc, offsets from object position on sky, in degrees
        #vel_bounds = [(-0.10, 0.10), (-3.0, 3.0), (-20.0, 20.0), (-20.0, 20.0)]
        #vel_bounds = [[0.01, self.z_range[1]], (0.0, 2.0), (-5.0, 5.0), (-5.0, 5.0)]
        vel_bounds = [[-10, 0.333], (0.0, 2.0), (-5.0, 5.0), (-5.0, 5.0)]
        #fit_central_zs = np.linspace(0.0, 1.0, 21)
        #fit_halo_scaling_powers = np.linspace(0.1, 4.0, 9) #In T M_{sun}
        #fit_r_scale_powers = np.linspace(0.1, 3, 10)
        #fit_impact_params = np.linspace(0.1, 3.0, 6)
        #fit_overdensities = [200]
        #mcmc_step_sizes = [0.005, 0.02, 0.015, 2.0, 0.05]
        mcmc_step_sizes = [0.005, 0.1, 0.05, 0.05]
        loaded_los_velocity_funct = lambda zs, RAs, Decs, field_center, central_z, unitless_overdensity, comoving_scale_radius_power, central_sky_RA_offset, central_sky_Dec_offset: self.getLOSVelocityFieldAlongLOSExponentialVoid(zs, RAs, Decs, field_center, central_z, unitless_overdensity, comoving_scale_radius_power, central_sky_RA_offset, central_sky_Dec_offset)
        sub_funct = loaded_los_velocity_funct
        muDiff_of_z_funct = self.getMuDiffOfVelocityField
        MCMC_chain_starts  = []
        seed_zs = np.arange(vel_bounds[0][0], vel_bounds[0][1], seed_z_step)
        for seed_z in seed_zs:
            MCMC_chain_starts = MCMC_chain_starts + [[seed_z, 0.0, 0.0, 0.0]] #Put one seed at the center of the field
            for j in range(n_seeds_around_sky):
                MCMC_chain_starts = MCMC_chain_starts  + [[seed_z, 0.0 * np.log10(2.0), np.cos(2.0 * np.pi / n_seeds_around_sky * (0.5 + j)) * 1.0, np.sin(2.0 * np.pi / n_seeds_around_sky * (0.5 + j)) * 1.0]]
        #MCMC_chain_starts = MCMC_chain_starts  + [[vel_bounds[0][1] / (n_seeds_in_z - 1) * i + vel_bounds[0][0], np.log10(2.0), 0.05, 0.05] for i in range(n_seeds_in_z)]
        #print ('MCMC_chain_starts = ' + str(MCMC_chain_starts))
        return [sub_funct, muDiff_of_z_funct, n_vel_fit_params, vel_bounds, MCMC_chain_starts, mcmc_step_sizes]



    def getFittingFunctVelUniformMass(self, seed_z_step = 0.01, n_seeds_around_sky = 4):
        #Which are the correct numbers of degrees of freedom?  The number of fit parameters on the sky, or the number run at each position on the grid?
        n_vel_fit_params = 5
        #Note: for scaling (second) and scale radius (third) params, the numbers are logarithmic
        #vel_bounds = [(0.0, 1.0), (0.0, 10.0), (0.5, 4.0), (0.0, 360.0), (0.01, 10.0)]
        #vel_bounds = [(0.0, 1.0), (0.0, 5.0), (0.0, 360.0), (0.01, 10.0)]
        #BOUNDS ARE: offsets from centered z, logarithm of mass, relative to 10 ^ 15 solar masses,
        #            scale radius, relative to XXXMpc, offsets from object position on sky, in degrees
        #vel_bounds = [(-0.10, 0.10), (-3.0, 3.0), (-20.0, 20.0), (-20.0, 20.0)]
        #vel_bounds = [[0.01, self.z_range[1]], (0.0, 2.0), (-5.0, 5.0), (-5.0, 5.0)]
        vel_bounds = [[0.322, 0.333], (0.0, 2.0), (-5.0, 5.0), (-5.0, 5.0)]
        #fit_central_zs = np.linspace(0.0, 1.0, 21)
        #fit_halo_scaling_powers = np.linspace(0.1, 4.0, 9) #In T M_{sun}
        #fit_r_scale_powers = np.linspace(0.1, 3, 10)
        #fit_impact_params = np.linspace(0.1, 3.0, 6)
        #fit_overdensities = [200]
        #mcmc_step_sizes = [0.005, 0.02, 0.015, 2.0, 0.05]
        mcmc_step_sizes = [0.005, 0.1, 0.01, 0.05, 0.05]
        sub_funct = self.getLOSVelocityFieldAlongLOSUniformMass
        muDiff_of_z_funct = self.getMuDiffOfVelocityField
        MCMC_chain_starts  = []
        seed_zs = np.arange(vel_bounds[0][0], vel_bounds[0][1], seed_z_step)
        for seed_z in seed_zs:
            MCMC_chain_starts = MCMC_chain_starts + [[seed_z, 0.0, 2.0, 0.0, 0.0]] #Put one seed at the center of the field
            for j in range(n_seeds_around_sky):
                MCMC_chain_starts = MCMC_chain_starts  + [[seed_z, 0.0 * np.log10(2.0), 2.0, np.cos(2.0 * np.pi / n_seeds_around_sky * (0.5 + j)) * 1.0, np.sin(2.0 * np.pi / n_seeds_around_sky * (0.5 + j)) * 1.0]]
        #MCMC_chain_starts = MCMC_chain_starts  + [[vel_bounds[0][1] / (n_seeds_in_z - 1) * i + vel_bounds[0][0], np.log10(2.0), 0.05, 0.05] for i in range(n_seeds_in_z)]
        #print ('MCMC_chain_starts = ' + str(MCMC_chain_starts))
        return [sub_funct, muDiff_of_z_funct, n_vel_fit_params, vel_bounds, MCMC_chain_starts, mcmc_step_sizes]
    '''

    def computeRChiSqr(self, all_zs, all_RAs, all_Decs, all_resids, all_errs, field_center, n_model_params, fit_params, print_params = 0):
        #print ('fit_params = ' + str(fit_params))
        fitted_resids = self.muDiff_of_z_funct(all_zs, all_RAs, all_Decs, all_zs, all_RAs, all_Decs, all_resids, all_errs, field_center, fit_params)
        #weights = np.array(all_errs) ** (-2.0)
        #weighted_mean_diff = np.sum(np.array(fitted_resids - all_resids)  * weights) / np.sum(weights)
        #weighted_mean_diff = 0.0
        total_dof = len(all_zs) - n_model_params
        chi_sqr = np.sum(((np.array(fitted_resids)) - np.array(all_resids)) ** 2.0 / (np.array(all_errs) ** 2.0))
        #print ('chi_sqr = ' + str(chi_sqr))
        r_chi_sqr = chi_sqr / total_dof
        #if print_params: print ('np.array(fit_params).tolist() + [r_chi_sqr] = ' + str(np.array(fit_params).tolist()+ [r_chi_sqr]))
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

    def adjustRAForSkyPlot(self, RA):
        x = np.remainder(RA+360-0,360) # shift RA values
        ind = x>180
        x[ind] -=360    # scale conversion to [-180, 180]
        x=-x    # reverse the scale: East to the left
        return x

    def makeSkyPlotOfSNe(self, sky_plot ):
        sn_by_survey = self.sn_by_survey
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
            RA = self.adjustRAForSkyPlot(RA)
            Dec = [(sn['Dec']) for sn in sn_in_survey]
            plots_for_legend[i] = sky_plot.scatter(np.radians(RA), np.radians(Dec) , color = [sn['color']  for sn in sn_in_survey], marker = 'o', s = 2.0 )
        return plots_for_legend


    def makePlotOfPS1MDFields(self, fields_to_plot = 'all', plot_fit = 0,
                              plot_delta_vs = 0,
                              save = 1, show = 0, fit_information = {'funct':'none'}, show_fit_label = 1, n_plot_zs = 1001, z_range_to_plot = [-0.1, 3.0],
                              mu_plot_lims = [-0.6, 0.6], vel_plot_lims = [-0.07, 0.07], z_plot_lims = [0.0, 1.0], archive_to_use = 'PS1MD' , figsize = [16.0, 8.0], save_name = 'SN_residuals_v_z_PS1MD_fields_.pdf',
                              xlabel = r'$z$', ylabel = r'$\Delta \mu$ (mags)',
                              n_z_bins = 10, super_plot_title = r'$\Delta \mu$ of Pan Stars 1 medium-deep fields', plot_param_round_to = 3, plot_chi_sqr_round_to = 3, min_fig_side = 3 ):

        master_start = time.time()
        #dir_archive = DirectoryArchive()
        plot_dir = self.dir_archive.getPlotDirectory()
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
        plot_indeces = [[fig_side - 1, i] for i in range(fig_side)] + [[0, i] for i in range(fig_side)] + can.flattenListOfLists( [[[i, 0], [i, fig_side - 1]] for i in range(1, fig_side-1)] )
        for plot_index in plot_indeces:
            fig.add_subplot(gs[plot_index[0], plot_index[1]])
        #Add the central sky plot
        sky_plot = fig.add_subplot(gs[1:fig_side-1, 1:fig_side-1], projection="aitoff")
        plots_for_legend = self.makeSkyPlotOfSNe(sky_plot)
        """
        sky_plot.grid(True)
        sky_plot.set_xlabel('R.A.')
        sky_plot.set_ylabel('Decl.')
        tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
        tick_labels = np.remainder(tick_labels+360+0,360)
        #fig = plt.figure(figsize=(10, 5))
        #ax = fig.add_subplot(111, projection=projection, axisbg ='LightCyan')
        sky_plot.set_xticklabels(tick_labels)     # we add the scale on the x axis
        for i in range(len(plots_for_legend)):
            sn_in_survey = sn_by_survey[i]
            RA = np.array([(sn['RA']) for sn in sn_in_survey])
            Dec = [(sn['Dec']) for sn in sn_in_survey]
            x = np.remainder(RA+360-0,360) # shift RA values
            ind = x>180
            x[ind] -=360    # scale conversion to [-180, 180]
            x=-x    # reverse the scale: East to the left
            plots_for_legend[i] = sky_plot.scatter(np.radians(x), np.radians(Dec) , color = [sn['color']  for sn in sn_in_survey], marker = 'o', s = 2.0 )
        """
        #plots_for_legend = [plt.scatter([], [], c = sn_in_survey[0]['color'], marker = 'o') for  sn_in_survey in sn_by_survey if len(sn_in_survey) > 0]
        legend_plot_index = plot_indeces[-1]
        ax = fig.add_subplot(gs[legend_plot_index[0], legend_plot_index[1]])
        ax.legend(plots_for_legend, self.included_surveys, ncol = 3, fontsize = 8)

        for j in range(len(field_nums)):
            field_num = field_nums[j]
        #for field_num in range(1):
            field_str = 'f' + str(field_num) + ' '
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

            all_zs = can.flattenListOfLists(zs_by_survey)
            all_RAs = can.flattenListOfLists(ras_by_survey)
            all_Decs = can.flattenListOfLists(decs_by_survey)
            all_resids  = can.flattenListOfLists(resids_by_survey)
            all_errs  = can.flattenListOfLists(mu_errs_by_survey)

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
                chi_sqr = self.computeRChiSqr(all_zs, all_RAs, all_Decs, all_resids, all_errs, field_center, self.n_fit_params, fit_params)
                null_chi_sqr = self.computeRChiSqr(all_zs, all_RAs, all_Decs, all_resids, all_errs, field_center, self.n_null_fit_params, null_params)
                label_str = label_str + 'params ' + str( [can.round_to_n(fit_param, plot_param_round_to) for fit_param in fit_params]) + r' with $\chi^2_{\nu}$ ' + str(can.round_to_n(chi_sqr, plot_chi_sqr_round_to)) + r'; $\chi^2_{\nu,0}$ ' + str(can.round_to_n(null_chi_sqr, plot_chi_sqr_round_to))
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
        r = (speedol) / self.H0 * integrate.quad( lambda z_int: 1.0 / self.H_of_z(z_int), 0, z_meas )[0]
        return r

    def dL_of_z(self, z_meas):
        comoving_dist = self.r_of_z(z_meas)
        dL = (1.0 + z_meas) * comoving_dist
        return dL

    def mu_of_z(self, z_meas):
        dL = self.dL_of_z(z_meas)
        mu = 25.0 + 5.0 * np.log10(dL)
        return mu

    def initialize_z_of_r_interp(self, interp_z_params):
        self.interp_zs = np.linspace(interp_z_params[0], interp_z_params[1], interp_z_params[2]).tolist() + [1000.0]
        self.z_of_r_interp = interpolate.interp1d([0.0] + [self.r_of_z(z) for z in self.interp_zs], [-1000.0] + self.interp_zs, fill_value = (0.0, self.r_of_z(1000)), bounds_error = False)
        return 1

    def initialize_r_of_z_interp(self, interp_z_params):
        self.interp_zs = np.linspace(interp_z_params[0], interp_z_params[1], interp_z_params[2]).tolist() + [1000.0]
        self.r_of_z_interp =interpolate.interp1d([-1000.0] + self.interp_zs, [0.0] + [self.r_of_z(z) for z in self.interp_zs], fill_value = (0.0, self.r_of_z(1000)), bounds_error = False)
        return 1

    #Randomly pair the SNe residuals and residual uncertainties with zs,
    def randomizeSNBySurvey(self, sn_to_randomize, rand_indeces):
        all_surveys = list(set([sn['survey'] for sn in sn_to_randomize]))
        randomized_sn = []
        sn_by_survey = {survey:[] for survey in all_surveys}
        for sn in sn_to_randomize:
            sn_by_survey[sn['survey']] = sn_by_survey[sn['survey']] + [sn]
        for survey in all_surveys:
            sn_in_survey = sn_by_survey[survey]
            sn_in_survey = self.randomizeAllSN(sn_in_survey)
            randomized_sn = randomized_sn + sn_in_survey
        #self.all_sn = randomized_sn
        return randomized_sn

    def randomizeAllSN(self, sn_to_randomize, rand_indeces):
        """
        Shuffle the number of mu sigmas by which a supernova is off, relative
            to the spatial coordinates (z, RA, Dec).  The mu uncertainties
            are kept, relative to z, RA, and Dec, but the numbers of such mu
            errors by which the mus of those SN are off from cosmology are
            randomly rearranged.
        The new numbers of sigma off are then use to redetermine the mus at
            each spatial coordinate.
        """
        all_zs = [sn['mu'] for sn in sn_to_randomize]
        all_mus = [sn['mu'] for sn in sn_to_randomize]
        all_muResids = [sn['mu'] - self.mu_of_z(sn['z']) for sn in sn_to_randomize]
        all_muErrs = [sn['muErr'] for sn in sn_to_randomize]
        all_nMuSigs = (np.array(all_muResids) / np.array(all_muErrs))
        #rand_mus, rand_muDiffs, rand_muErrs = rsd.randomShuffleListOfLists([all_mus, all_muDiffs, all_muErrs])
        rand_nMuSigs = rsd.randomShuffleListOfLists([all_nMuSigs])[0]
        rand_nMuSigs = [all_nMuSigs[i] for i in rand_indeces]
        for j in range(len(sn_to_randomize)):
            sn_to_randomize[j]['muTrue'] = sn_to_randomize[j]['mu']
            sn_to_randomize[j]['mu'] = self.mu_of_z(sn_to_randomize[j]['z']) + rand_nMuSigs[j] * all_muErrs[j]
            #sn_to_randomize[j]['muDiffTrue'] = sn_to_randomize[j]['muDiff']
            #sn_to_randomize[j]['muDiff'] = rand_muDiffs[j]
            #sn_to_randomize[j]['muErrTrue'] = sn_to_randomize[j]['muErr']
            #sn_to_randomize[j]['muErr'] = rand_muErrs[j]
        return sn_to_randomize

    def generateRandomizedSNe(self, sn_to_randomize, rand_sn_data_list, rand_sn_data_dir,  randomize_by_survey, randomize_all_sn):
        n_sn = len(sn_to_randomize)
        orig_indeces = list(range(n_sn))
        orig_paired_indeces = []
        rand_indeces = []
        if randomize_all_sn:
            rand_indeces = rsd.randomShuffleListOfLists([orig_indeces])[0]
        elif randomize_by_survey:
            all_surveys = list(set([sn['survey'] for sn in sn_to_randomize]))
            sn_indeces_by_survey = {survey:[] for survey in all_surveys}
            for i in range(len(sn_to_randomize)):
                sn_indeces_by_survey[sn['survey']] = sn_by_survey[sn['survey']] + [i]
            for survey in all_surveys:
                sn_indeces_in_survey = sn_by_survey[survey]
                orig_paired_indeces = orig_paired_indeces + sn_indeces_in_survey
                rand_sn_indeces_in_survey = rsd.randomShuffleListOfLists([sn_indeces_in_survey])[0]
                rang_indeces = rand_indeces + rand_sn_indeces_in_survey
            orig_paired_indeces, rand_indeces = can.safeSortOneListByAnother(orig_paired_indeces, [orig_paired_indeces, rand_indeces])
        else:
            rand_indeces = orig_indeces
        can.saveListsToColumns([orig_indeces, rand_indeces], rand_sn_data_list, rand_sn_data_dir, header = 'orig, rand', sep = ' ')
        return [orig_indeces, rand_indeces]

    def loadRandomizedSN(self, all_sns, surveys_to_include, randomized_sn_number, randomize_by_survey, randomize_all):
        rand_sn_data_dir = self.dir_archive.getRandomizedDrawResultsDir()
        rand_sn_data_list = self.rand_sn_file_prefix + str(randomized_sn_number) + '.txt'
        if os.path.exists(rand_sn_data_dir + rand_sn_data_list):
            randomize_matching = can.readInColumnsToList(rand_sn_data_dir + rand_sn_data_list, delimiter = ' ', n_ignore = 1, convert_to_int = 1)
        else:
            randomize_matching = self.generateRandomizedSNe(all_sns, rand_sn_data_list, rand_sn_data_dir, randomize_by_survey, randomize_all_sn)
        if randomize_all:
            all_sns = self.randomizeAllSN(all_sns, randomize_matching[1])
        elif randomize_by_survey:
            all_sns = self.randomizeSNBySurvey(all_sns, randomize_matching[1])
        return all_sns


    def initializeSN(self, z_range, surveys_to_include, surveys_to_ignore, randomized_sn_number = 0, randomize_by_survey = 0, randomize_all = 0):

        #print ('[self.data_set, surveys_to_include, self.sn_data_type, self.pull_extinctions, self.zHD, self.OmM, self.OmL, self.Om0, self.OmR, self.H0] = ' + str([self.data_set, surveys_to_include, self.sn_data_type, self.pull_extinctions, self.zHD, self.OmM, self.OmL, self.Om0, self.OmR, self.H0]))
        all_sns = lsn.loadSN(self.data_set, surveys_to_include, data_type = self.sn_data_type, pull_extinctions = self.pull_extinctions, zHD = self.zHD, OmM = self.OmM, OmL = self.OmL, Om0 = self.Om0, OmR = self.OmR, H0 = self.H0, )
        print ('randomized_sn_number = ' + str(randomized_sn_number ))
        if randomized_sn_number > 0:
            print ('Loading randomize SNe... ')
            all_sns = self.loadRandomizedSN(all_sns, surveys_to_include, randomized_sn_number, randomize_by_survey, randomize_all_sn)

        if surveys_to_include[0] == 'all':
            surveys_to_include = [sn['survey'] for sn in all_sns]
            surveys_to_include = list(set(surveys_to_include))
        self.surveys_to_include = surveys_to_include
        self.surveys_to_ignore = surveys_to_ignore
        all_sns = [sn for sn in all_sns if not(sn['survey'] in self.surveys_to_ignore)]
        all_sns = [sn for sn in all_sns if (sn['z'] >= z_range[0] and sn['z'] <= z_range[1]) ]
        self.all_sns = all_sns
        self.all_zs = [sn['z'] for sn in all_sns]
        self.all_RAs = [sn['RA'] for sn in all_sns]
        self.all_Decs = [sn['Dec'] for sn in all_sns]
        self.all_mus = [sn['mu'] for sn in all_sns]
        print ('self.all_mus[0:10] = ' + str(self.all_mus[0:10]))
        self.all_mu_errs = [sn['muErr'] for sn in all_sns]
        self.all_plot_colors = [sn['color'] for sn in all_sns]
        self.all_surveys = [sn['survey'] for sn in all_sns]
        self.survey_to_color_dict = {self.all_surveys[i]:self.all_plot_colors[i] for i in range(len(self.all_surveys))}
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
        #print('self.fields = ' + str(self.fields))

        #print ('self.surveys_to_include = ' + str(self.surveys_to_include ) )
        self.sn_by_survey = [[sn for sn in self.all_sns if sn['survey'] == survey] for survey in self.surveys_to_include if len([sn for sn in self.all_sns if sn['survey'] == survey]) > 0]
        self.included_surveys = [sn_by_single_survey[0]['survey'] for sn_by_single_survey in self.sn_by_survey]
        #print ('self.included_surveys = ' + str(self.included_surveys ))
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
            all_zs = can.flattenListOfLists(all_zs_in_field_by_survey)
            all_muResids = can.flattenListOfLists(all_muResids_in_field_by_survey)
            all_muErrs = can.flattenListOfLists(all_muErrs_in_field_by_survey)
            MCMC_fit_funct = lambda params: self.computeRChiSqr(all_zs, all_muResids, all_muErrs, field_center, params)
            new_MCMC_fitter = mcmc.MCMCFitter(MCMC_fit_funct, mcmc_start_params, mcmc_step_sizes, n_mcmc_steps, bounds = self.bounds, likelihood_from_chisqr  = 1)
            null_chi_sqr = MCMC_fit_funct(self.null_params)
            MCMC_fitters_by_field [field_num] = {'fit':new_MCMC_fitter, 'null_chi_square': null_chi_sqr}
        self.MCMC_fits = MCMC_fitters_by_field

        return 1

    #This is the cumulative density function, according to Wikipedia:
    #   https://en.wikipedia.org/wiki/Chi-squared_distribution
    def computeChiSqrProb(self, all_zs, all_RAs, all_Decs, all_muResids, all_muErrs, field_center, n_fit_params, params, print_params = 0):
        dof = len(all_zs) - n_fit_params
        rchi_sqr = self.computeRChiSqr(all_zs, all_RAs, all_Decs, all_muResids, all_muErrs, field_center, n_fit_params, params, print_params = print_params)
        #print ('chi_sqr = ' + str(chi_sqr))
        chi_sqr_prob = 1.0 - special.gammainc(dof / 2.0, rchi_sqr * dof / 2.0) # chi_sqr ** ((dof - 2) / 2) * np.exp(-chi_sqr / 2.0)
        if chi_sqr_prob > 0.0:
            log_prob = np.log10(chi_sqr_prob)
        else:
            #fitted_resids = self.muDiff_of_z_funct(all_zs, all_RAs, all_Decs, all_zs, all_RAs, all_Decs, all_muResids, all_muErrs, field_center, params)
            #print ('(np.array(fitted_resids) - np.array(all_muResids)).tolist() = ' + str((np.array(fitted_resids) - np.array(all_muResids)).tolist()))
            #chi_sqr = np.sum(((np.array(fitted_resids)) - np.array(all_muResids)) ** 2.0 / (np.array(all_muErrs) ** 2.0))
            #print ('calculated chi_sqr = ' + str(chi_sqr))
            #print ('[dof, rchi_sqr, chi_sqr_prob] = ' + str([dof, rchi_sqr, chi_sqr_prob]))
            #log_prob = -np.inf
            log_prob = self.low_log_prob_val
        if print_params: print ('[dof, rchi_sqr, chi_sqr_prob, log_prob] = ' + str([dof, rchi_sqr, chi_sqr_prob, log_prob] ))
        #print ('[chi_sqr, chi_sqr_prob, log_prob] = ' + str([chi_sqr, chi_sqr_prob, log_prob]))
        #print ('params = ' + str(params) + ' => log_prob = ' + str(log_prob))
        #print ('log_prob = ' + str(log_prob))
        return rchi_sqr, log_prob

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
        n_gal_density, n_gal_density_err = self.archive.measureRelativeGalDensityAroundCoord(targ_z, obj_center[0], obj_center[1], inner_comoving_rad, annulus_inner_comoving_rad, annulus_outer_comoving_rad, verbose = 1)
        #print ('[n_gal_density, n_gal_density_err] = ' + str([n_gal_density, n_gal_density_err]))
        if halo_excess_mass > 0.0:
            normal_gal_dens_prob = 1.0 - 0.5 * (1.0 + special.erf((n_gal_density - 1.0) / (n_gal_density_err * np.sqrt(2.0))))
        else:
            normal_gal_dens_prob = 0.5 * (1.0 + special.erf((n_gal_density - 1.0) / (n_gal_density_err * np.sqrt(2.0))))

        return n_gal_density, np.log10(n_gal_density_err), np.log10(normal_gal_dens_prob)

    def overall_log_prob_funct(self, all_zs, all_RAs, all_Decs, all_muResids, all_muErrs, field_center, ext_params, print_params = 1, gal_dens_weight = 0, n_model_params = None, null_r_chi_sqr = None, minimize_chi_sqr_ratio = 0) :

        if n_model_params is None:
            n_model_params = self.n_fit_params
        sne_rchi_sqr, sne_fit = self.computeChiSqrProb(all_zs, all_RAs, all_Decs, all_muResids, all_muErrs, field_center, n_model_params, ext_params, print_params = print_params)
        sne_fit = sne_fit * (1.0 - gal_dens_weight)
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
            if sphere_center[1] < 180.0: sphere_center[1] = -360.0 - sphere_center[1]
            sphere_center[1] = -90 - sphere_center[1]

        if gal_dens_weight > 0.0:
            gal_dens_fit = self.computeGalDensProb(ext_params[0], sphere_center[0], sphere_center[1], ext_params[1])[-1] * gal_dens_weight
        else:
            gal_dens_fit = 0.0
        #print ('gal_dens_prob = ' + str(gal_dens_prob))
        overall_fit = sne_fit + gal_dens_fit
        #overall_prob = sne_probf
        if print_params:
            print ('[ext_params, sne_fit, gal_dens_fit, overall_fit] = ' + str([list(ext_params), sne_fit, gal_dens_fit, overall_fit]))
        if minimize_chi_sqr_ratio and null_r_chi_sqr != None:
            return -np.log10(sne_rchi_sqr / null_chi_sqr)
        else:
            return overall_fit

    def mcmc_prior_funct(self, params, print_params = 0):
        if np.all([params[i] > self.bounds[i][0] and params[i] < self.bounds[i][1] for i in range(len(params))]):
            prior_prob = 0.0
        else:
            prior_prob = -np.inf
        if print_params: print ('[params, prior_prob] = ' + str([params, prior_prob]))
        return prior_prob

    def centerRAs(self, RAs_to_center, RA_offset_to_center, undo_centering = 0):
        if undo_centering:
            uncentered_RAs = [(RA - RA_offset_to_center) % 360 for RA in RAs_to_center]
            return uncentered_RAs
        else:
            centered_RAs = [ RAs_to_center[j] + RA_offset_to_center for j in range(len(RAs_to_center)) ]
            centered_RAs = [RA - 360.0 if RA > 360.0 else RA + 360.0 if RA < 0.0 else RA for RA in centered_RAs ]
            return centered_RAs

    def getBoundsGivingMaxComovingRange(self, comoving_bound, central_coord, nominal_bounds):
         if not(comoving_bound == None):
             #If we are given a comoving bound, we change the bounds of the fitter to not move the fit point outside of those bounds.
             #  But the bounds need to be updated based on teh redshift - different comovign distances subtend different angular sizes at different cosmologies!
             comoving_dist_of_z = self.r_of_z_interp(central_coord[0])
             deg_sep_of_comoving_at_redshift = comoving_bound / comoving_dist_of_z / self.deg_to_rad
             local_bounds = [ (np.max([self.z_of_r_interp(comoving_dist_of_z - comoving_bound) , 0.0]), self.z_of_r_interp(comoving_dist_of_z + comoving_bound) )  ,
                                   nominal_bounds[1],
                                   (- deg_sep_of_comoving_at_redshift, deg_sep_of_comoving_at_redshift), (- deg_sep_of_comoving_at_redshift, deg_sep_of_comoving_at_redshift) ]
         else:
             local_bounds = [[central_coord[0] + nominal_bounds[0][0], central_coord[0] + nominal_bounds[0][1]], nominal_bounds[1], nominal_bounds[2], nominal_bounds[3] ]
         print ('local_bounds = ' + str(local_bounds))
         return local_bounds

    def functToMinimizeMCMC(self, params, funct, center = None, bounds = None, penalization_dist_from_params = None, penalization_scaling_from_params = None, ):
        val = funct(params)
        val = - 10.0 ** (-val)
        penalty = 0.0
        if center != None and penalization_dist_from_params != None and penalization_scaling_from_params != None:
            penalty = penalty - np.sum([0.0 if abs(params[i] - center[i]) < penalization_dist_from_params[i] else np.exp(abs(params[i] - center[0]) * penalization_scaling_from_params[i]) for i in range(len(params))] )
        if bounds != None:
            penalty = penalty - np.sum([-1 + np.exp(abs(min(0, params[i] - bounds[i][0])) * penalization_scaling_from_params[i])
                                -1 + np.exp(abs(max(0, params[i] - bounds[i][1])) * penalization_scaling_from_params[i]) for i in range(len(params))] )
        if penalty < -10:
            penalty = -np.inf
        else:
            penalty = -10.0 ** (-penalty)
        val = val + penalty
        #print ('[params] = ' + str(params) + ' => vals = ' + str(val))
        return val

    def doMCMCFitNearParams(self, fit_funct, center, bounds = None,
                            mcmc_step_size = 1e-5,
                            n_mcmc_chains = 5, n_mcmc_steps = None, seed_scaling = 10.0 ** -4.0,
                            burn_in_frac = 0.25, thinning = 5, corner_plot_labels = None,
                            penalization_dist_from_params = None, penalization_scaling_from_params = None,
                            show_mcmc = 1, mcmc_corner_save_file_name = 'MCMC_corner_plot.pdf', mcmc_chains_save_file_name = 'MCMC_chains_plot.pdf' ):
        """
        Do an MCMC minimization (using emcee) near an initial point.  Ideally, this is the
            reined minimization after a crude, grid based minimization has already been
            done.
        """
        if show_mcmc == None:
            show_mcmc = self.show_mcmc
        if n_mcmc_steps == None:
            n_mcmc_steps = self.n_mcmc_steps
        n_params = len(center)
        print ('center = ' + str(center))
        n_mcmc_chains = max(n_mcmc_chains, 2 * n_params + 1)
        loaded_funct_for_mcmc = lambda params: self.functToMinimizeMCMC(params, fit_funct, center = center, penalization_dist_from_params = penalization_dist_from_params, penalization_scaling_from_params = penalization_scaling_from_params)
        pos = center + seed_scaling * np.random.randn(n_mcmc_chains, n_params)
        nwalkers, ndim = pos.shape
        self.mcmc_sampler = emcee.EnsembleSampler(nwalkers, ndim, loaded_funct_for_mcmc, moves=[
        (emcee.moves.GaussianMove(mcmc_step_size), 1.0),
    ], )
        self.mcmc_sampler.run_mcmc(pos, n_mcmc_steps, progress = True)
        print ('burn_in_frac, n_mcmc_steps = ' + str([burn_in_frac, n_mcmc_steps ]))
        samples = self.mcmc_sampler.get_chain(discard = int(burn_in_frac * n_mcmc_steps), thin = thinning, flat = False)
        print ('np.shape(samples) = ' + str(np.shape(samples)))
        flat_samples = self.mcmc_sampler.get_chain(discard = int(burn_in_frac * n_mcmc_steps), thin = thinning, flat = True)
        param_medians = [np.percentile(flat_samples[:, i], [50])[0] for i in range(n_params)]
        print ('param_medians = ' + str(param_medians))
        fit = {'x':param_medians, 'fun':fit_funct(np.array(param_medians))}
        print ('fit = ' + str(fit))
        if show_mcmc:
            fig_corner = corner.corner( flat_samples, truths = param_medians, labels = corner_plot_labels)
            plt.savefig(self.mcmc_plot_save_dir + mcmc_corner_save_file_name)
            plt.close('all')
            f, axarr = plt.subplots(ndim,1)
            for i in range(ndim):
                ax = axarr[i]
                ax.plot(samples[:,:,i], c = 'k', alpha = 0.3)
                ax.set_xlim(0, len(samples))
                ax.set_ylabel(corner_plot_labels[i])
                ax.yaxis.set_label_coords(-0.1, 0.5)
            axarr[-1].set_xlabel('step number')
            plt.savefig(self.mcmc_plot_save_dir + mcmc_chains_save_file_name)
            plt.close('all')
        return fit

    def doStaticFitOnSNeSubset(self, central_coord, nearby_sn_indeces,
                             method = 'mcmc', init_guess = None, bounds = None, one_d_fit = 0, n_grid_samples = None, null_param_vals = [0.0, 2.0], print_params = 0,
                             show_coarse_fit = 0):
        '''
        Fit a cosmic over/under density to a subset of SNe, given with a set of indeces.
            The fit keeps the location of the overdensity fixed, changing only its
            fundumental parameters (size and mass, for example).
        '''
        if bounds == None:
            bounds = self.bounds
        if init_guess == None:
            init_guess = self.param_init_guess
        central_RA, central_Dec = central_coord[1:]
        #print ('[init_guess, bounds, n_grid_samples] = ' + str([init_guess, bounds, n_grid_samples]))
        fitting_zs, fitting_RAs, fitting_Decs, fitting_mu_resids, fitting_muErrs, fitting_surveys = [[arr[i] for i in nearby_sn_indeces] for arr in [self.all_zs, self.all_RAs, self.all_Decs, self.all_muResids, self.all_mu_errs, self.all_surveys]]
        #Are we running a minimization algorithm or an emcee?  The determines if we want the - here or not
        fit_funct_no_gal_dens = lambda params: -self.overall_log_prob_funct(fitting_zs, fitting_RAs, fitting_Decs, fitting_mu_resids, fitting_muErrs, [central_RA, central_Dec], [central_coord[0], 0.0, 0.0, params[0], params[1]], print_params = print_params, gal_dens_weight = 0.0)

        if one_d_fit:
            min_res = optimize.minimize_scalar(fit_funct_no_gal_dens, bounds = bounds)
        elif (init_guess == None) and (bounds == None and n_grid_samples == None):
            print ('I am doing a greater than 1d fit and need an initial guess.  Returning None... ')
            return [None, None, None, None]
        else:
            if bounds != None and n_grid_samples != None:
                n_params = len(bounds[0])
                param_arrs = [np.linspace(bounds[i][0], bounds[i][1], n_grid_samples[i]) for i in range(len(bounds[0]))]
                full_param_points = can.cartesian(param_arrs)
                fits_on_grid = [[] for param in full_param_points]
                for i in range(len(full_param_points)):
                    params = full_param_points[i]
                    #print ('params = ' + str(params))
                    fits_on_grid[i] = fit_funct_no_gal_dens(params)
                    #if not(print_params): print('\r' + 'Doing coarse (pre-mcmc) fit: ' + str( can.round_to_n(i / len(full_param_points) * 100, 2)) + '% done...', sep=' ', end='', flush=True)
                #fits_on_grid = [fit_funct_no_gal_dens(params) for params in full_param_points]
                fits_on_grid = np.transpose(np.reshape(fits_on_grid, n_grid_samples))
                if show_coarse_fit:
                    x_mesh = np.transpose(np.reshape([param_point[0] for param_point in full_param_points], n_grid_samples))
                    y_mesh = np.transpose(np.reshape([param_point[1] for param_point in full_param_points], n_grid_samples))
                    fig1, ax2 = plt.subplots(constrained_layout = True)
                    CS = ax2.contourf(y_mesh, x_mesh, fits_on_grid)
                    cbar = fig1.colorbar(CS)
                    plt.savefig(self.plot_save_dir + 'coarse_fits_on_grid.pdf')
                    plt.show()
                #Are we running a minimization algorithm or an emcee?  The determines if we want the max or the min
                grid_best_fit_coord = can.niceReverse(list(can.find2DMin(fits_on_grid) ) )
                grid_best_fit_params = [param_arrs[i][grid_best_fit_coord[i]] for i in range(n_params)]
                init_guess = grid_best_fit_params
            print ('init_guess = ' + str(init_guess))
            if method.lower() == 'mcmc':
                #I we run the emcee, we need to add a - sign to the fit, since the mcmc function maximizes, rather than minimizing.
                min_res = self.doMCMCFitNearParams(lambda params: -fit_funct_no_gal_dens(params), init_guess, bounds = bounds, penalization_dist_from_params = [100, 1.0], penalization_scaling_from_params = [1.0, 1.0], corner_plot_labels = [ r"$\delta_0$", r"$ \log_{10}$" + ' ' + r"$(r_c $" +' Mpc' + r"$^{-1})$"],
                                                   mcmc_corner_save_file_name = 'MCMC_corner_plot_' + '_'.join([ ['z', 'RA','Dec'][i] + str(can.round_to_n(central_coord[i], 3)) for i in range(len(central_coord)) ]) + '.pdf',
                                                   mcmc_chains_save_file_name = 'MCMC_chains_plot_' + '_'.join([ ['z', 'RA','Dec'][i] + str(can.round_to_n(central_coord[i], 3)) for i in range(len(central_coord)) ]) + '.pdf')
                min_res['fun'] = -min_res['fun']
            elif method.lower() == 'grid':
                min_res = {'x':init_guess, 'fun':fit_funct_no_gal_dens(init_guess)}
            else:
                min_res = optimize.minimize(fit_funct_no_gal_dens, init_guess, bounds = bounds, method = method  )
            for i in range(len(min_res['x'])):
                if min_res['x'][i] < bounds[i][0]: min_res['x'][i] = bounds[i][0]
                if min_res['x'][i] > bounds[i][1]: min_res['x'][i] = bounds[i][1]
        print ('min_res = ' + str(min_res))
        #fitted_resids = self.muDiff_of_z_funct(fitting_zs, fitting_RAs, fitting_Decs, fitting_zs, fitting_RAs, fitting_Decs, fitting_mu_resids, fitting_muErrs, [central_RA, central_Dec], [central_coord[0], min_res['x'][0], min_res['x'][1], 0.0, 0.0], self.fit_funct)
        null_fit_val = fit_funct_no_gal_dens(null_param_vals)
        fit_r_chi_square = self.computeRChiSqr(fitting_zs, fitting_RAs, fitting_Decs, fitting_mu_resids, fitting_muErrs, [central_RA, central_Dec], self.n_fit_params, [central_coord[0], 0.0, 0.0] + min_res['x'], print_params = 0)
        full_r_chi_square = ( np.sum([ self.all_muResids[i] ** 2.0 / self.all_mu_errs[i] ** 2.0 for i in range(len(self.all_zs)) if not (i in nearby_sn_indeces) ]) + fit_r_chi_square * (len(nearby_sn_indeces) - self.n_fit_params)) / (len(self.all_zs) - self.n_fit_params)
        null_r_chi_square = self.computeRChiSqr(fitting_zs, fitting_RAs, fitting_Decs, fitting_mu_resids, fitting_muErrs, [central_RA, central_Dec], self.n_null_fit_params, [central_coord[0], 0.0, 0.0] + null_param_vals, print_params = 0)
        #full_null_r_chi_square = np.sum([ self.all_muResids[i] ** 2.0 / self.all_mu_errs[i] ** 2.0 for i in range(len(self.all_zs)) ]) / (len(self.all_zs) - 1)
        print ('null_fit_val = ' + str(null_fit_val))
        if null_fit_val > 99:
            print ('!!! null_param_vals = ' + str(null_param_vals) + ' at coordinate ' + str(central_coord) + ' yielded null_fit_val = ' + str(null_fit_val) + ' !!!!!')
            print ('Let us redo the calculation...')
            -self.overall_log_prob_funct(fitting_zs, fitting_RAs, fitting_Decs, fitting_mu_resids, fitting_muErrs, [central_RA, central_Dec], [central_coord[0], null_param_vals[0], null_param_vals[1], 0.0, 0.0], print_params = 1, gal_dens_weight = 0.0)
        return min_res, fit_r_chi_square, null_r_chi_square, null_fit_val

    def fullSNeFitFunct(self, params, sky_start, comoving_bin, min_n_sn):
        central_coord = [params[0], sky_start[0] + params[1] / np.cos(np.deg2rad(sky_start[1] + params[2])), sky_start[1] + params[2]]
        central_RA, central_Dec = central_coord[1:]
        #print ('central_coord = ' + str(central_coord))
        nearby_sn_indeces = self.getSNWithinComovingDistance(central_coord[0], [central_coord[1], central_coord[2]], comoving_bin)
        #print ('len(nearby_sn_indeces) = ' + str(len(nearby_sn_indeces)) + ' : ' + str(nearby_sn_indeces) )
        if len(nearby_sn_indeces) < min_n_sn:
            results = self.low_log_prob_val
        else:
            #print ('nearby_sn_indeces = ' + str(nearby_sn_indeces))
            fitting_zs, fitting_RAs, fitting_Decs, fitting_mu_resids, fitting_muErrs, fitting_surveys = [[arr[i] for i in nearby_sn_indeces] for arr in [self.all_zs, self.all_RAs, self.all_Decs, self.all_muResids, self.all_mu_errs, self.all_surveys]]
            #RA_offset_to_center = (180.0 - central_coord[1])
            #fitting_RAs = [ fitting_RAs[j] + (180.0 - fitting_RAs[0]) for j in range(len(fitting_RAs)) ]
            #fitting_RAs = [RA - 360.0 if RA > 360.0 else RA + 360.0 if RA < 0.0 else RA for RA in fitting_RAs]
            #fitting_RAs = self.centerRAs(fitting_RAs,  RA_offset_to_center )
            #central_RA, central_Dec = [self.centerRAs([central_coord[1]],  RA_offset_to_center )[0],  central_coord[2]]
            #print ('[[central_RA, central_Dec],[fitting_RAs[0], fitting_Decs[0]] ]=  ' + str([[central_RA, central_Dec],[fitting_RAs[0], fitting_Decs[0]] ] ) )
            results = self.overall_log_prob_funct(fitting_zs, fitting_RAs, fitting_Decs, fitting_mu_resids, fitting_muErrs, [central_RA, central_Dec], params, print_params = 0, gal_dens_weight = 0.0)
            #fit_funct_no_gal_dens = lambda ext_params: -self.overall_log_prob_funct(fitting_zs, fitting_RAs, fitting_Decs, fitting_mu_resids, fitting_muErrs, [central_RA, central_Dec], ext_params, print_params = 0, gal_dens_weight = 0.0)
        #if results <= -100:
        #    print ('len(nearby_sn_indeces) = ' + str(len(nearby_sn_indeces)) + ' : ' + str(nearby_sn_indeces) )
        #print ('params: ' + str(params.tolist()) + ' => results: ' + str(results))
        return results

    def doFullFitAroundPoint(self, start_min_seed, sky_start, comoving_bin, min_n_sn, method = 'mcmc', bounds = None,  mcmc_save_suffix = '', null_param_vals = [0.0, 2.0]):
        #start_min_seed = [best_fit_point]
        mcmc_seed_rand_scalings = [0.0001, 0.0001, 0.0001, 0.0001]
        if method == 'mcmc':
            #local_fit = self.doMCMCFitNearParams(lambda params: self.fullSNeFitFunct(params, sky_start, comoving_bin, min_n_sn), start_min_seed, bounds = None, penalization_dist_from_params = [np.inf, np.inf, np.inf, 20, 20], penalization_scaling_from_params = [1.0, 1.0, 1.0, 1.0, 1.0],  corner_plot_labels = [r'$z$', r"$\log_{10}$" + ' ' + r"$(M$ M$_{\cdot}^{-12})$", r"$ \log_{10}$" + ' ' + r"$(r_c $" +' Mpc' + r"$^{-1})$" , r'$\Delta$ RA', r'$\Delta$ Dec'])
            local_fit = self.doMCMCFitNearParams(lambda params: self.fullSNeFitFunct(params, sky_start, comoving_bin, min_n_sn), start_min_seed, bounds = None, penalization_dist_from_params = [np.inf, np.inf, np.inf, 20, 20], penalization_scaling_from_params = [1.0, 1.0, 1.0, 1.0, 1.0],  corner_plot_labels = [r'$z$', r'$\Delta$ RA (deg)', r'$\Delta$ Dec (deg)', r'$\delta_0$', r"$ \log_{10}$" + ' ' + r"$(r_c $" +' Mpc' + r"$^{-1})$" ],mcmc_corner_save_file_name = 'MCMC_corner_plot' + mcmc_save_suffix + '.pdf', mcmc_chains_save_file_name = 'MCMC_chains_plot' + mcmc_save_suffix + '.pdf' )
            local_fit['fun'] = -local_fit['fun']

        fitted_params = local_fit['x']
        print ('fitted_params = ' + str(fitted_params))
        central_coord = [fitted_params[0], sky_start[0] + fitted_params[1] / np.cos(np.deg2rad(sky_start[1] + fitted_params[2])), sky_start[1] + fitted_params[2]]
        central_RA, central_Dec = central_coord[1:]
        #print ('central_coord = ' + str(central_coord))
        nearby_sn_indeces = self.getSNWithinComovingDistance(central_coord[0], [central_coord[1], central_coord[2]], comoving_bin)
        fitting_zs, fitting_RAs, fitting_Decs, fitting_mu_resids, fitting_muErrs, fitting_surveys = [[arr[i] for i in nearby_sn_indeces] for arr in [self.all_zs, self.all_RAs, self.all_Decs, self.all_muResids, self.all_mu_errs, self.all_surveys]]
        fit_r_chi_square = self.computeRChiSqr(fitting_zs, fitting_RAs, fitting_Decs, fitting_mu_resids, fitting_muErrs, [central_RA, central_Dec], self.n_fit_params, fitted_params, print_params = 0)
        null_r_chi_square = self.computeRChiSqr(fitting_zs, fitting_RAs, fitting_Decs, fitting_mu_resids, fitting_muErrs, [central_RA, central_Dec], self.n_null_fit_params, fitted_params[0:3] + null_param_vals , print_params = 0, )

        #null_fit_val = fit_funct_no_gal_dens(start_min_seed)
        null_fit_val = -self.overall_log_prob_funct(fitting_zs, fitting_RAs, fitting_Decs, fitting_mu_resids, fitting_muErrs, [central_RA, central_Dec], fitted_params[0:3] + [0.0, fitted_params[1]], print_params = 0, gal_dens_weight = 0.0, n_model_params = 0)
        return local_fit, fit_r_chi_square, null_r_chi_square, null_fit_val

    def doFitOnSNeSubset(self, central_coord, nearby_sn_indeces,
                         method = 'CG', n_mcmc_chains = 8, n_mcmc_steps = 2000, mcmc_seed_rand_scalings = [0.01, 1.0, 1.0, 1.0], comoving_bound = None, start_min_seed = None, null_param_vals = [0.0, 2.0], bounds = None, print_params = 0):
        local_bounds = self.getBoundsGivingMaxComovingRange(comoving_bound, central_coord, self.bounds)
        mcmc_seed_rand_scalings = [0.0001, 0.0001, 0.0001, 0.0001]
        fitting_zs, fitting_RAs, fitting_Decs, fitting_mu_resids, fitting_muErrs, fitting_surveys = [[arr[i] for i in nearby_sn_indeces] for arr in [self.all_zs, self.all_RAs, self.all_Decs, self.all_muResids, self.all_mu_errs, self.all_surveys]]
        RA_offset_to_center = (180.0 - central_coord[1])
        #fitting_RAs = [ fitting_RAs[j] + (180.0 - fitting_RAs[0]) for j in range(len(fitting_RAs)) ]
        #fitting_RAs = [RA - 360.0 if RA > 360.0 else RA + 360.0 if RA < 0.0 else RA for RA in fitting_RAs]
        fitting_RAs = self.centerRAs(fitting_RAs,  RA_offset_to_center )
        #RA_center, Dec_center = [np.mean(fitting_RAs), np.mean(fitting_Decs)]
        #start_min_seed = [np.mean(fitting_zs), 0.0, 0.0, 0.0]
        #RA_center, Dec_center = [central_coord[1], central_coord[2]]
        if start_min_seed == None:
            start_min_seed = [central_coord[0], 0.0, 0.0, 0.0]
        central_RA, central_Dec = [self.centerRAs([central_coord[1]],  RA_offset_to_center )[0],  central_coord[2]]

        fit_funct_no_gal_dens = lambda ext_params: -self.overall_log_prob_funct(fitting_zs, fitting_RAs, fitting_Decs, fitting_mu_resids, fitting_muErrs, [central_RA, central_Dec], ext_params, print_params = print_params, gal_dens_weight = 0.0)
        #mcmc_log_prob_funct = lambda params: fit_funct_no_gal_dens(params) + self.mcmc_prior_funct([params[0] - start_min_seed[0] if params[0] > 0.0 else -np.inf ] + list(params[1:]), print_params = 0)
        if method == 'mcmc':
            local_fit = self.doMCMCFitNearParams(lambda params: -fit_funct_no_gal_dens(params), start_min_seed, bounds = local_bounds, penalization_dist_from_params = [np.inf, np.inf, np.inf, 20, 20], penalization_scaling_from_params = [1.0, 1.0, 1.0, 1.0, 1.0],  corner_plot_labels = [r'$z$', r"$\log_{10}$" + ' ' + r"$(M$ M$_{\cdot}^{-12})$", r"$ \log_{10}$" + ' ' + r"$(r_c $" +' Mpc' + r"$^{-1})$" , r'$\Delta$ RA', r'$\Delta$ Dec'])
            local_fit['fun'] = -local_fit['fun']
        else:
            local_fit = optimize.minimize(fit_funct_no_gal_dens, start_min_seed, method = method, tol = 10.0 ** -5.0, bounds = local_bounds)

        fit_r_chi_square = self.computeRChiSqr(fitting_zs, fitting_RAs, fitting_Decs, fitting_mu_resids, fitting_muErrs, [central_RA, central_Dec], self.n_fit_params, [local_fit['x'][0], local_fit['x'][1], local_fit['x'][2], local_fit['x'][3], local_fit['x'][4]], print_params = 0)
        null_r_chi_square = self.computeRChiSqr(fitting_zs, fitting_RAs, fitting_Decs, fitting_mu_resids, fitting_muErrs, [central_RA, central_Dec], self.n_null_fit_params, [local_fit['x'][0], null_param_vals[0], null_param_vals[1], local_fit['x'][3], local_fit['x'][4]], print_params = 0)

        null_fit_val = fit_funct_no_gal_dens(start_min_seed)
        null_fit_val = -self.overall_log_prob_funct(fitting_zs, fitting_RAs, fitting_Decs, fitting_mu_resids, fitting_muErrs, [central_RA, central_Dec], start_min_seed, print_params = 0, gal_dens_weight = 0.0, n_model_params = 0)
        return local_fit, fit_r_chi_square, null_r_chi_square, null_fit_val

    def plotSNeFitOnSNeSubset(self, fit_params, sn_indeces_to_plot, central_3d_coord,
                             fit_bounds = None, plot_pause_time = 0.5, n_plot_fill_points = 51, figsize = [9,9],
                             full_title = '', save_name = 'SingleFieldFit_SN_.png', save_fig = 1, show_fig = 0, ref_axes = None, show_legend = 1):
        central_z = central_3d_coord[0]
        central_coord = central_3d_coord[1:]
        RA_offset_to_center = (180.0 - central_coord[0])
        central_RA, central_Dec =[self.centerRAs([central_coord[0]],  RA_offset_to_center )[0],  central_coord[1]]
        fitting_RAs, fitting_Decs, fitting_zs, fitting_mu_resids, fitting_mu_errs, fitting_mu_colors = [ [arr[i] for i in sn_indeces_to_plot] for arr in [self.all_RAs, self.all_Decs, self.all_zs, self.all_muResids, self.all_mu_errs, self.all_plot_colors] ]
        zs_to_plot, RAs_to_plot, Decs_to_plot = [ np.linspace(min(fitting_zs), max(fitting_zs), n_plot_fill_points), [central_RA for i in range(n_plot_fill_points)], [central_Dec for i in range(n_plot_fill_points)] ]
        #fitting_RAs = [ fitting_RAs[j] + (180.0 - fitting_RAs[0]) for j in range(len(fitting_RAs)) ]
        #fitting_RAs = [RA - 360.0 if RA > 360.0 else RA + 360.0 if RA < 0.0 else RA for RA in fitting_RAs]
        fitting_RAs = self.centerRAs(fitting_RAs,  RA_offset_to_center )
        RAs_to_plot = self.centerRAs(RAs_to_plot,  RA_offset_to_center )
        muDiffs_to_plot = self.muDiff_of_z_funct(zs_to_plot, RAs_to_plot, Decs_to_plot, fitting_zs, fitting_RAs, fitting_Decs, fitting_mu_resids, fitting_mu_errs, [ central_RA, central_Dec], fit_params, self.fit_funct)
        fitted_muDiffs = self.muDiff_of_z_funct(fitting_zs, fitting_RAs, fitting_Decs, fitting_zs, fitting_RAs, fitting_Decs, fitting_mu_resids, fitting_mu_errs, [ central_RA, central_Dec], fit_params, self.fit_funct)
        sne_rchi_sqr, sne_fit = self.computeChiSqrProb(fitting_zs, fitting_RAs, fitting_Decs, fitting_mu_resids, fitting_mu_errs, [ central_RA, central_Dec], 5, fit_params[0:4], print_params = 0)
        #                self.muDiff_of_z_funct(all_zs, all_RAs, all_Decs, all_zs, all_RAs, all_Decs, all_resids, all_errs, field_center, fit_params, self.fit_funct)
        #fitting_RAs = [RA - RA_offset_to_center for RA in fitting_RAs]
        if ref_axes == None:
            plt.close('all')
            fig = plt.figure(constrained_layout=True, figsize = figsize)
            gs = fig.add_gridspec(3,4)
            full_sky_ax = fig.add_subplot(gs[1:, :], projection="aitoff")
            #f, axarr = plt.subplots(1,3, figsize = figsize )
            plots_for_legend = self.makeSkyPlotOfSNe(full_sky_ax)
            if show_legend: plt.legend(plots_for_legend, self.included_surveys, ncol = 5)
            partial_sky_ax = fig.add_subplot(gs[0, 0:2])
            fit_ax = fig.add_subplot(gs[0, 2:])
            plt.suptitle( full_title  )
            subfield_coords = [[rect_coord * self.deg_to_rad for rect_coord in self.adjustRAForSkyPlot(np.array(partial_sky_ax.get_xlim())) ], [rect_coord * self.deg_to_rad for rect_coord in partial_sky_ax.get_ylim()]]
            if subfield_coords[0][1] < subfield_coords[0][0]: #indicates we wrap through 0 RA
                subfield_coords = [ [[subfield_coords[0][1], 2.0 * np.pi], subfield_coords[1]], [[0.0, subfield_coords[0][0]], subfield_coords[1]] ]
            else:
                subfield_coords = [subfield_coords]
            field_rects = [Rectangle([subfield_coord[0][0], subfield_coord[1][0]], (subfield_coord[0][1] - subfield_coord[0][0]), subfield_coord[1][1] - subfield_coord[1][0], edgecolor='black', facecolor='black', alpha = 0.25) for subfield_coord in subfield_coords ]
            [full_sky_ax.add_patch(field_rect) for field_rect in field_rects]
            full_sky_ax.scatter(np.radians(self.adjustRAForSkyPlot(np.array(self.centerRAs( [central_RA + fit_params[-2]], RA_offset_to_center, undo_centering = 1)))), np.radians(central_Dec + fit_params[-1]), marker = '+', c = 'r')
            full_sky_ax.scatter(np.radians(self.adjustRAForSkyPlot(np.array(self.centerRAs([central_RA], RA_offset_to_center, undo_centering = 1)))), np.radians(central_Dec), marker = 'x', c = 'k')
        else:
            partial_sky_ax, fit_ax = ref_axes
        sne_on_sky = partial_sky_ax.scatter(self.centerRAs(fitting_RAs,  RA_offset_to_center, undo_centering = 1 ) , fitting_Decs, c = fitting_mu_colors)
        init_z = fit_ax.axvline(central_z, linestyle = '--', color = 'k')
        fit_z = fit_ax.axvline(fit_params[0], linestyle = '--', color = 'r')
        sne_in_z = fit_ax.scatter(fitting_zs, fitting_mu_resids , c = fitting_mu_colors)
        fit_ax.errorbar(fitting_zs, fitting_mu_resids , yerr = fitting_mu_errs, fmt = 'none', color = fitting_mu_colors)
        fit_start = partial_sky_ax.scatter(self.centerRAs( [central_RA], RA_offset_to_center, undo_centering = 1), central_Dec, marker= '+', c = 'k')
        fit_end = partial_sky_ax.scatter( self.centerRAs( [central_RA + fit_params[-2]], RA_offset_to_center, undo_centering = 1), central_Dec + fit_params[-1], marker = 'x', c = 'r')
        if not(fit_bounds == None):
            print ('fit_bounds = ' + str(fit_bounds))
            print ('central_3d_coord = ' + str(central_3d_coord) )
            z_lims = fit_ax.get_xlim()
            z_lims = [min(z_lims[0], fit_bounds[0][0] - (z_lims[1] - z_lims[0]) * 0.1), max(z_lims[1], fit_bounds[0][1] + (z_lims[1] - z_lims[0]) * 0.1)]
            fit_ax.axvspan(z_lims[0], fit_bounds[0][0], color = 'gray', alpha = 0.2)
            fit_ax.axvspan(fit_bounds[0][1], z_lims[1], color = 'gray', alpha = 0.2)
            fit_ax.set_xlim(z_lims)
            RA_lims, Dec_lims = [partial_sky_ax.get_xlim(), partial_sky_ax.get_ylim()]
            RA_center, Dec_center = [central_coord[0], central_coord[1]]
            print ('[RA_lims, Dec_lims] = ' + str([RA_lims, Dec_lims] ))
            print ('fit_bounds = ' + str(fit_bounds))
            print ('[np.min([RA_lims[0], fit_bounds[2][0] + RA_center]), np.max([RA_lims[1], fit_bounds[2][1] + RA_center])] = ' + str([np.min([RA_lims[0], fit_bounds[2][0] + RA_center]), np.max([RA_lims[1], fit_bounds[2][1] + RA_center])]))
            RA_lims = [np.min([RA_lims[0], (fit_bounds[2][0] + RA_center) - (RA_lims[1] - RA_lims[0]) * 0.1]), np.max([RA_lims[1], fit_bounds[2][1] + RA_center + (RA_lims[1] - RA_lims[0]) * 0.1])]
            Dec_lims = [np.min([Dec_lims[0], fit_bounds[3][0] + Dec_center - (Dec_lims[1] - Dec_lims[0]) * 0.1]), np.max([Dec_lims[1], fit_bounds[3][1] + Dec_center + (Dec_lims[1] - Dec_lims[0]) * 0.1])]
            partial_sky_ax.axvspan(RA_lims[0], fit_bounds[2][0] + RA_center, color = 'gray', alpha = 0.2)
            partial_sky_ax.axvspan(fit_bounds[2][1] + RA_center, RA_lims[1], color = 'gray', alpha = 0.2)
            partial_sky_ax.set_xlim(RA_lims)
            partial_sky_ax.axhspan(Dec_lims[0], fit_bounds[3][0] + Dec_center, color = 'gray', alpha = 0.2)
            partial_sky_ax.axhspan(fit_bounds[3][1] + Dec_center, Dec_lims[1], color = 'gray', alpha = 0.2)
            partial_sky_ax.set_ylim(Dec_lims)

        #full_sky_ax.scatter(self.adjustRAForSkyPlot(np.radians(self.centerRAs([center_RA + fit_params[-2]], RA_offset_to_center, undo_centering = 1))), np.radians(center_Dec + fit_params[-1]), marker = '+', c = 'r')
        #fit_ax.plot(zs_to_plot, muDiffs_to_plot, c = 'k')
        fitted_mu_resids = fit_ax.scatter(fitting_zs, fitted_muDiffs, c = 'k', marker = 'x')
        fit_ax.set_xlabel(r'$z$')
        fit_ax.set_ylabel(r'$\Delta \mu$')
        partial_sky_ax.set_xlabel('RA (deg)')
        partial_sky_ax.set_ylabel('Dec (deg)')
        if show_legend:
            partial_sky_ax.legend([fit_start, fit_end, sne_on_sky], [r'$\delta$ cent. init guess', r'fit $\delta$ center', 'True SNe'], ncol = 1)
            fit_ax.legend([init_z, fit_z, sne_in_z, fitted_mu_resids], [r'$\delta$ cent. init guess', r'fit $\delta$ center', 'True SNe', r'$\Delta \mu$ from $\delta$'], ncol = 1)
        #shifted_RA_ticks = partial_sky_ax.get_xticks()
        #true_RA_tick_labels = [round(x_tick * 10) // 10 for x_tick in self.centerRAs(shifted_RA_ticks, RA_offset_to_center, undo_centering = 1) ]
        #[round(((x_tick - RA_offset_to_center) % 360) * 10) // 10 for x_tick in shifted_RA_ticks]
        #true_RA_ticks = self.centerRAs(true_RA_tick_labels, RA_offset_to_center)
        #partial_sky_ax.set_xticks(true_RA_ticks)
        #partial_sky_ax.set_xticklabels(true_RA_tick_labels)
        partial_sky_ax.invert_xaxis()
        if save_fig: plt.savefig(self.plot_save_dir + save_name)
        if show_fig: plt.pause(plot_pause_time)
        return [fit_start, fit_end, init_z, fit_z, fitted_mu_resids]


    def computeBestFitsAroundSN(self, max_comoving_sep = 300, min_sne_for_fit = 5, plot_pause_time = 0.5, n_plot_fill_points = 51, figsize = [9,9], n_mcmc_chains = 8, n_mcmc_steps = 10000, sn_indeces_to_fit_around = 'all', save_figs = 1):
        #We want to run a minimization centered around each supernova, looking only at those SNe
        # that are within some comoving radius of of the target supernova.
        fits_by_sn = [None for i in range(len(self.all_zs))]
        if sn_indeces_to_fit_around == 'all':
            sn_indeces_to_fit_around = list(range(len(self.all_zs)))
        for i in sn_indeces_to_fit_around :
        #for i in range(710, 712):
            nearby_sn_indeces = [ j for  j in list(range(0, i)) + list(range(i+1, len(self.all_zs))) if self.comoving_sne_seps[(i, j)] < max_comoving_sep ]
            n_nearby_sn = len(nearby_sn_indeces)
            print ('Within ' + str(max_comoving_sep) + 'Mpc of the ' + str(i) + 'th sn, there are the following ' + str(n_nearby_sn) + ' sn: ' + str(nearby_sn_indeces))
            if n_nearby_sn >= min_sne_for_fit:
                fitting_sn = [ self.all_sns[j] for j in [i] + nearby_sn_indeces ]
                fitting_zs = [self.all_zs[j] for j in [i] + nearby_sn_indeces]
                #if the RAs wrap around RA = 360, we need to use negative RAs
                fitting_RAs = [self.all_RAs[j] for j in [i] + nearby_sn_indeces]
                fitting_RAs = self.centerRAs(fitting_RAs, RA_offset_to_center)
                #fitting_RAs = [ fitting_RAs[j] + (180.0 - fitting_RAs[0]) for j in range(len(fitting_RAs)) ]
                #fitting_RAs = [RA - 360.0 if RA > 360.0 else RA + 360.0 if RA < 0.0 else RA for RA in fitting_RAs]
                #fitting_RAs = [self.all_RAs[i] for i in [i] + nearby_sn_indeces]
                fitting_Decs = [self.all_Decs[j] for j in [i] + nearby_sn_indeces]
                fitting_mu_resids = [self.all_muResids[j] for j in [i] + nearby_sn_indeces]
                fitting_mu_errs = [self.all_mu_errs[j] for j in [i] + nearby_sn_indeces]
                fitting_mu_colors = [self.all_sns[j]['color'] for j in [i] + nearby_sn_indeces]
                center_z, center_sn, center_RA, center_Dec, center_mu_resid, center_mu_err, center_mu_color = [ np.mean(fitting_zs), fitting_sn[0], np.mean(fitting_RAs), np.mean(fitting_Decs), fitting_mu_resids[0], fitting_mu_errs[0], fitting_mu_colors[0] ]
                fitting_zs, fitting_sn, fitting_RAs, fitting_Decs, fitting_mu_resids, fitting_mu_errs, fitting_mu_colors = can.safeSortOneListByAnother(fitting_zs, [fitting_zs, fitting_sn, fitting_RAs, fitting_Decs, fitting_mu_resids, fitting_mu_errs, fitting_mu_colors ])
                #Do fit on only the nearby supernovae
                print ('Doing fit on ' + str(i) + 'th SNe at [z, RA, Dec]: '+ str([self.all_zs[i], self.all_RAs[i], self.all_Decs[i]]) )
                #local_mcmc_sampler = self.doFitOnSNeSubset(i, nearby_sn_indeces, n_mcmc_chains = n_mcmc_chains, n_mcmc_steps = n_mcmc_steps )
                #flat_samples = local_mcmc_sampler.get_chain(discard=100, thin=15, flat=True)
                #local_res = [float(np.percentile(flat_samples[:, i], [50])) for i in range(4)]
                local_fit, null_fit_val =  self.doFitOnSNeSubset([center_z, center_RA, center_Dec], nearby_sn_indeces, n_mcmc_chains = n_mcmc_chains, n_mcmc_steps = n_mcmc_steps )
                local_res, local_fit_val = [local_fit['x'], local_fit['fun']]
                local_fit_val, null_fit_val = [ -local_fit_val, -null_fit_val ]
                computeProbImprovement  = self.computeProbImprovement( local_fit_val , null_fit_val )

                fits_by_sn[i] = [local_res, local_fit_val, null_fit_val ]
                print ('center_RA + local_res[-2] = ' + str(center_RA + local_res[-2]))
                print ('self.centerRAs([center_RA + local_res[-2]], RA_offset_to_center, undo_centering = 1) = ' + str(self.centerRAs([center_RA + local_res[-2]], RA_offset_to_center, undo_centering = 1)))
                print ('Done! - best fit has params: ' + str(local_res) + ' with fit res ' + str(local_fit_val) + ' VS null fit of ' + str(null_fit_val))

        return fits_by_sn

    def determineComovingSep(self, center_coord, ref_coord, min_sep = np.inf ):
        center_comoving = self.r_of_z_interp(center_coord[0])
        ref_comoving = self.r_of_z_interp(ref_coord[0])
        radial_comoving_sep = abs(ref_comoving - center_comoving)
        dec_ang_sep = self.deg_to_rad * abs(ref_coord[2] - center_coord[2])
        min_dec_comoving_sep = min(dec_ang_sep * center_comoving, dec_ang_sep * ref_comoving)
        #Is a lower bound on the separation greater than the minimum separation?
        # We can just return infinity then without having to do an expensive calculation.
        #print ('[radial_comoving_sep, min_sep, min_dec_comoving_sep] = ' + str([radial_comoving_sep, min_sep, min_dec_comoving_sep]))
        if radial_comoving_sep > min_sep or min_dec_comoving_sep > min_sep:
            return np.inf
        #Otherwise, we need to actually compute
        comoving_sep = self.computeComovingSepsOfCoords(center_coord, ref_coord)
        return comoving_sep

    def downselectComovingGridByOctant(self, cartesian_points, hemisphere, angle_divs, angle_slice):
        """
        Pick as seeds only points in only 1 of eight cartesian quadrant.
            Used to parallellize lengthy fitting procedures. Quadrants
            because it makes the breaking up of the cartesian grid very easy.
        """
        """
        if search_octant == 1:
            cartesian_points = [ point for point in cartesian_points if (point[0] >= 0 and point[1] >= 0 and point[2] >= 0) ]
        elif search_octant == 2:
            cartesian_points = [ point for point in cartesian_points if (point[0] < 0 and point[1] >= 0 and point[2] >= 0) ]
        elif search_octant == 3:
            cartesian_points = [ point for point in cartesian_points if (point[0] < 0 and point[1] < 0 and point[2] >= 0) ]
        elif search_octant == 4:
            cartesian_points = [ point for point in cartesian_points if (point[0] >= 0 and point[1] < 0 and point[2] >= 0) ]
        elif search_octant == 5:
            cartesian_points = [ point for point in cartesian_points if (point[0] >= 0 and point[1] < 0 and point[2] < 0) ]
        elif search_octant == 6:
            cartesian_points = [ point for point in cartesian_points if (point[0] < 0 and point[1] < 0 and point[2] < 0) ]
        elif search_octant == 7:
            cartesian_points = [ point for point in cartesian_points if (point[0] < 0 and point[1] >= 0 and point[2] < 0) ]
        elif search_octant == 8:
            cartesian_points = [ point for point in cartesian_points if (point[0] >= 0 and point[1] >= 0 and point[2] < 0) ]
        """
        if hemisphere == 0:
            cartesian_points = [point for point in cartesian_points if point[2] >= 0 ]
        elif hemisphere == 1:
            cartesian_points = [point for point in cartesian_points if point[2] < 0 ]
        angles = [(np.arctan2(point[0], point[1]) + np.pi) % (2.0 * np.pi) for point in cartesian_points]
        slice_size = 2.0 * np.pi / angle_divs
        #print ('angles = ' + str(angles))
        #print ('slice_size = ' + str(slice_size))

        cartesian_points = [cartesian_points[i] for i in range(len(cartesian_points)) if (angles[i] >= slice_size * (angle_slice - 1) and angles[i] < slice_size * angle_slice) ]

        return cartesian_points


    def determineComovingGrid(self, comoving_grid_sep_Mpc = 100, z_lims = None, hemisphere = 'both', angle_divs = 1, angle_slice = 1 ):
        """
        Make a cartesian grid of points in comoving space, returning those that lie within a spherical
            annulus of minimum/maximum redshift radius.
        Returned in spherical (z, RA, Dec) coords.
        """
        if z_lims == None:
            z_lims = self.z_range
        if hemisphere == None:
            hemisphere = self.hemisphere
        if angle_divs == None:
            angle_divs = self.angle_divs
        if angle_slice == None:
            angle_slice = self.angle_slice
        #if search_octant == None:
        #    search_octant = self.search_octant
        r_lims = [self.r_of_z_interp(z_lims[0]), self.r_of_z_interp(z_lims[1])]
        print ('r_lims = ' + str(r_lims))
        #one_d_comoving_axis = np.arange(r_lims[0] + comoving_grid_sep_Mpc * 0.1, r_lims[1], comoving_grid_sep_Mpc )
        one_d_comoving_axis = np.arange(0.0, r_lims[1], comoving_grid_sep_Mpc )
        one_d_comoving_axis = [elem for elem in one_d_comoving_axis if (elem > r_lims[0] and elem < r_lims[1])]
        one_d_comoving_axis = can.niceReverse([-elem for elem in one_d_comoving_axis]) + one_d_comoving_axis
        cartesian_points = can.flattenListOfLists(can.flattenListOfLists([[[[x,y,z] for x in one_d_comoving_axis] for y in one_d_comoving_axis] for z in one_d_comoving_axis]))
        cartesian_points = [point for point in cartesian_points if (point[0] ** 2.0 + point[1] ** 2.0 + point[2] ** 2.0 > r_lims[0] ** 2.0) and (point[0] ** 2.0 + point[1] ** 2.0 + point[2] ** 2.0 < r_lims[1] ** 2.0)]
        cartesian_points = self.downselectComovingGridByOctant(cartesian_points, hemisphere, angle_divs, angle_slice)
        rs = [(np.sqrt(point[0] ** 2.0 + point[1] ** 2.0 + point[2] ** 2.0)) for point in cartesian_points]

        xs, ys, zs_cart = [[point[0] for point in cartesian_points], [point[1] for point in cartesian_points], [point[2] for point in cartesian_points]]
        Rs = [(np.sqrt(point[0] ** 2.0 + point[1] ** 2.0)) for point in cartesian_points]
        thetas = [math.atan2(Rs[i], zs_cart[i]) for i in range(len(Rs))]
        phis = [math.atan2(xs[i], ys[i]) for i in range(len(Rs))]
        phis = [-phi + np.pi if phi < 0.0 else phi for phi in phis]
        spherical_coords = [[float(self.z_of_r_interp(rs[i])), 360.0 - phis[i] / self.astro_arch.getDegToRad(), 90.0 - thetas[i] / self.astro_arch.getDegToRad()] for i in range(len(cartesian_points))]

        return spherical_coords

    def getSNWithinComovingDistance(self, redshift, sky_coord, max_comoving_dist):
        comoving_dists = [ self.determineComovingSep([redshift] + sky_coord, [self.all_zs[i], self.all_RAs[i], self.all_Decs[i]], min_sep = max_comoving_dist)
                              for i in range(len(self.all_sns)) ]
        nearby_sn_indeces = [ i for i in range(len(comoving_dists)) if comoving_dists[i] < max_comoving_dist ]
        return nearby_sn_indeces

    def doSingleFitAroundCoord(self, seed_redshift, shell_coord, min_n_sn, angular_scale = None, comoving_bin = None, comoving_bound = None, init_guess = None, bounds = None, n_grid_samples = None, method = 'mcmc', one_d_fit = 0, print_params = 0, show_coarse_fit = 0):
        if comoving_bin == None:
            if angular_scale == None:
                angular_scale = self.binning_angular_scale
            comoving_bin = self.getComovingCrossSectionOfAngularScaleAtRedshift(seed_redshift, angular_scale)
        seed_comoving = self.r_of_z_interp(seed_redshift)
        #print ('[comoving_bin, seed_comoving] = ' + str([comoving_bin, seed_comoving] ))
        comoving_bin = min([comoving_bin, seed_comoving]) #Don't sweep up SNe that are on the other side of us as observers; that complicates the analysis.
        #print ('comoving_bing = ' + str(comoving_bin))
        nearby_sn_indeces = self.getSNWithinComovingDistance(seed_redshift, shell_coord, comoving_bin)
        if len(nearby_sn_indeces) >= min_n_sn:
            #print ( 'Within ' + str(comoving_bin) + ' Mpc of the coordinate: ' + str([seed_redshift] + shell_coord) + ', we have the following ' + str(len(nearby_sn_indeces)) + ' sn ' + str(nearby_sn_indeces) )
            fit_around_coord, full_r_chi_square, full_null_r_chi_square, null_prob = self.doStaticFitOnSNeSubset([seed_redshift] + shell_coord, nearby_sn_indeces, one_d_fit = one_d_fit, init_guess = init_guess, bounds = bounds, method = method, n_grid_samples = n_grid_samples, print_params = print_params, show_coarse_fit = show_coarse_fit)
            fit_param, fit_prob = [fit_around_coord['x'], fit_around_coord['fun']]
            prob_improvement = self.computeProbImprovement( -fit_prob , -null_prob)
            print ('[seed_redshift, shell_coord] = ' + str([seed_redshift, shell_coord]))
            print ('[fit_param, fit_prob, null_prob, prob_improvement] = ' + str([fit_param, fit_prob, null_prob, prob_improvement]))
            return [fit_param, full_r_chi_square, full_null_r_chi_square, prob_improvement]
        else:
            #print ('For [z, RA, Dec] = ' + str([seed_redshift] + shell_coord) + ', there are ' + str(len(nearby_sn_indeces)) + ' < min N SN ' + str(min_n_sn) + ' within ' + str(comoving_bin) + ' Mpc ')
            return None



    def doFullFitAroundCoord(self, seed_redshift, shell_coord, comoving_bin, min_n_sn, comoving_bound = None, method = 'CG'):
        nearby_sn_indeces = self.getSNWithinComovingDistance(seed_redshift, shell_coord, comoving_bin)
        if len(nearby_sn_indeces) >= min_n_sn:
            #print ( 'Within ' + str(comoving_bin) + 'Mpc of the coordinate: ' + str([seed_redshift] + shell_coord) + ', we have the following ' + str(len(nearby_sn_indeces)) + ' sn ' + str(nearby_sn_indeces) )

            fit_around_coord, null_prob = self.doFitOnSNeSubset([seed_redshift] + shell_coord, nearby_sn_indeces, method = method, comoving_bound = comoving_bound)
            fit_params, fit_prob = [fit_around_coord['x'], fit_around_coord['fun']]
            fit_bounds = self.getBoundsGivingMaxComovingRange(comoving_bound, [seed_redshift] + shell_coord, self.bounds)
            prob_improvement = self.computeProbImprovement( -fit_prob , -null_prob)
            print ('For [z, RA, Dec] = ' + str([seed_redshift] + shell_coord) + ', best fit params and prob improvement were: = ' + str([fit_params, prob_improvement]))
            fitted_plots = self.plotSNeFitOnSNeSubset(fit_params, nearby_sn_indeces, [seed_redshift] + shell_coord,
                                     fit_bounds = fit_bounds, full_title = 'Best fit params: ' + str([can.round_to_n(param, 3) for param in fit_params]) +  '; prob ' + str(can.round_to_n(prob_improvement, 3))  + ' times more likely than null',
                                     save_name = 'SingleFieldFit_z' + str(can.round_to_n(seed_redshift, 3)) + '_RA' + str(can.round_to_n(shell_coord[0], 4)) + '_RA' + str(can.round_to_n(shell_coord[1], 4)) +  '_ComovRadius' + str(comoving_bin) + '.pdf' )
            return [fit_params, -fit_prob, -null_prob]
        else:
            #print ('For [z, RA, Dec] = ' + str([seed_redshift] + shell_coord) + ', there are ' + str(len(nearby_sn_indeces)) + ' < min N SN ' + str(min_n_sn))
            return None

    def computeProbImprovement(self, best_log_prob, null_log_prob):
        prob_improvement = 10.0 ** (best_log_prob - null_log_prob )
        return prob_improvement

    def computeFitsOnGrid(self, grid_sep_comoving = 150, comoving_bin = 250, min_n_sn = 10, redshift_edges = [0.0, 1.0], coord_edges = [[-360, 720], [-90.0, 90.0]], method = 'CG',
                          save_fits = 1, save_fits_file_name = None, save_dir =  None,
                          init_param_guess = None ):
        #To run a grid on a field, do something like:
        # >>> fitter = mpfc.PanStarsFieldManager(1, sn_data_type = 'real', randomize_each_field = 1)
        # >>> fitter.computeFitsOnGrid(grid_sep_comoving = 150, comoving_bin = 250, min_n_sn = 10, redshift_edges = [0.0, 1.0], coord_edges = [[-360, 720], [-90.0, 90.0]])
        if save_dir == None:
            save_dir = self.dir_archive.getFieldFitsDir()
        if redshift_edges[0] <= 0.0:
            redshift_edges[0] = self.z_of_r_interp(grid_sep_comoving)
        comoving_sphere_radii = [self.r_of_z_interp(redshift_edges[0])]
        while comoving_sphere_radii[-1] <  self.r_of_z_interp(redshift_edges[1]):
            comoving_sphere_radii = comoving_sphere_radii + [ comoving_sphere_radii[-1] + grid_sep_comoving ]
        comoving_shell_areas = [4.0 * np.pi * rad ** 2.0 for rad in comoving_sphere_radii]
        comoving_sky_area = np.pi * grid_sep_comoving ** 2.0
        n_points_by_shells = [ int(area / comoving_sky_area) for area in comoving_shell_areas ]
        shell_coords = [can.getPointsOnSphere(n_points, return_astro_coords = 1) for n_points in n_points_by_shells]
        print ('n_points_by_shells = ' + str(n_points_by_shells))
        sky_fits = [ [ [] for j in range(len(shell_coords[i])) ] for i in range(len(comoving_sphere_radii)) ]
        n_points = np.sum(n_points_by_shells )
        n_iteration = 0
        start_time = time.time()
        for i in range(len(comoving_sphere_radii)):
            comoving_rad = comoving_sphere_radii[i]
            seed_redshift = float(self.z_of_r_interp(comoving_rad))
            for j in range(len(shell_coords[i])):
                if n_iteration % 25 == 25 - 1:
                    print ('!!!!! Working on coordinate ' + str(n_iteration) + ' of ' + str(n_points))
                    current_time = time.time()
                    elapsed_time = (current_time - start_time)/3600
                    elapsed_rate = elapsed_time / n_iteration
                    print ('We have spent ' + str(elapsed_time ) + 'h so far.  At that rate, we have ' + str((n_points - n_iteration) * elapsed_rate) + 'h to go...' )
                shell_coord = shell_coords[i][j]
                if ( (shell_coord[0] > coord_edges[0][0]) & (shell_coord[0] < coord_edges[0][1]) & (shell_coord[1] > coord_edges[1][0]) & (shell_coord[1] < coord_edges[1][1]) ) :
                    local_fit = self.doFullFitAroundCoord( seed_redshift, shell_coord, comoving_bin, min_n_sn, method = method, comoving_bound = grid_sep_comoving, init_guess = init_param_guess )
                    sky_fits[i][j] = [[seed_redshift] + shell_coord, local_fit]
                n_iteration = n_iteration + 1
        plt.close('all')
        end_time = time.time()
        print ('Full fit across all sky took ' + str((end_time - start_time) / 3600) + 'h')
        flat_sky_fits  = can.flattenListOfLists(sky_fits)
        flat_sky_fits = [fit for fit in flat_sky_fits if fit != []]
        self.sky_fits = flat_sky_fits
        if save_fits:
            if save_fits_file_name == None:
                save_fits_file_name = 'PantheonSkyGridFits_Randomize' + str(self.randomize_all) + '_ComovSep' + str(grid_sep_comoving) + '_ComovBin' + str(comoving_bin) + '_MinSNe' + str(min_n_sn) + '_RedshiftEdges' + str(redshift_edges[0]) + '_' + str(redshift_edges[1]) +  '_RAEdges' + str(coord_edges[0][0]) + '_' + str(coord_edges[0][1]) + '_DecEdges' + str(coord_edges[1][0]) + '_' + str(coord_edges[1][1]) + '.txt'
            self.saveSkyFits(save_fits_file_name, save_dir = save_dir, fits_to_save = self.sky_fits,
                           sky_fit_super_params = [grid_sep_comoving, comoving_bin, min_n_sn, redshift_edges, coord_edges ],
                           sky_fit_super_param_strs = ['grid_sep_comoving', 'comoving_bin', 'min_n_sn', 'redshift_edges', 'coord_edges' ])

        return flat_sky_fits

    def plotSkyFits(self, save_file_name_root, sky_fits = None, sky_fits_file = None,
                    load_dir = None, save_dir = None,
                    fig_unit_size = [3,3], n_cols = 8, max_fits_per_plot = 31, save_fig = 1, show_fig = 0, save_file_type = '.pdf'):
        if save_dir == None:
            save_dir = self.dir_archive.getDensityFitsDir()
        if load_dir == None:
            load_dir = self.dir_archive.getDensityFitsDir()
        if sky_fits == None and sky_fits_file == None:
            print ('I need either fits to the sky or a file with sky fits to read in! ')
            return 0
        elif sky_fits == None:
            sky_fits = self.loadSkyFits(sky_fits_file, load_dir = load_dir)
        fit_super_params = sky_fits[0]
        grid_sep_comoving, comoving_bin, min_n_sn, redshift_edges, coord_edges = fit_super_params[1:]
        grid_sep_comoving, comoving_bin, min_n_sn, redshift_edges, coord_edges = [float(grid_sep_comoving), float(comoving_bin), float(min_n_sn), can.recursiveStrToListOfLists(redshift_edges), can.recursiveStrToListOfLists(coord_edges)]
        sky_fit_vals = sky_fits[1]
        good_fits = [fit for fit in sky_fit_vals if fit[3].strip() != 'None']
        fit_prob_improvements = [float(fit[-1]) for fit in good_fits]
        print ('fit_prob_improvements = ' + str(fit_prob_improvements))
        sorted_prob_improvements, sorted_good_fits = can.safeSortOneListByAnother(fit_prob_improvements, [fit_prob_improvements, good_fits])
        sorted_prob_improvements, sorted_good_fits = [can.niceReverse(sorted_prob_improvements), can.niceReverse(sorted_good_fits)]
        print ('sorted_prob_improvements = ' + str(sorted_prob_improvements))
        n_good_fits = len(good_fits)
        n_figs = int(np.ceil((n_good_fits + 1) / max_fits_per_plot))
        all_plotted_surveys = []
        for fig_num in range(n_figs):
            print ('Making plot ' + str(fig_num + 1) + ' of ' + str(n_figs))
            set_indeces = [fig_num * max_fits_per_plot, min((fig_num + 1) * max_fits_per_plot, n_good_fits)]
            n_fits_to_plot = set_indeces[1] - set_indeces[0]
            good_fits_to_plot = sorted_good_fits[set_indeces[0]:set_indeces[1]]
            prob_improvements_to_plot = sorted_prob_improvements[set_indeces[0]:set_indeces[1]]
            n_rows = int(np.ceil(((n_fits_to_plot) * 2 + 1)/ n_cols))
            figsize = [fig_unit_size[0] * n_cols, fig_unit_size[1] * n_rows]
            f_fits, axarr_fits = plt.subplots( n_rows, n_cols, squeeze = 'False', figsize = figsize )
            for i in range(n_fits_to_plot):
                row_index, col_index = [(2 * i) // n_cols, (2 * i) % n_cols]
                sky_sect_axis, redshift_sect_axis = [axarr_fits[row_index, col_index], axarr_fits[row_index, col_index + 1]]
                good_fit = good_fits_to_plot[i]
                center_z, center_RA, center_Dec = [float(good_fit[0]), float(good_fit[1]), float(good_fit[2])]

                nearby_sn_indeces = self.getSNWithinComovingDistance(center_z, [center_RA, center_Dec], comoving_bin)
                plotted_surveys = [self.all_surveys[j] for j in nearby_sn_indeces]
                nearby_RAs, nearby_Decs, nearby_zs = [ [self.all_RAs[j] for j in nearby_sn_indeces], [self.all_Decs[j] for j in nearby_sn_indeces], [self.all_zs[j] for j in nearby_sn_indeces] ]
                good_fit_params = can.recursiveStrToListOfLists(good_fit[3], elem_type_cast = float)
                fit_prob, fit_null_prob, fit_prob_improvement = [float(good_fit[4]), float(good_fit[5]), float(good_fit[6])]
                redshift_sect_axis.text(0.05, 0.95, 'Fit prob ' + str(can.round_to_n(fit_prob, 3)) + '\n' +
                                                          'Null prob ' + str(can.round_to_n(fit_null_prob, 3)) + '\n' +
                                                          'Prob improvement ' + str(can.round_to_n(fit_prob_improvement, 3)),
                                                           horizontalalignment='left', verticalalignment='top', transform = redshift_sect_axis.transAxes )
                fit_plots_for_legend = self.plotSNeFitOnSNeSubset(good_fit_params, nearby_sn_indeces, [center_z, center_RA, center_Dec],
                                           ref_axes = [sky_sect_axis, redshift_sect_axis], show_fig = 0, save_fig = 0, show_legend = 0  )
                all_plotted_surveys = np.unique([all_plotted_surveys + plotted_surveys]).tolist()
            axarr_fits[-1,-1].legend(fit_plots_for_legend, ['Fit sky coord start', 'Fit sky coord end', 'Fit redshift start', 'Fit redshift start', 'Mu resids from fit'])
            axarr_fits[-1,-1].set_xticks([])
            axarr_fits[-1,-1].set_yticks([])
            plots_of_included_sn = [plt.scatter([], [], c = self.survey_to_color_dict[survey], marker = 'o') for survey in all_plotted_surveys]
            axarr_fits[-1,-2].legend(plots_of_included_sn, all_plotted_surveys)
            axarr_fits[-1,-2].set_xticks([])
            axarr_fits[-1,-2].set_yticks([])
            plt.tight_layout()
            if save_fig:
                plt.savefig( save_dir + save_file_name_root + '_' + str(fig_num) + save_file_type)
            if show_fig:
                plt.show()
            else:
                plt.close('all')
        return 1


    def loadSkyFits(self, sky_fits_file, delimiter = '|', n_header_lines = 3, super_fit_params = 2, load_dir = None):
        if load_dir == None:
            load_dir = self.dir_archive.density_fits_dir()
        all_sky_lines = can.readInColumnsToList(sky_fits_file, file_dir = load_dir, n_ignore = 0, delimiter = delimiter)
        sky_fits = all_sky_lines[n_header_lines:]
        super_fit_params = all_sky_lines[super_fit_params - 1]
        return [super_fit_params, sky_fits]

    def saveSkyFits(self, save_file, fits_to_save = None, sky_fit_super_params = None, sky_fit_super_param_strs = None, save_dir = None, delimiter = ' | '):
        if fits_to_save == None:
            fits_to_save = self.sky_fits
        if load_dir == None:
            load_dir = self.dir_archive.density_fits_dir()
        print ('fits_to_save = ' + str(fits_to_save) )
        if fits_to_save == None:
            print ('No fits to save! ')
            return 0
        else:

            print (('None' if sky_fit_super_param_strs == None else delimiter.join(sky_fit_super_param_strs)))
            print ( ('None' if sky_fit_super_param_strs == None else (delimiter.join([str(param) for param in sky_fit_super_params]))) )
            header = ['Sky Fitting Param IDs' + delimiter + ('None' if sky_fit_super_param_strs == None else delimiter.join(sky_fit_super_param_strs)),
                      'Sky Fitting Param Vals' + delimiter + ('None' if sky_fit_super_param_strs == None else (delimiter.join([str(param) for param in sky_fit_super_params]))),
                      delimiter.join([param_str for param_str in ['zCenter', 'RACenter', 'DecCenter', 'FitParams', 'BestFitVal', 'NullFitVal', 'ProbImprovement']]) ]
            z_seed_col = [can.round_to_n(fit[0][0], 3) for fit in fits_to_save ]
            RA_seed_col = [can.round_to_n(fit[0][1], 8) for fit in fits_to_save ]
            Dec_seed_col = [can.round_to_n(fit[0][2], 8) for fit in fits_to_save ]
            fit_params_col = [ '[' + ','.join([str(elem) for elem in fit[1][0]]) + ']' if fit[1] != None else None for fit in fits_to_save ]
            fit_prob_col = [ fit[1][1] if fit[1] != None else None for fit in fits_to_save ]
            null_prob_col = [ fit[1][2] if fit[1] != None else None for fit in fits_to_save ]
            prob_improvement_col = [self.computeProbImprovement(fit_prob_col[i], null_prob_col[i]) if fit_prob_col[i] != None and null_prob_col[i] != None else None for i in range(len(fit_prob_col))]
            can.saveListsToColumns([z_seed_col, RA_seed_col, Dec_seed_col, fit_params_col, fit_prob_col, null_prob_col, prob_improvement_col], save_file, save_dir, sep = delimiter, header = header)
            return 1


    def computeByHandFitsByField(self, field_nums, n_mcmc_steps = None, mcmc_step_sizes = None, n_chains = 'all', method = 'CG', add_to_master_min_dict = 1, update_field_fits_as_best_fit = 1, update_null_chi_by_field = 1, start_comoving_rad_frac = None,  frac_seeds_to_search_for_gals = 0.1, max_min_iterations = 200):
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
            all_zs = can.flattenListOfLists(all_zs_in_field_by_survey)
            print ('all_zs = '+ str(all_zs))
            min_z, max_z = [min(all_zs), max(all_zs)]
            all_muResids = can.flattenListOfLists(all_muResids_in_field_by_survey)
            all_muErrs = can.flattenListOfLists(all_muErrs_in_field_by_survey)
            all_RAs = can.flattenListOfLists(all_RAs_in_field_by_survey)
            all_Decs = can.flattenListOfLists(all_Decs_in_field_by_survey)
            field = self.fields[field_num]
            field_center = self.field_centers[field_num]
            print ('!!!! Starting working on [field_num, field, field_center] = ' + str([field_num, field, field_center]))
            #Note: when doing the fit based on a logarithmic probability, we need to minimize by the negative of the probability (largest probs give best fit)
            fit_funct_no_gal_dens = lambda ext_params: -self.overall_log_prob_funct(all_zs, all_RAs, all_Decs, all_muResids, all_muErrs, field_center, ext_params, print_params = 0, gal_dens_weight = 0.0)
            null_fit_funct_no_gal_dens = lambda ext_params: -self.overall_log_prob_funct(all_zs, all_RAs, all_Decs, all_muResids, all_muErrs, field_center, ext_params, print_params = 0, gal_dens_weight = 0.0, n_model_params = 0)
            fit_funct_with_gal_dens = lambda ext_params: -self.overall_log_prob_funct(all_zs, all_RAs, all_Decs, all_muResids, all_muErrs, field_center, ext_params, print_params = 0, gal_dens_weight = self.gal_dens_weighting )
            no_gal_min_res = [0 for seed in range(len(mcmc_start_params)) ]
            #Make a figure to show the best fit after every iteration
            fits_fig = plt.figure(figsize = [5,10])
            axarr = fits_fig.subplots(2)
            null_prob = fit_funct_no_gal_dens(mcmc_start_params[0])
            null_chi_sqr = self.computeRChiSqr(all_zs, all_RAs, all_Decs, all_muResids, all_muErrs, field_center, self.n_fit_params , mcmc_start_params[0], print_params = 0)
            best_chi_sqr = null_chi_sqr
            best_prob = null_prob

            for min_seed_index in range(len(mcmc_start_params)):
                seed_start = time.time()
                start_min_seed = mcmc_start_params[min_seed_index]
                #start_comoving_rad = start_comoving_rad_frac * self.r_of_z_interp(max(min_seed[0], 0.001))
                #min_seed = min_seed + [start_comoving_rad]
                if (start_min_seed[0] > min_z and start_min_seed[0] < max_z):
                    print ('Starting minimization with min_seed: ' + str(start_min_seed))
                    #no_gal_new_min_res = optimize.minimize(fit_funct_no_gal_dens, start_min_seed, method = method, tol = 10.0 ** -4.0, bounds = self.bounds)
                    no_gal_new_min_res = optimize.minimize(fit_funct_no_gal_dens, start_min_seed, method = method, tol = 10.0 ** -4.0 , options = {'maxiter':max_min_iterations})
                    #print ('new_min_res = ' + str(new_min_res))
                    fit_prob = no_gal_new_min_res['fun']
                    fit_chi_sqr = self.computeRChiSqr(all_zs, all_RAs, all_Decs, all_muResids, all_muErrs, field_center, self.n_fit_params, no_gal_new_min_res['x'], print_params = 0)
                    print ('Found minimizing params without including galaxy density: ' + str(no_gal_new_min_res['x']) + ' with chi_sqr, chi_sqr likelihood: ' + str(best_chi_sqr) + ', ' + str(fit_funct_no_gal_dens(no_gal_new_min_res['x']))  )
                    no_gal_min_res[min_seed_index] = [no_gal_new_min_res['x'], no_gal_new_min_res['fun']]
                    print ('no_gal_min_res[min_seed_index] = ' + str(no_gal_min_res[min_seed_index]))
                    #if self.gal_dens_weighting > 0.0:
                    #    print ('Starting minimization with that min_seed, now including galaxy density: ' + str(no_gal_new_min_res['x']))
                    #    new_min_res = optimize.minimize(fit_funct_with_gal_dens, no_gal_new_min_res['x'], method = method, tol = 10.0 ** -4.0, bounds = self.bounds)
                    #    print ('Found minimizing params with including galaxy density: ' + str(new_min_res['x']) + ' with chi_sqr_likelihood, gal_n_dens_likelihood, and overall likelihoods: ' + str([fit_funct_no_gal_dens(new_min_res['x']) * (1.0 - self.gal_dens_weighting), fit_funct_with_gal_dens(new_min_res['x']) - fit_funct_no_gal_dens(new_min_res['x']) * (1.0 - self.gal_dens_weighting) , fit_funct_with_gal_dens(new_min_res['x'])]))
                    #else:
                    #    new_min_res = no_gal_new_min_res
                    #min_res[min_seed_index] = [new_min_res['x'], new_min_res['fun']]
                else:
                    print ('Not going to do fit for seed ' + str(start_min_seed) + ' because seed central z is outside the z range in this field. ')
                    fit_chi_sqr = null_chi_sqr
                    fit_prob = null_prob
                    no_gal_min_res[min_seed_index] = [start_min_seed, fit_funct_no_gal_dens(start_min_seed)]
                if fit_prob < best_prob:
                    print ('[fit_prob, best_prob] = ' + str([fit_prob, best_prob]) + '.  Reassigning...')
                    best_prob = fit_prob

                seed_end = time.time()
                print ('Single seed no-galaxy  fit took ' + str(seed_end - seed_start) + 's')
                axarr[0].clear()
                axarr[1].clear()
                zs_to_plot = np.linspace(self.min_z, self.max_z, 1001).tolist( )
                RAs_to_plot =  [field_center[0] for z in zs_to_plot]
                Decs_to_plot = [field_center[1] for z in zs_to_plot]
                muDiffs_to_plot = self.muDiff_of_z_funct(zs_to_plot, RAs_to_plot, Decs_to_plot, all_zs, all_RAs, all_Decs, all_muResids, all_muErrs, field_center, no_gal_min_res[min_seed_index][0], self.fit_funct)
                computed_muDiffs = self.muDiff_of_z_funct(all_zs, all_RAs, all_Decs, all_zs, all_RAs, all_Decs, all_muResids, all_muErrs, field_center, no_gal_min_res[min_seed_index][0], self.fit_funct)
                #fitted_resids = self.muDiff_of_z_funct(plot_zs, plot_zs * 0.0, plot_zs * 0.0, plot_zs, plot_zs * 0.0, plot_zs * 0.0, plot_zs * 0.0, plot_zs * 0.0, field_center, no_gal_min_res[min_seed_index][0], self.fit_funct)
                print ('muDiffs_to_plot[0:10] = ' + str(muDiffs_to_plot[0:10]))
                print ('[zs_to_plot, RAs_to_plot, Decs_to_plot, all_zs, all_RAs, all_Decs, all_muResids, all_muErrs, field_center, no_gal_min_res[min_seed_index][0]] = ' + str([zs_to_plot, RAs_to_plot, Decs_to_plot, all_zs, all_RAs, all_Decs, all_muResids, all_muErrs, field_center, no_gal_min_res[min_seed_index][0]]))
                axarr[0].scatter(all_zs, all_muResids, c = 'blue')
                axarr[0].errorbar(all_zs, all_muResids, yerr = all_muErrs, fmt = 'none', c = 'blue')
                axarr[0].plot(zs_to_plot, muDiffs_to_plot, c = 'k')
                axarr[0].scatter(all_zs, computed_muDiffs , c = 'k', marker = 'x')
                axarr[0].axvline(start_min_seed[0], c = 'r', linestyle = '--')
                axarr[0].axvline(no_gal_min_res[min_seed_index][0][0], c = 'g', linestyle = '--')
                axarr[0].set_xlabel(r'$z$')
                axarr[0].set_ylabel(r'$\Delta \mu$ (mag)')
                axarr[0].set_ylim(-0.5, 0.5)
                axarr[0].text(0.1, 0.85, r'log $L_0=$' + str(can.round_to_n(-null_prob, 3)) + '\n' + r'log $L_{fit}=$' + str(can.round_to_n(-fit_prob, 3)) + '\n' + r'log $L_{best}=$' + str(can.round_to_n(-best_prob, 3)), transform=axarr[0].transAxes, )
                axarr[1].scatter(all_RAs, all_Decs, c = 'blue')
                axarr[1].scatter(field_center[0] + start_min_seed[2], field_center[1] + start_min_seed[3], marker = 'x', c = 'r')
                axarr[1].scatter(field_center[0] + no_gal_min_res[min_seed_index][0][2], field_center[1] + no_gal_min_res[min_seed_index][0][3], marker = 'x', c = 'g')
                axarr[1].set_xlabel('RA (deg)')
                axarr[1].set_ylabel('Dec (deg)')
                fits_fig.canvas.draw()
                plt.pause(1.0)

            if self.gal_dens_weighting > 0.0:
                min_res = [0 for seed_num in range(int(len(no_gal_min_res) * frac_seeds_to_search_for_gals))]
                no_gal_min_values = [single_seed_res[1] for single_seed_res in no_gal_min_res]
                indeces_sorted_by_min_val, no_gal_sorted_min_values = can.safeSortOneListByAnother(no_gal_min_values, [list(range(len(no_gal_min_values))), no_gal_min_values])
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
            #print ('min_res = ' + str(min_res))

            best_fit_min_res = min_res[np.argmin( [min_res[j][1] for j in range(len(min_res)) ]) ]
            best_fit_params = best_fit_min_res[0]
            best_fit_val = best_fit_min_res[1]
            r_chi_square = self.computeRChiSqr(all_zs, all_RAs, all_Decs, all_muResids, all_muErrs, field_center, self.n_fit_params , best_fit_params, print_params = 0)
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
            #print ('self.fit_params_by_field = ' + str(self.fit_params_by_field))
            for key in self.min_res.keys():
                print('key = ' + str(key))
                #print ('self.min_res[key] = ' + str(self.min_res[key]))
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
    def saveFitInformation(self, save_file_name, save_dir = None):
        if save_dir == None:
            save_dir = self.dir_archive.getFieldFitsDir()
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
        can.saveListToFile(rows_to_save, save_file_name, save_dir = save_dir, sep = ', ', append = False, header = header)
        return 1

    def fit_funct(self, H0, OmM):
        H_of_z_fitted = lambda z, H0, OmM: (H0 * np.sqrt(OmM * (1.0 + z) ** 3.0 + (self.Om0 - OmM) + self.OmR * (1.0 + z) ** 4.0 ))
        speedol = self.astro_arch.getc()
        fit_model = [25.0 + 5.0 * np.log10((1.0 + z) * speedol * integrate.quad(lambda z_int: 1.0 / H_of_z_fitted(z_int, H0, OmM), 0, z)[0]) for z in self.all_zs]
        chi_square = np.sum((np.array(fit_model) - np.array(self.all_mus)) ** 2.0 / np.array(self.all_mu_errs) ** 2.0)
        if np.isnan(chi_square): chi_square = np.inf
        return chi_square

    def doCosmicFit(self):
        fit_res = optimize.minimize(lambda params: self.fit_funct(*params), [100.0, 0.5 ], bounds = [[10.0, 200.0], [0.0, 1.0]])
        best_H0, best_OmM = fit_res['x']
        best_OmL = self.Om0 - best_OmM
        return [best_OmM, best_OmL, best_H0]

    """
    def computeComovingSepsOfCoords(self, coord1, coord2):
        zs = [coord1[0], coord2[0]]
        RAs = [coord1[1], coord2[1]]
        Decs = [coord1[2], coord2[2]]
        comoving_dists = [self.r_of_z(zs[0]), self.r_of_z(zs[1])]
        ang_sep = can.measureAngularSeparationOnSky([RAs[0], Decs[0]], [RAs[1], Decs[1]], return_radian = 1)
        comoving_center =  0.5 * (comoving_dists[0] ** 2.0 + 2.0 * comoving_dists[0] * comoving_dists[1] * np.cos(ang_sep) + comoving_dists[1] ** 2.0) ** 0.5
        comoving_parallel = 0.5 * abs(comoving_dists[0] ** 2.0 - comoving_dists[1] ** 2.0) / comoving_center
        comoving_perp = np.abs(np.sin(ang_sep)) * comoving_dists[0] * comoving_dists[1] / comoving_center
        comoving_sep = np.sqrt(comoving_parallel ** 2.0 + comoving_perp ** 2.0)

        return comoving_sep
    """

    #I think the previous calculation of comoving separation was
    #   unnecessarily complicated - comoving coordinates are
    #   fixed, so we can use normal geometric laws (since we
    #   assume a flat cosmology)...right?
    def computeComovingSepsOfCoords(self, coord1, coord2):
        zs = [coord1[0], coord2[0]]
        RAs = [coord1[1], coord2[1]]
        Decs = [coord1[2], coord2[2]]
        comoving_dists = [self.r_of_z(zs[0]), self.r_of_z(zs[1])]
        ang_sep = can.measureAngularSeparationOnSky([RAs[0], Decs[0]], [RAs[1], Decs[1]], return_radian = 1)
        comoving_sep = self.r_well_from_geometry(comoving_dists[0], comoving_dists[1], ang_sep)
        return comoving_sep

    def computeComovingSepsOfSNe(self, sn_indeces):
        zs, RAs, Decs = [ [self.all_zs[sn_indeces[0]], self.all_zs[sn_indeces[1]]], [self.all_RAs[sn_indeces[0]], self.all_RAs[sn_indeces[1]]], [self.all_Decs[sn_indeces[0]], self.all_Decs[sn_indeces[1]]] ]
        comoving_sep = self.computeComovingSepsOfCoords([zs[0], RAs[0], Decs[0]], [zs[1], RAs[1], Decs[1]])
        return comoving_sep

    def determineMinimumSpacingComovingBinningByNSNe(self, min_n_sn, verbose = 1):
        all_pairs = can.flattenListOfLists([ [(i, j) for j in range(i+1, len(self.all_sns)) ] for i in range(len(self.all_sns)) ])
        comoving_seps_by_pairs = {}
        if verbose:
            print ('Computing comoving sep between all SNe pairs...')
        for i in range(len(all_pairs)):
           if i % 1000 == 999 and verbose: print('\r' + str( can.round_to_n(i / len(all_pairs) * 100, 3)) + '% done...', sep=' ', end='', flush=True)
           pair = all_pairs[i]
           coord1 = [self.all_zs[pair[0]], self.all_RAs[pair[0]], self.all_Decs[pair[0]]]
           coord2 = [self.all_zs[pair[1]], self.all_RAs[pair[1]], self.all_Decs[pair[1]]]
           comoving_seps_by_pairs[pair] = self.computeComovingSepsOfCoords(coord1, coord2)
        if verbose:
            print ('Done computing comoving sep between all SNe pairs!')
        tightest_cluster = can.findTightestGroupingFromPairs(comoving_seps_by_pairs, min_n_sn )
        return tightest_cluster[1]

    def doFitsOnGrid(self, fit_comoving_density, min_n_sn, z_range = None, angular_scale = None, comoving_bin = None, one_d_fit = 0, init_guess = None, bounds = None, n_grid_samples = None, hemisphere = None, angle_divs = None, angle_slice = None ):
        """
        Do a fit in one dimension (over/under density mass) at a grid of points
            in comoving space.  The three passable paramters are: the density with
            which seed points are seeded in comoving space (fit_comoving_density),
            the size of the comoving sphere around each point to include SNe in
            the fit (comoving_bin), and the minimum number of SNe in the comoving
            bin to make a fit worth running (min_n_sn).
        The fits done at each grid location are performed only if some minimum
            number of SNe are sufficiently close to that seed location.
        """
        if angular_scale == None:
            angular_scale = self.binning_angular_scale
        if init_guess == None:
            init_guess = self.param_init_guess
        if bounds == None:
            bounds = self.bounds
        if hemisphere == None:
            hemisphere = self.hemisphere
        if angle_divs == None:
            angle_divs = self.angle_divs
        if angle_slice == None:
            angle_slice = self.angle_slice
        #if search_octant == None:
        #    search_octant = self.search_octant
        spherical_points = self.determineComovingGrid(comoving_grid_sep_Mpc = fit_comoving_density, z_lims = z_range, hemisphere = hemisphere, angle_divs = angle_divs, angle_slice = angle_slice)
        #print ('spherical_points = ' + str(spherical_points))
        fits = [None for i in spherical_points]
        #print ('spherical_points = ' + str(spherical_points) )
        start_time = time.time()
        prev_time = time.time()
        for i in range(len(spherical_points)):
            spherical_point = spherical_points[i]
            fits[i] = self.doSingleFitAroundCoord(spherical_point[0], spherical_point[1:], min_n_sn, angular_scale = angular_scale, comoving_bin = comoving_bin, one_d_fit = one_d_fit, init_guess = init_guess, bounds = bounds, n_grid_samples = n_grid_samples, method = 'mcmc')
            curr_time = time.time()
            if i % 100 == 100 - 1:
                n_past = i + 1
                n_to_go = len(spherical_points) - n_past
                delta_t = curr_time - start_time
                print ('!!! We are ' + str(can.round_to_n(i / len(spherical_points) * 100, 5)) + '% done.  Approximately ' + str(can.round_to_n(delta_t * n_to_go / n_past, 5)) + 's to go.')
                prev_time = curr_time
        return [spherical_points, fits]

    def saveOneDFits(self, fitted_points, fits, save_file, header = None):
        """
        Take a series of fits run on a comoving grid in 3D and save those
            fits to a plaintext file.
        """
        good_fit_indeces = [i for i in range(len(fits)) if fits[i] != None]
        cols_to_save = [ [fitted_points[i][0] for i in good_fit_indeces],
                         [fitted_points[i][1] for i in good_fit_indeces],
                         [fitted_points[i][2] for i in good_fit_indeces]]
        n_good_fits = len(cols_to_save[0])
        #print ('fits = ' + str(fits))
        print ('n_good_fits = ' + str(n_good_fits))
        if n_good_fits > 0:
            n_fit_params = len(fits[good_fit_indeces[0]][0])
        else:
            n_fit_params = 0
        cols_to_save = cols_to_save + [ [fits[i][0][j] for i in good_fit_indeces ] for j in range(n_fit_params ) ]
        cols_to_save = cols_to_save + [ [fits[i][1] for i in good_fit_indeces],
                                        [fits[i][2] for i in good_fit_indeces],
                                        [fits[i][3] for i in good_fit_indeces] ]
        results_dir = self.dir_archive.getFieldFitsDir()
        can.saveListsToColumns(cols_to_save, save_file, results_dir, header = header, sep = ', ')
        return 1

    def makePlotOfOneDFits(self, fitted_points_spher_coords, fits, fig_size_unit = 0.5, n_fits_x = 6, n_fits_y = 2,
                           legend_fontsize = 7, ticksize = 8, labelsize = 10, comoving_bin = None, angular_scale = None, show = 0,
                           save_plot_prefix = '', threeD_plot_suffix = '_FitsIn3Sky.pdf', gridspec_density_plot_suffix = '_BestDensityFits.pdf', vel_plot_density_plot_suffix = '_BestDensityVelocities.pdf',
                           hist_of_fits = '_HistOfFitImprovements.pdf', vel_lim_min = 0.0005):
        """
        Take a series of fits run on a comoving grid in 3D and show the best
            n_good_fits_to_show of those fits.  In the center of these plots,
            we show the distribution of SNe on the sky.
        """
        if comoving_bin == None:
            if angular_scale == None:
                angular_scale = self.binning_angular_scale
        plot_dir = self.dir_archive.getPlotDirectory()
        good_indeces  = [i for i in range(len(fitted_points_spher_coords)) if fits[i] != None]
        fitted_points = [fitted_points_spher_coords[i] for i in good_indeces]
        fitted_param_vals = [fits[i][0] for i in good_indeces]
        fitted_r_chi_squares = [fits[i][1] for i in good_indeces]
        fitted_null_chi_squares = [fits[i][2] for i in good_indeces]
        fitted_improvements = [fits[i][3] for i in good_indeces]
        """
        ordered_points, ordered_param_vals, ordered_r_chi_squares, ordered_null_chi_squares, ordered_improvements = can.safeSortOneListByAnother(fitted_improvements, [fitted_points, fitted_param_vals, fitted_r_chi_squares, fitted_null_chi_squares, fitted_improvements])
        ordered_points = can.niceReverse(ordered_points)
        ordered_param_vals = can.niceReverse(ordered_param_vals)
        ordered_r_chi_squares = can.niceReverse( ordered_r_chi_squares)
        ordered_null_chi_squares = can.niceReverse( ordered_null_chi_squares)
        ordered_improvements = can.niceReverse(ordered_improvements)
        """
        ordered_points, ordered_param_vals, ordered_r_chi_squares, ordered_null_chi_squares, ordered_improvements = can.safeSortOneListByAnother(fitted_r_chi_squares, [fitted_points, fitted_param_vals, fitted_r_chi_squares, fitted_null_chi_squares, fitted_improvements])

        nearby_sn_indeces = [self.getSNWithinComovingDistance(ordered_points[i][0], ordered_points[i][1:], (self.getComovingCrossSectionOfAngularScaleAtRedshift(ordered_points[i][0], angular_scale) if comoving_bin == None else comoving_bin)) for i in range(len(ordered_points))]
        deg_to_rad = self.astro_arch.getDegToRad()
        cart_points = [point[0] * np.array([np.cos(point[1] * deg_to_rad ) * np.sin((90.0 - point[2]) * deg_to_rad ), np.sin(point[1] * deg_to_rad ) * np.sin((90.0 - point[2]) * deg_to_rad ), np.cos((90.0 - point[2]) * deg_to_rad ) ]) for point in ordered_points]
        xs_cart, ys_cart, zs_cart = [[point[0] for point in cart_points], [point[1] for point in cart_points], [point[2] for point in cart_points] ]
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.scatter(xs_cart, ys_cart, zs_cart, c = np.log(ordered_improvements))
        plt.tight_layout()
        plt.savefig(plot_dir + save_plot_prefix + threeD_plot_suffix)
        plt.close('all')

        max_n_good_fits_to_show = n_fits_x * (n_fits_y - 1) + 2
        figsize = ( n_fits_x * fig_size_unit , n_fits_y * 3 * fig_size_unit)
        fig = plt.figure(constrained_layout=True, figsize = figsize)
        gs = fig.add_gridspec(n_fits_y * 2 * 3 + 1, n_fits_x * 2)
        gs_indeces_for_single_fits = [ [ [[(i // n_fits_x) * 4, (i // n_fits_x) * 4 + 1], [(i % n_fits_x) * 2, (i % n_fits_x) * 2 + 1]],
                                         [[(i // n_fits_x) * 4 + 2, (i // n_fits_x) * 4 + 2 + 1], [(i % n_fits_x) * 2, (i % n_fits_x) * 2 + 1]],
                                         [[(i // n_fits_x) * 4 + 4, (i // n_fits_x) * 4 + 4 + 1], [(i % n_fits_x) * 2, (i % n_fits_x) * 2 + 1]] ]
                                        for i in range(max_n_good_fits_to_show - 2)]
        gs_indeces_for_single_fits = gs_indeces_for_single_fits + [ [ [[((max_n_good_fits_to_show-2) // n_fits_x) * 4 + 2, ((max_n_good_fits_to_show-2) // n_fits_x) * 4 + 2 + 1] ,[0, 1]],
                                                                      [[((max_n_good_fits_to_show-2) // n_fits_x) * 4 + 2 + 2, ((max_n_good_fits_to_show-2) // n_fits_x) * 4 + 2 + 2 + 1] ,[0, 1]],
                                                                      [[((max_n_good_fits_to_show-2) // n_fits_x) * 4 + 4 + 2, ((max_n_good_fits_to_show-2) // n_fits_x) * 4 + 4 + 2 + 1] ,[0, 1]]],
                                                                    [ [[((max_n_good_fits_to_show-2) // n_fits_x) * 4 + 2, ((max_n_good_fits_to_show-2) // n_fits_x) * 4 + 2 + 1] ,[ (n_fits_x - 1) * 2, (n_fits_x - 1) * 2 + 1]],
                                                                      [[((max_n_good_fits_to_show-2) // n_fits_x) * 4 + 2 + 2, ((max_n_good_fits_to_show-2) // n_fits_x) * 4 + 2 + 2 + 1] ,[(n_fits_x - 1) * 2, (n_fits_x - 1) * 2 + 1]],
                                                                      [[((max_n_good_fits_to_show-2) // n_fits_x) * 4 + 4 + 2, ((max_n_good_fits_to_show-2) // n_fits_x) * 4 + 4 + 2 + 1] ,[(n_fits_x - 1) * 2, (n_fits_x - 1) * 2 + 1]] ] ]

        #f, axarr = plt.subplots(2, max_n_good_fits_to_show, figsize = )
        #axarr[0,0].set_ylabel(r'$\Delta \mu$ (mag)')
        #axarr[1,0].set_ylabel(r'Dec (deg)')
        #max_n_good_fits_to_show = 7
        #print ('gs_indeces_for_single_fits = '  +str(gs_indeces_for_single_fits))
        for i in range(min(len(nearby_sn_indeces), max_n_good_fits_to_show, len(gs_indeces_for_single_fits))):
            gs_indeces = gs_indeces_for_single_fits[i]
            #print ('gs_indeces = ' + str(gs_indeces))
            #print ('ordered_param_vals = ' + str(ordered_param_vals))

            fitting_coord = ordered_points[i]
            print ('fitting_coord = ' + str(fitting_coord))
            fitted_params = ordered_param_vals[i]
            improvement = ordered_improvements[i]
            r_chi_sqr = ordered_r_chi_squares[i]
            null_r_chi_sqr = ordered_null_chi_squares[i]
            #fitted_mass = self.getHaloMassFromVariedParam(fitted_mass_power)
            #fitted_radius = self.getRadiusFromVariedParam(fitted_radius_power)
            sn_indeces_to_plot = nearby_sn_indeces[i]
            zs_to_plot = [self.all_zs[index] for index in sn_indeces_to_plot]
            muResids_to_plot = [self.all_muResids[index] for index in sn_indeces_to_plot]
            mu_errs_to_plot = [self.all_mu_errs[index] for index in sn_indeces_to_plot]
            RAs_to_plot = [self.all_RAs[index] for index in sn_indeces_to_plot]
            Decs_to_plot = [self.all_Decs[index] for index in sn_indeces_to_plot]
            surveys_to_plot = [self.all_surveys[index] for index in sn_indeces_to_plot]
            colors_to_plot = [self.survey_to_color_dict[survey] for survey in surveys_to_plot]
            #RA_offset_to_center = (180.0 - fitting_coord[1])
            #RAs_for_resids = self.centerRAs(RAs_to_plot, RA_offset_to_center)
            #central_RA, central_Dec = [self.centerRAs([fitting_coord[1]],  RA_offset_to_center )[0],  fitting_coord[2]]
            fitted_resids = self.muDiff_of_z_funct(zs_to_plot, RAs_to_plot, Decs_to_plot, zs_to_plot, RAs_to_plot, Decs_to_plot, muResids_to_plot, mu_errs_to_plot, [fitting_coord[1], fitting_coord[2]], [fitting_coord[0], 0.0, 0.0]  + fitted_params)
            null_resids = self.muDiff_of_z_funct(zs_to_plot, RAs_to_plot, Decs_to_plot, zs_to_plot, RAs_to_plot, Decs_to_plot, muResids_to_plot, mu_errs_to_plot, [fitting_coord[1], fitting_coord[2]], [fitting_coord[0], 0.0, 0.0] + self.param_init_guess)
            fitted_velocities =  self.getVelocityField(zs_to_plot, RAs_to_plot, Decs_to_plot, fitting_coord, 0.0, 0.0, fitted_params)[0]
            null_velocities = self.getVelocityField(zs_to_plot, RAs_to_plot, Decs_to_plot, fitting_coord, 0.0, 0.0, self.param_init_guess)[0]

            ax0 = fig.add_subplot(gs[gs_indeces[0][0][0]:gs_indeces[0][0][1] + 1, gs_indeces[0][1][0]:gs_indeces[0][1][1] + 1])
            ax1 = fig.add_subplot(gs[gs_indeces[1][0][0]:gs_indeces[1][0][1] + 1, gs_indeces[1][1][0]:gs_indeces[1][1][1] + 1])
            ax2 = fig.add_subplot(gs[gs_indeces[2][0][0]:gs_indeces[2][0][1] + 1, gs_indeces[2][1][0]:gs_indeces[2][1][1] + 1])
            ax0.scatter(zs_to_plot, muResids_to_plot, c = colors_to_plot, marker = 'o')
            ax0.scatter(zs_to_plot, fitted_resids, c = 'k', marker = 'x', alpha = 0.5)
            ax1.scatter(zs_to_plot, np.abs(fitted_velocities), c = 'k', marker = 'x')
            [ ax0.annotate(str(j + 1), (zs_to_plot[j], muResids_to_plot[j]), color = 'k', fontsize = ticksize, verticalalignment = 'center', horizontalalignment = 'center') for j in range(len(sn_indeces_to_plot)) ]
            ax0.errorbar(zs_to_plot, muResids_to_plot, yerr = mu_errs_to_plot, colors = colors_to_plot, fmt = 'none', ecolor = colors_to_plot)
            ax0.axvline(ordered_points[i][0], color = 'k', linestyle = '--')
            ax2.scatter(RAs_to_plot, Decs_to_plot, c = colors_to_plot, marker = 'o')
            [ ax2.annotate(str(j + 1), (RAs_to_plot[j], Decs_to_plot[j]), color = 'k', fontsize = ticksize, verticalalignment = 'center', horizontalalignment = 'center' ) for j in range(len(sn_indeces_to_plot)) ]
            ax2.scatter(fitting_coord[1], fitting_coord[2], c = 'k', marker = 'x', s = 10 )
            ax2.axvline(fitting_coord[1], linestyle = '--', color = 'k')
            ax2.axhline(fitting_coord[2], linestyle = '--', color = 'k')
            ax0.set_xlabel(r'$z$', fontsize = labelsize)
            ax1.set_xlabel(r'$z$', fontsize = labelsize)
            ax2.set_xlabel(r'RA (deg)', fontsize = labelsize)
            ax0.set_ylabel(r'$\Delta \mu$ (mag)', fontsize = labelsize)
            ax1.set_ylabel(r'$|v_{pec}| / c$', fontsize = labelsize)
            ax2.set_ylabel(r'Dec (deg)', fontsize = labelsize)
            ax0.tick_params(axis='both', labelsize= ticksize)
            ax1.tick_params(axis='both', labelsize= ticksize)
            ax2.tick_params(axis='both', labelsize= ticksize)
            default_ylim = ax0.get_ylim()
            ax0.set_ylim([default_ylim[0] - (default_ylim[1] - default_ylim[0]) * 0.1, default_ylim[1] + (default_ylim[1] - default_ylim[0]) * 0.1])
            default_vel_lim = ax1.get_ylim()
            if max([abs(default_vel_lim[0]), abs(default_vel_lim[1])]) < vel_lim_min:
                ax1.set_ylim([0.0, vel_lim_min])

            print ('[r_chi_sqr, null_r_chi_sqr] = ' + str([r_chi_sqr, null_r_chi_sqr] ))
            ax0.text(0.05, 0.97,
                                 '' # 'fit ' + str(can.round_to_n(improvement, 3)) + 'X better' + '\n'
                                + r'$\chi_\nu^2 / \chi_{\nu, 0}^2=$' + str(can.round_to_n(r_chi_sqr / null_r_chi_sqr, 6)) + '; '
                                #+ r'$\chi_\nu^2=$' + str(can.round_to_n(r_chi_sqr, 3)) + '; '
                                #+ r'$\chi_{\nu, 0}^2=$' + str(can.round_to_n(null_r_chi_sqr, 3))
                                , fontsize = labelsize, transform = ax0.transAxes, verticalalignment = 'top', horizontalalignment = 'left', color = 'grey')
            label_strs = self.param_info_dict['label_strs']
            unit_strs = self.param_info_dict['units']
            print ('[len(label_strs), len(fitted_params), len(unit_strs)] = ' + str([len(label_strs), len(fitted_params), len(unit_strs)]))
            val_text = '\n'.join([label_strs[i] + r'$=$' + str(can.round_to_n(fitted_params[i], 3)) + ' ' + unit_strs[i]  for i in range(len(fitted_params))])
            ax0.text(0.05, 0.03, val_text, fontsize = labelsize, transform = ax0.transAxes, verticalalignment = 'bottom', horizontalalignment = 'left', color = 'grey')
            #ax0.text(0.05, 0.03, r'$M=$' + str(can.round_to_n(fitted_mass, 3)) + r' T$M_{\odot}$' +  '\n' + r'$r_C=$' + str(can.round_to_n(fitted_radius, 3)) + r' Mpc', fontsize = labelsize, transform = ax0.transAxes, verticalalignment = 'bottom', horizontalalignment = 'left', color = 'grey')

        legend_gs_indeces = [(n_fits_y - 1) * 6, [2,-2]]
        sky_gs_indeces = [[(n_fits_y - 1) * 6 + 1, (n_fits_y - 1) * 6 + 7], [2,-2]]
        sky_ax = fig.add_subplot(gs[sky_gs_indeces[0][0]:sky_gs_indeces[0][1], sky_gs_indeces[1][0]:sky_gs_indeces[1][1]], projection="aitoff")
        plots_for_legend = self.makeSkyPlotOfSNe(sky_ax)
        legend_ax = fig.add_subplot(gs[legend_gs_indeces[0], legend_gs_indeces[1][0]:legend_gs_indeces[1][1]])
        legend_ax.legend(plots_for_legend, self.included_surveys, ncol = 7, fontsize = legend_fontsize)
        legend_ax.set_xticks([])
        legend_ax.set_yticks([])
        fig.tight_layout()
        fig.savefig(plot_dir + save_plot_prefix + gridspec_density_plot_suffix)
        if show: plt.show()
        plt.close('all')

        plt.hist(np.log10(ordered_improvements), bins = 20, edgecolor = 'k', color = 'white')
        plt.xlabel('Log10 of improvement over null fit')
        plt.ylabel(r'Log10($N$)')
        plt.tight_layout()
        plt.savefig(plot_dir + save_plot_prefix + hist_of_fits)

        return 1


    def __init__(self, data_set,
                      randomize_each_survey = 0, randomize_each_field = 0, randomize_all_sn = 0, randomized_sn_number = 0,
                      n_z_bins = 10,
                      zHD = 1, OmM = 0.3, OmL = 0.7, Om0 = 1.0, OmR = 0.0, H0 = 70.0, do_cosmic_fit = 0,
                      archive_healpix_sides = 32, annulus_inner_comoving_rad_scaling = 1.0, annulus_outer_comoving_rad_scaling = 2.0, gal_dens_weighting = 0.5,
                      full_sdss_gal_data_file = 'SDSS_PSSuperField3_SDSSGals_Allpz.csv', preloaded_sdss_gals = None, start_spherical_comoving_rad_frac = 0.03,
                      interp_z_params = [0.0, 100.0, 1001], z_range = [-0.1, 3.0], pull_extinctions = 0, surveys_to_include = ['all'], surveys_to_ignore = [],
                      n_mcmc_steps = 2000, archive_to_use = 'ps1md', min_v_scale_power = -5,
                      resid_from_grav = 0, resid_from_vel = 1, resid_profile_funct = 'exp_void',
                      plot_save_subdir = '',mcmc_extra_save_dir = 'mcmcs/',
                      fixed_comoving_scale_radius_power = 1.0, sn_data_type = 'pantheon_plus', low_log_prob_val = -100,
                      cutoff_radius_Mpc = 200, NFW_overdensity_param = 200, angular_scale_to_bin = 3.3 * 2  ,
                      hemisphere = 'both', angle_divs = 1, angle_slice = 1, rand_sn_file_prefix = 'Rand_SNe_Sequence_' #search_octant = 'all',
                      ):

        self.cutoff_radius = cutoff_radius_Mpc
        self.show_mcmc = (1 + randomize_all_sn) % 2
        self.n_mcmc_steps  = n_mcmc_steps
        self.overdensity_param = NFW_overdensity_param
        self.data_set = data_set
        self.astro_arch = AstronomicalParameterArchive()
        self.dir_archive = DirectoryArchive()
        self.plot_save_dir = self.dir_archive.getPlotDirectory() + plot_save_subdir
        self.rand_sn_file_prefix = rand_sn_file_prefix
        self.mcmc_plot_save_dir = self.plot_save_dir + mcmc_extra_save_dir
        self.deg_to_rad = self.astro_arch.getDegToRad()
        self.sn_data_type = sn_data_type
        self.n_z_bins = n_z_bins
        self.pull_extinctions = pull_extinctions
        self.z_range = z_range
        self.low_log_prob_val = low_log_prob_val
        self.zHD = zHD
        self.Om0 = Om0
        self.OmR = OmR
        self.OmM = OmM
        self.OmL = OmL
        self.H0 = H0
        self.surveys_to_ignore = surveys_to_ignore
        self.surveys_to_include = surveys_to_include
        self.randomize_all = randomize_all_sn
        self.randomize_by_survey = randomize_each_survey
        self.randomize_sn_number = randomized_sn_number
        self.hemisphere = 'both'
        self.angle_divs = 1
        self.angle_slice = 1
        self.hemisphere = 'both'
        self.angle_divs = 1
        self.angle_slice = 1
        print ('initializing SNe...')
        self.initializeSN(self.z_range, surveys_to_include, surveys_to_ignore, randomize_all = self.randomize_all, randomize_by_survey = self.randomize_by_survey, randomized_sn_number = randomized_sn_number)
        print ('SNe initialized.')
        if do_cosmic_fit:
            print ('Doing cosmic fit. ')
            self.OmM, self.OmL, self.H0 = self.doCosmicFit()
        else:
            print ('Not doing cosmic fit. ')
        #print ('self.mu_of_z = ' + str(self.mu_of_z))
        self.all_muResids = [sn['mu'] - self.mu_of_z(sn['z']) for sn in self.all_sns]
        self.interp_z_params = interp_z_params
        self.initialize_r_of_z_interp(interp_z_params)
        self.initialize_z_of_r_interp(interp_z_params)
        self.n_mcmc_steps = n_mcmc_steps
        #self.min_v_scale_power = min_v_scale_power #10^(%this number) = smallest velocity, as a fraction of c, that we consider # if number goes below this value, we flip it around and look at negative velocities
        self.fixed_comoving_scale_radius_power = fixed_comoving_scale_radius_power
        self.archive_healpix_sides =  archive_healpix_sides
        self.annulus_inner_comoving_rad_scaling = annulus_inner_comoving_rad_scaling
        self.annulus_outer_comoving_rad_scaling = annulus_outer_comoving_rad_scaling
        self.full_sdss_gal_data_file = full_sdss_gal_data_file
        self.preloaded_sdss_gals = preloaded_sdss_gals
        self.start_spherical_comoving_rad_frac = start_spherical_comoving_rad_frac
        self.gal_dens_weighting = gal_dens_weighting
        self.randomize = (randomize_each_survey)
        self.binning_angular_scale = angular_scale_to_bin

        #if randomize_each_survey:
        #    self.randomizeSNBySurvey()
        #elif randomize_all:
        #    self.randomizeAllSN()

        self.archive_to_use = archive_to_use
        print ('initializing fields...')
        self.initializeSNByField(archive_to_use, randomize = randomize_each_field)
        print ('fields initialized...')

        #if resid_profile_funct == 'redshift_from_grav':
        #    #For a redshift due to a gravitational well
        #    self.fit_funct, self.muDiff_of_z_funct, self.n_fit_params, self.bounds, self.mcmc_start_params, self.mcmc_step_sizes  = self.getFittingFunctG()
        #How do we calculate the total redshift: gravitational, doppler, or both
        if resid_from_grav and resid_from_vel:
            #!!! Both not written yet !!!
            self.muDiff_of_z_funct = self.ResidFunctFromGravAndVelocity
        elif resid_from_grav:
            #Gravitational
            self.muDiff_of_z_funct = self.redshiftDueToGravPotential
        else:
            #Peculiar velocity
            self.muDiff_of_z_funct = self.getMuDiffOfVelocityField

        self.resid_profile_funct = resid_profile_funct
        self.CosmicDensityProfile = cdp.CosmicDensityProfile(density_profile_type = self.resid_profile_funct)
        self.radial_mass_funct, model_free_params, self.bounds, self.param_conversion_funct, self.param_info_dict, self.param_init_guess = self.CosmicDensityProfile.getProfileFittingPieces()

        self.n_fit_params = model_free_params + 1 #An additional fit parameter, in the form of the subtracted away mean of SNe in a field
        self.n_null_fit_params = 1
        self.fit_params_by_field = {field_key:[] for field_key in self.fields.keys()}
        self.prob_fits_by_field = {field_key:{'sne_chi_sqr':0.0, 'sne_null_chi_sqr':0.0, 'sne_chi_sqr_prob':0.0, 'n_sig_gal_dens':[0.0, 0.0], 'n_sig_gal_prob':0.0, 'sne_overall_prob':0.0} for field_key in self.fields.keys()}
        print ('Starting to assign best fit curves...')
        self.assignBestFitCurves()
        print ('Done assigning best fit curves...')
        self.min_res = { }

        self.null_params = [0.5, 0.0, 0.1, 0.1]
        self.comoving_sne_seps = {}
        #Compute the comoving separations between supernovae so we can select which ones are worth including in overdensity fits.
        measure_seps = 0
        if measure_seps:
            print ('Computing comoving separations of supernovae...')
            for i in range(len(self.all_zs)):
                start_time = time.time()
                for j in range(i + 1, len(self.all_zs)):
                    comoving_sep = self.computeComovingSepsOfSNe([i,j])
                    self.comoving_sne_seps[(i,j)] = comoving_sep
                    self.comoving_sne_seps[(j,i)] = comoving_sep
                end_time = time.time()

def loadPickledPlotter(file_to_load):
    loaded_fielder = pickle.load(open(file_to_load, 'rb'))

    return loaded_fielder


if __name__ == "__main__":
    # $ python makePlotOfPS1MDFieldsClass.py 1 28 150 14 0.8 0 both 180 1
    # $ python makePlotOfPS1MDFieldsClass.py 1 28 150 14 0.8 1 both 180 1
    line_args = sys.argv[1:]
    #fitter_id = line_args[0]
    #field_ids = line_args[1]
    fit_id = line_args[0]
    fit_id = int(fit_id)
    fit_comoving_density = int(line_args[1]) #50
    comoving_bin = line_args[2]
    if comoving_bin == 'None':
        comoving_bin = None
    comoving_bin = int(comoving_bin)
    min_n_sn = int(line_args[3]) #15
    max_redshift = float(line_args[4])
    randomize_all_sn = int(line_args[5])
    #search_octant = line_args[6]
    #if search_octant.isdigit():
    #    search_octant = int(search_octant)
    hemisphere = line_args[6]
    if hemisphere.isdigit():
        hemisphere = int(hemisphere)
    angle_divs = int(line_args[7])
    angle_slice = int(line_args[8])
    print ('[fit_id, fit_comoving_density, comoving_bin, min_n_sn, max_redshift, randomize_all_sn, hemisphere, angle_divs, angle_slice] = ' + str([fit_id, fit_comoving_density, comoving_bin, min_n_sn, max_redshift, randomize_all_sn, hemisphere, angle_divs, angle_slice]))

    #To calculate this, do:
    # field_plotter = PanStarsFieldManager(1, full_sdss_gal_data_file = 'SDSS_fullCoverage_SDSSGals_pzAll.csv', preloaded_sdss_gals = fulldata_fastRead, gal_dens_weighting = gal_dens_weighting, z_range = z_range, sn_data_type = sn_data_type, zHD = zHD, cutoff_radius_Mpc = init_comoving_bin_guess, NFW_overdensity_param = overdensity_param, resid_from_grav = 0, resid_from_vel = 1, resid_profile_funct = resid_profile_funct, randomize_all_sn = randomize_all_sn, n_mcmc_steps = n_mcmc_steps, search_octant = 'all', angular_scale_to_bin = angular_scale_to_bin )
    # smallest_cluster_by_n_sn[min_n_sn] = field_plotter.determineMinimumSpacingComovingBinningByNSNe(min_n_sn)
    #
    smallest_cluster_by_n_sn = {2: 0.0, 3: 0.0, 4:0.0, 5: 0.4706012544945625, 6: 1.291879585616636, 7: 5.273206273547085, 8: 5.873918626542572, 9:8.088450520567488, 10: 8.871590795028839, 11: 8.871590795028839, 12: 13.384796636581399, 13:15.474809950037185, 14:15.47826263944724, 15:17.073129458455856, 16:17.946806338740117}


    init_comoving_bin_guess = 200
    angular_scale_PS1MD = 3.3   #The angular scale of PS1 is 3.3 degree in diameter (see Abstract of: https://iopscience.iop.org/article/10.1088/0004-637X/745/1/42/pdf).
    angular_scale_to_bin = angular_scale_PS1MD * 2
    z_range = [-0.1, 3.0]
    min_redshift = 0.01
    n_mcmc_steps = 2000
    overdensity_param = 200 #Delta, the fraction of background mass density that the average halo mass density must be
    n_good_fits_to_show = 5
    resid_profile_funct = 'exp_void'
    #randomize_all_sn = 0
    #print ('field_ids = ' + str(field_ids))
    #field_ids[1:-1].split(',')
    #field_ids = [int(id) for id in field_ids[1:-1].split(',')]
    #print('field_ids = ' + str(field_ids))
    do_randomization_by_field = 0
    do_randomization_by_survey = 0
    gal_dens_weighting = 0.0 #can be 0.5
    save_plot = 1
    #init_guess = [0.0, 2]
    #bounds = [[-500.0, 500.0], [0.0, 3.0]]
    n_grid_samples = [101, 31]
    sn_data_type = 'pantheon_plus' #'pantheon_plus' #real - for old Pantheon data
    zHD = 1 # Was 1 for old Pantheon data
    #fit_params_file =  ('RandField' if do_randomization_by_field else 'RandSurvey' if do_randomization_by_survey else 'True') + 'PS1MDFieldFitter_field' + '_'.join([str(elem) for elem in field_ids])  + '_GalWeight0p' + str(int(10 * gal_dens_weighting)) + '_N' + str(fitter_id) + '.csv'
    #plot_file_name = ('RandField' if do_randomization_by_field else 'RandSurvey' if do_randomization_by_survey else 'True') + 'PS1MDFieldFitter_' + sn_data_type + '_field' + '_'.join([str(elem) for elem in field_ids]) + '_GalWeight0p' + str(int(10 * gal_dens_weighting)) + 'z_range' + str(z_range[0]) + '_' + str(z_range[1]) +'_N' + str(fitter_id) + '.pdf'
    sdssdir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/SNIsotropyProject/SDSSGalaxies/'

    read_in_galaxy_data = 0
    if read_in_galaxy_data:
        fulldata_fastRead = can.readInColumnsToList('SDSS_fullCoverage_SDSSGals_pzAll.csv', sdssdir, delimiter = ',', n_ignore = 2, all_np_readable = 1)
    else:
        fulldata_fastRead = None

    #print('fit_id = ' + str(fit_id))
    #print('fit_params_file =  ' + str(fit_params_file))
    field_plotter = PanStarsFieldManager(1, full_sdss_gal_data_file = 'SDSS_fullCoverage_SDSSGals_pzAll.csv', preloaded_sdss_gals = fulldata_fastRead, gal_dens_weighting = gal_dens_weighting, z_range = z_range, sn_data_type = sn_data_type, zHD = zHD, cutoff_radius_Mpc = init_comoving_bin_guess, NFW_overdensity_param = overdensity_param, resid_from_grav = 0, resid_from_vel = 1, resid_profile_funct = resid_profile_funct, randomize_all_sn = randomize_all_sn, randomized_sn_number = (fit_id if randomize_all_sn else 0), n_mcmc_steps = n_mcmc_steps, hemisphere = hemisphere, angle_divs = angle_divs, angle_slice = angle_slice, angular_scale_to_bin = angular_scale_to_bin )# , surveys_to_include = ['PS1MD' ,'SDSS', 'SNLS'])
    #comoving_bin = field_plotter.getComovingCrossSectionOfAngularScaleAtRedshift(max_redshift, angular_scale_to_bin )
    #comoving_bin = None
    print ('comoving_bin = ' + str(comoving_bin) + ' Mpc')
    field_plotter.cutoff_radius = comoving_bin
    #min_redshift = field_plotter.z_of_r_interp(comoving_bin)
    print ('[min_redshift, max_redshift] = ' + str([min_redshift, max_redshift] ))
    #field_fitter =  PanStarsFieldManager(1, randomize_each_field = do_randomization_by_field, randomize_each_survey = do_randomization_by_survey, surveys_to_include = ['PS1MD' ,'SDSS', 'SNLS'])
    #print ('field_plotter.mcmc_start_params = ' + str(field_plotter.mcmc_start_params))
    spherical_points, fits = field_plotter.doFitsOnGrid( fit_comoving_density, min_n_sn, z_range = [min_redshift, max_redshift], init_guess = None, bounds = None, n_grid_samples = n_grid_samples, hemisphere = hemisphere, angle_divs = angle_divs, angle_slice = angle_slice, comoving_bin = comoving_bin )
    save_prefix = 'OverdensityFit' + str(fit_id) + '_MinNSN_' + str(min_n_sn) + '_GridDens_' + str(fit_comoving_density) + '_BinSize_' + str(comoving_bin) + '_Hemisphere_' + str(hemisphere) + '_RArange_' + str(angle_slice) + 'of' + str(angle_divs) + '_Rand' + str(randomize_all_sn)
    #field_plotter.makePlotOfOneDFits(spherical_points, fits, fig_size_unit = 2.5, save_plot_prefix = save_prefix)


    #Finally, do a big fit centered at the best fit location, which includes all supernovae

    if len([fit for fit in fits if fit != None]) > 0:
        final_n_mcmc_steps = 10000
        #final_n_mcmc_steps = 4000
        field_plotter.n_mcmc_steps = final_n_mcmc_steps
        best_fit_index = np.argmin([fit[1] / fit[2] if fit != None else np.inf for fit in fits ])
        best_fit_point, best_fit_params = [spherical_points[best_fit_index], fits[best_fit_index][0]]
        print ('best_fit_point = ' + str(best_fit_point))
        print ('best_fit_params = ' + str(best_fit_params))
        final_fit, final_r_chi_square, final_null_r_chi_square, final_null_prob = field_plotter.doFullFitAroundPoint([best_fit_point[0]] + [0.0, 0.0] + best_fit_params, best_fit_point[1:], comoving_bin, min_n_sn, method = 'mcmc', mcmc_save_suffix = save_prefix )
        print ('final_fit = ' + str(final_fit))
        final_params = final_fit['x']

        spherical_points = spherical_points + [[final_params[0], best_fit_point[1] + final_params[1] / np.cos(np.deg2rad(best_fit_point[2] + final_params[2])), best_fit_point[2] + final_params[2]]]
        fits = fits + [[final_params[3:]] + [final_r_chi_square, final_null_r_chi_square, final_null_prob]]

    field_plotter.saveOneDFits(spherical_points, fits, save_file= save_prefix + '_fits.txt', header = 'Grid Density ' + str(fit_comoving_density) + ' Mpc' + '\n'                                                                                        + 'Min Nearby SN ' + str(min_n_sn) + '\n'
                                                                                                   + 'Comoving bin ' + str(comoving_bin) + '\n'
                                                                                                    + 'z, RA (deg), Dec (deg), delta_0, r_0 (Mpc), z_c, Delta phi (deg), Delta theta (deg), Full R Chi^2, Full Null R Chi^2, Prob Improvement'  )
    field_plotter.makePlotOfOneDFits(spherical_points, fits, fig_size_unit = 2.5, save_plot_prefix = save_prefix, comoving_bin = comoving_bin)

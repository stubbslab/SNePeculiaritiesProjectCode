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
from FitStorer import FitStorer
from SDSSArchive import SDSSArchive
from checkIfInCyclicRange import checkIfInCyclicRange
from AstronomicalParameterArchive import AstronomicalParameterArchive
from matplotlib.patches import Rectangle
import loadSN as lsn


class PanStarsFieldManager:

def getErroneousRedshiftPlot(zg_of_z_funct, z_space,
                             r_of_z_interp = None,
                             zHD = 1, z_include_range = [0.0, np.inf], surveys_of_interest = ['all'], surveys_to_excise = [''],
                             pull_extinctions = 0, OmM = 0.3, OmL = 0.7, H0 = 70.0, Om0 = 1.0, OmR = 0.0, astro_arch = None,
                             ):

    if astro_arch is None:
        astro_arch = AstronomicalParameterArchive()
    c = astro_arch.getc() # in km/s

    #Note H0 is given in km/s/Mpc.  c is in km/s.  So c/H0 gives this ratio in Mpc.
    if r_of_z_interp is None:
        dLInMpc = lambda zf, zg: c / H0 * (1 + zf) * scipy.integrate.quad(lambda z_int: 1/np.sqrt((1 + z_int) ** 3 * OmM + (1 + z_int) ** 0 * OmL + (1 + z_int) ** 4 * OmR), 0, (zf - zg)/(zg + 1))[0]
    else:
        #dLInMpc = lambda zf, zg: c / H0 * (1 + zf) * r_of_z_interp((zf - zg)/(zg + 1))
        dLInMpc = lambda zf, zg: (1 + zf) * r_of_z_interp((zf - zg)/(zg + 1))
    mu =  lambda zf, zg: 5 * np.log10(dLInMpc(zf, zg)) + 25

    #z_space = np.linspace(min_z, max_z, n_zs)
    #plt.plot(z_space, [mu(z, 0.0) for z in z_space], c = 'k')
    #plt.plot(z_space, [mu(z, zg_of_z(z, 0.0, 0.5, 0.05)) - mu(z, 0.0) for z in z_space], c = 'k')
    #plt.plot(z_space, [mu(z, zg_of_z(z, 0.1, 0.5, 0.05)) - mu(z, 0.0) for z in z_space], c = 'r')
    #plt.scatter(all_zs, all_muDiffs, marker = '.')
    #plt.errorbar(all_zs, all_muDiffs, yerr = all_muErrs, fmt = 'none')
    #plt.show( )

    #print ('[min([(z - zg_of_z_funct(z))/(zg_of_z_funct(z) + 1) for z in z_space]), max([(z - zg_of_z_funct(z))/(zg_of_z_funct(z) + 1) for z in z_space])] = ' + str([min([(z - zg_of_z_funct(z))/(zg_of_z_funct(z) + 1) for z in z_space]), max([(z - zg_of_z_funct(z))/(zg_of_z_funct(z) + 1) for z in z_space])]))
    #print ('[[z, zg_of_z_funct(z)] for z in z_space[0:5]] = ' + str([[z, zg_of_z_funct(z)] for z in z_space[0:5]]))
    #print ('[[z, dLInMpc(z, zg_of_z_funct(z)), dLInMpc(z, 0.0)] for z in z_space[0:5]] = ' + str([[z, dLInMpc(z, zg_of_z_funct(z)), dLInMpc(z, 0.0)] for z in z_space[0:5]]))
    #print ('[[z, mu(z, zg_of_z_funct(z)), mu(z, 0.0)] for z in z_space[0:5]] = ' + str([[z, mu(z, zg_of_z_funct(z)), mu(z, 0.0)] for z in z_space[0:5]]))
    #print ('[min(z_space), max(z_space) = ' + str([min(z_space), max(z_space[-1]]))
    vals = [mu(z, zg_of_z_funct(z)) - mu(z, 0.0) for z in z_space]

    return vals

def plot_mwd(RA, Dec, ax, org=0,title='Aitoff Sky', projection='aitoff'):
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

def plotPantheonSNOnSky(data_set,
                        projection = 'aitoff', save = 1, show = 0, surveys_to_display = 'all', z_range = [-0.1, 3.0],
                        pull_extinctions = 0, figsize = [10.0, 8.0], zHD = 1):

    dir_archive = DirectoryArchive()
    plot_dir = dir_archive.getPlotDirectory()
    astro_arch = AstronomicalParameterArchive()
    deg_to_rad = astro_arch.getDegToRad()

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


def printingFitFunct(fit_funct, zs, muDiffs, muErrs, params, print_params = 0):
    fit_vals = fit_funct(zs, *params)
    diffs_from_fit = (np.array(fit_vals) - np.array(muDiffs))
    weighted_mean = np.sum(diffs_from_fit * np.array(muErrs) ** (-2.0)) / np.sum( np.array(muErrs) ** (-2.0))
    chi_sqr = np.sum((diffs_from_fit - weighted_mean) ** 2.0 / np.array(muErrs) ** 2.0)
    if print_params: print ('params.tolist() + [chi_sqr] = ' +str(params.tolist() + [chi_sqr]) )
    return fit_vals, chi_sqr, weighted_mean


def doMinimization(fit_funct, zs, resids, errs, fit_params_A, fit_params_B, fit_params_C, dof, fit_params_D = None, fit_params_E = None, bounds= None):

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
                        this_fit_resids, this_fit_chisqr, this_fit_weighted_mean = printingFitFunct(fit_funct, zs, resids, errs, np.array(this_fit))
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
        final_minimization = scipy.optimize.minimize(lambda params: printingFitFunct(fit_funct, zs, resids, errs, params, print_params = 1)[1], brute_minimization, bounds = bounds)#, method='Nelder-Mead')
        #fit_res = optimize.curve_fit(lambda zs, A, mu: fit_funct(zs, A, mu, sig, 0.2, 0.0), all_zs, all_resids, p0 = [0.0, max(all_zs) / 2.0], sigma=all_errs, maxfev = maxfev)
        this_funct_fit = final_minimization['x']
    except RuntimeError:
        print("Curve_fit failed!.  Plotting initial guess. ")
        this_funct_fit = np.array(brute_minimization)
    final_chisqr = printingFitFunct(fit_funct, zs, resids, errs, np.array(this_funct_fit))[1]
    print ('[this_funct_fit, final_chisqr] = ' + str([this_funct_fit, final_chisqr]))
    final_rchisqr = final_chisqr / dof

    return [this_funct_fit, final_chisqr]


def d_H_d_z (single_z):


def redshift_of_v_funct(beta):
    redshift = np.sign(beta) * (np.sqrt((1+abs(beta))/(1-abs(beta))) - 1)
    return redshift

def d_z_D_growth_over_D(single_z):


def computeVelocityFieldFromNFWHalo(zs_to_calc, central_z, critical_mass, concentration, impact_param,
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
    #delta_portion_funct = lambda single_z, central_z, delta_s, rs, c, b_incidence: rs ** 1.0 * (rs / r_well_from_geometry(r_of_z_interp(single_z), r_of_z_interp(central_z), b_incidence * rs)) ** 2.0 * delta_s * (int_of_rsqr_delta_nfw( min(c, r_well_from_geometry(r_of_z_interp(single_z), r_of_z_interp(central_z), b_incidence * rs) / rs), c) ) if (single_z != central_z or abs(b_incidence) > 0.0) else 0.0
    delta_portion_funct = lambda single_z, central_z, delta_s, rs, c, b_incidence: rs ** 1.0 * (rs / r_well_from_geometry(r_of_z_interp(single_z), r_of_z_interp(central_z), b_incidence * rs)) ** 2.0 * delta_s * (int_of_rsqr_delta_nfw( min(c, r_well_from_geometry(r_of_z_interp(single_z), r_of_z_interp(central_z), b_incidence * rs) / rs), c) ) if (single_z != central_z or abs(b_incidence) > 0.0) else 0.0





def makePlotOfPS1MDFields(data_set,
                          n_bins = 20, bin_scheme = 'bin_size', save = 1, show = 0, surveys_to_display = 'all', surveys_to_ignore = [], fit_information = {'funct':'none'}, show_fit_label = 1, z_range = [-0.1, 3.0],
                          res_range = [-0.8, 0.8], binned_res_range = [-0.4, 0.4], archive_to_use = 'PS1MD' , pull_extinctions = 0, figsize = [16.0, 8.0], separate_fields_by_plot = 0,
                          n_randomizations = 0, n_z_bins = 10, zHD = 0, OmM = 0.3, OmL = 0.7, Om0 = 1.0, OmR = 0.0, H0 = 70.0,
                          interp_z_params = [0.0, 100.0, 1001]):

    master_start = time.time()
    dir_archive = DirectoryArchive()
    plot_dir = dir_archive.getPlotDirectory()
    astro_arch = AstronomicalParameterArchive()
    deg_to_rad = astro_arch.getDegToRad()
    pc_to_m = astro_arch.getParsecToM()
    Mpc_to_km = pc_to_m * 10.0 ** 6.0 * 10.0 ** -3.0

    speedol = astro_arch.getc()

    PSArch = PanStars1Archive()
    SDSSArch = SDSSArchive()

    if archive_to_use.lower() == 'sdss':
        archive = SDSSArch
    else:
        archive = PSArch
    fields = archive.fields

    all_sns = loadSN(1, ['all'], pull_extinctions = pull_extinctions, zHD = zHD, OmM = OmM, OmL = OmL, Om0 = Om0, OmR = OmR, H0 = H0)
    all_sns = [sn for sn in all_sns if not(sn['survey'] in surveys_to_ignore)]
    H_of_z = lambda z_for_H: np.sqrt((1 + z_for_H) ** 3 * OmM + (1 + z_for_H) ** 0 * OmL + (1 + z_for_H) ** 4 * OmR)
    d_H_d_z = lambda z_for_H: 0.5 * 1.0 / np.sqrt((1 + z_for_H) ** 3 * OmM + (1 + z_for_H) ** 0 * OmL + (1 + z_for_H) ** 4 * OmR) * (3.0 * (1 + z_for_H) ** 2.0 * OmM + 4.0 * (1 + z_for_H) ** 3 * OmR)
    r_of_z = lambda z_meas: (speedol) / H0 * scipy.integrate.quad( lambda z_int: 1.0 / H_of_z(z_int), 0, z_meas )[0]
    interp_zs = np.linspace(interp_z_params[0], interp_z_params[1], interp_z_params[2]).tolist() + [1000.0]
    r_of_z_interp = scipy.interpolate.interp1d([-1000.0] + interp_zs, [0.0] + [r_of_z(z) for z in interp_zs], fill_value = (0.0, r_of_z(1000)), bounds_error = False)
    all_sns = [sn for sn in all_sns if (sn['z'] >= z_range[0] and sn['z'] <= z_range[1]) ]
    all_zs = [sn['z'] for sn in all_sns]
    min_z, max_z = [np.min(all_zs), np.max(all_zs)]
    n_zs_to_display_fit = 1001
    zs_to_plot = np.linspace(min_z, max_z, n_zs_to_display_fit)
    if surveys_to_display in ['all','All','ALL']:
        surveys_to_display = [sn['survey'] for sn in all_sns]
        surveys_to_display = np.unique(surveys_to_display).tolist()
    print ('surveys_to_display = ' + str(surveys_to_display ) )
    sn_by_survey = [[sn for sn in all_sns if sn['survey'] == survey] for survey in surveys_to_display if len([sn for sn in all_sns if sn['survey'] == survey]) > 0]
    #print ('sn_by_survey = '  +str(sn_by_survey))
    colors = [sn_set[0]['color'] for sn_set in sn_by_survey]

    sn_in_survey_by_field = []

    for sn_set in sn_by_survey:
        sn_by_field = {}
        for key in fields:
            field = fields[key]
            #sn_by_field[key] = [sn for sn in sn_set if (sn['RA']>field[0] and sn['RA'] < field[1] and sn['Dec'] > field[2] and sn['Dec'] < field[3] )]
            sn_by_field[key] = [sn for sn in sn_set if (sn['RA']>field[0] and sn['RA'] < field[1] and sn['Dec'] > field[2] and sn['Dec'] < field[3] )]
        sn_in_survey_by_field = sn_in_survey_by_field + [sn_by_field]

    n_fields = len(fields)
    fig_side = 0
    while fig_side * 4 - 4 < n_fields:
        fig_side = fig_side + 1

    fig = plt.figure(constrained_layout=True, figsize = figsize)
    if not(separate_fields_by_plot):
        gs = fig.add_gridspec(fig_side, fig_side)
        plot_indeces = [[fig_side - 1, i] for i in range(fig_side)] + [[0, i] for i in range(fig_side)] + cant.flattenListOfLists( [[[i, 0], [i, fig_side - 1]] for i in range(1, fig_side-1)] )
        for plot_index in plot_indeces:
            fig.add_subplot(gs[plot_index[0], plot_index[1]])
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
        ax.legend(plots_for_legend, surveys_to_display, ncol = 3, fontsize = 8)

    field_nums_to_display = list(range(n_fields))
    #field_nums_to_display = [0,1,2]
    for field_num in field_nums_to_display:
    #for field_num in range(1):
        field = fields[field_num]
        field_str = 'f' + str(field_num)
        if not(separate_fields_by_plot):
            plot_index = plot_indeces[field_num]
            ax = fig.add_subplot(gs[plot_index[0], plot_index[1]])
            rect_coords1 =[ -(field[1] if field[1] < 180.0 else field[1] - 360), field[2]]
            rect_coords2 =[ -(field[0] if field[0] < 180.0 else field[0] - 360) - rect_coords1[0] , field[3] - field[2]]
            field_rect = Rectangle([rect_coord * deg_to_rad for rect_coord in rect_coords1], *[rect_coord * deg_to_rad for rect_coord in rect_coords2], edgecolor='black', facecolor='black', alpha = 0.25)
            sky_plot.add_patch(field_rect)
            sky_plot.text((rect_coords1[0] + rect_coords2[0]) * deg_to_rad, (rect_coords1[1] + rect_coords2[1]) * deg_to_rad, field_str)
        else:
            fig, ax = plt.subplots(1, sharex = True, figsize = figsize)
        print ('Displaying field number ' + str(field_num))
        survey_plots = []
        fitted_plots = []
        fit_position_index = 0
        min_z = np.inf
        max_z = -np.inf
        min_mu = np.inf
        max_mu = -np.inf
        all_zs = []
        all_dls = []
        all_resids = []
        all_errs = []
        fitters = [ [FitStorer(fit_information) for survey in surveys_to_display ] for field in fields]
        for i in range(len(surveys_to_display)):
            sns_in_field_in_survey = sn_in_survey_by_field[i][field_num]
            zs = [sn['z'] for sn in sns_in_field_in_survey]
            dls = [sn['dl'] for sn in sns_in_field_in_survey]
            rs = [dls[j] / (1.0 + zs[j]) for j in range(len(zs))]
            all_zs = all_zs + zs
            muDiffs = [sn['muDiff'] - sn['muDiffWMean'] for sn in sns_in_field_in_survey]
            all_resids = all_resids + muDiffs
            muErrs = [sn['muErr'] for sn in sn_in_survey_by_field[i][field_num]]
            all_errs = all_errs + muErrs
            #if len(zs) > 0:
            #    min_z = min(min_z, min(zs))
            #    max_z = max(max_z, max(zs))
            #if len(muDiffs):
            #    min_mu = min(min_mu, np.min(np.array(muDiffs) + np.array(muErrs)))
            #    max_mu = max(max_mu, np.max(np.array(muDiffs) + np.array(muErrs)))
            color = colors[i]
            survey = surveys_to_display[i]

            if len(zs) > 0:
                survey_plots = survey_plots + [ax.scatter(zs, muDiffs, c = color) ]
                ax.errorbar(zs, muDiffs, yerr = muErrs, ecolor = color, fmt = 'none')
            else:
                survey_plots = survey_plots + [ax.scatter([],[])]


            #if len(sns_in_field_in_survey) > 0:
            #    z_bin_centers, binned_Sn_data = binData(zs, muDiffs, y_errs = muErrs, n_bins = n_bins, bin_scheme = bin_scheme)

            #    axarr[1].scatter(z_bin_centers, binned_Sn_data[0], c = color)
            #    axarr[1].errorbar(z_bin_centers, binned_Sn_data[0], yerr = binned_Sn_data[1], fmt = None, ecolor = color)
            if fitters[field_num][i].getDimension() == 1 and len(zs) > 2:
                print ('Doing fit.')
                fitters[field_num][i].generateFit(zs, muDiffs, muErrs)
                #print ('For survey ' + surveys_to_display[i] + ', fit params are: ')
                #print (fitters[field_num][i].fit_params)
            #    #print [ [[ sn['z'] for sn in sns if sn['survey'] == surveys_to_display[i] ], fitters[i].getFitValues([ sn['z'] for sn in sns if sn['survey'] == surveys_to_display[i] ])] for i in range(len(surveys_to_display)) ]
                z_step = 0.001
                extra_points_for_fit = np.arange(z_range[0],z_range[1], z_step).tolist()
                fitted_plots = fitted_plots + [ ax.plot( sorted(zs + extra_points_for_fit ), fitters[field_num][i].getFitValues(sorted(zs + extra_points_for_fit)), c =color ) ]
                if show_fit_label:
                    fit_string = fitters[field_num][i].fit_string
                    #print ('Showing text: ' + fit_string)
                    ax.text(0.0,-0.5 + 0.1 * (fit_position_index) ,fit_string, color = color)
                    fit_position_index = fit_position_index + 1



        fit_funct = lambda zs, A, mean, sig, shift: A * np.exp(-(mean - zs) ** 2.0 / (2.0 * sig ** 2.0)) + shift
        extra_z_of_zMeas_funct =  lambda single_z, A, mean, sig: A * np.exp(-(single_z - mean) ** 2.0 / (2.0 * sig ** 2.0))

        #For a redshift due to a gravitational well
        well_depth = lambda delta, R: -delta * OmM * H0 ** 2.0 * R ** 2.0 / 2.0
        phi_of_delta = lambda single_z, central_r, delta, R: well_depth(delta, R) * (2 - (abs(r_of_z_interp(single_z) - central_r) / R) ** 2.0 if abs(r_of_z_interp(single_z) - central_r) < R else R / abs(r_of_z_interp(single_z) - central_r))
        extra_z_of_zMeas_funct_G = lambda single_z, central_z, delta, R: (phi_of_delta(single_z, float(r_of_z_interp(central_z)), delta, R) - phi_of_delta(0.0, float(r_of_z_interp(central_z)), delta, R)) / (speedol) ** 2.0
        fit_funct_G = lambda zs, central_z, delta, R: np.array(getErroneousRedshiftPlot(lambda single_z: extra_z_of_zMeas_funct_G(single_z, central_z, delta, R), zs, zHD = zHD, astro_arch = astro_arch, r_of_z_interp = r_of_z_interp ))

        n_G_fit_params = 3
        G_bounds = [(0.0, 4.0), (-1.0, 10.0), (50.0, 5000.0)]
        fit_central_zs = np.linspace(min_z, 1.0, 15)
        fit_deltas = np.linspace(-1.0, 1.0, 6)
        fit_Rs = np.linspace(100.0, 1000.0, 6) #in Mpc
        #G_method = 'bounded'

        #for a redshift due to a velocity profile (v given in km/s)
        redshift_of_v_funct = lambda beta: np.sign(beta) * (np.sqrt((1+abs(beta))/(1-abs(beta))) - 1)
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
        #delta_portion_funct = lambda single_z, central_z, delta_s, rs, c, b_incidence: rs ** 1.0 * (rs / r_well_from_geometry(r_of_z_interp(single_z), r_of_z_interp(central_z), b_incidence * rs)) ** 2.0 * delta_s * (int_of_rsqr_delta_nfw( min(c, r_well_from_geometry(r_of_z_interp(single_z), r_of_z_interp(central_z), b_incidence * rs) / rs), c) ) if (single_z != central_z or abs(b_incidence) > 0.0) else 0.0
        delta_portion_funct = lambda single_z, central_z, delta_s, rs, c, b_incidence: rs ** 1.0 * (rs / r_well_from_geometry(r_of_z_interp(single_z), r_of_z_interp(central_z), b_incidence * rs)) ** 2.0 * delta_s * (int_of_rsqr_delta_nfw( min(c, r_well_from_geometry(r_of_z_interp(single_z), r_of_z_interp(central_z), b_incidence * rs) / rs), c) ) if (single_z != central_z or abs(b_incidence) > 0.0) else 0.0

        #for Burkert; r_s should be given in Mpc
        #burkert_unitless_int = lambda s: 0.25 * (np.log(s ** 2.0 + 1) + 2.0 * np.log(s + 1.0) - 2.0 * np.arctan(s))
        #int_of_rsqr_delta = lambda single_r, delta_s, rs, c: delta_s * rs ** 3.0 * (burkert_unitless_int(single_r / rs) - 1.0 / 3.0 * (single_r / rs)** 3.0 / ((1.0 + c) * (1.0 + c ** 2.0)) )
        #delta_portion_funct = lambda single_z, central_z, delta_s, rs, c, b_incidence: rs ** 1.0 * (rs / r_well_from_geometry(r_of_z_interp(single_z), r_of_z_interp(central_z), b_incidence * rs)) ** 2.0 * delta_s * (burkert_unitless_int(min(c, r_well_from_geometry(r_of_z_interp(single_z), r_of_z_interp(central_z), b_incidence * rs) / rs)) - min(r_well_from_geometry(r_of_z_interp(single_z), r_of_z_interp(central_z), b_incidence * rs) / rs, c) ** 3.0 / (3.0 * (1.0 + c) * (1.0 + c ** 2.0))) if (single_z != central_z or abs(b_incidence) > 0.0) else 0.0

        #the velocity field should be in km/s => so I need to add a conversion from Mpc to km
        velocity_field_funct = lambda single_z, central_z, delta_s, rs, c, b_incidence: - H0 / speedol * H_of_z(single_z) * d_z_D_growth_over_D(single_z) * delta_portion_funct(single_z, central_z, delta_s, rs, c, b_incidence) * sin_r_angle_from_geometry(r_of_z_interp(single_z), r_of_z_interp(central_z), rs * b_incidence)
        #velocity_field_funct = lambda single_z, central_z, beta_peak, beta_peak_ratio, r_width: (  beta_peak * (1 if single_z < central_z else beta_peak_ratio) * (- r_of_z_interp(single_z) + r_of_z_interp(central_z)) / (r_width * np.exp(-0.5)) * np.exp(-abs((r_of_z_interp(single_z) - r_of_z_interp(central_z)) / (np.sqrt(2.0) * r_width)) ** 2.0 ))
        #extra_z_of_zMeas_funct = lambda single_z, central_z, v, r_width, power: redshift_of_v_funct(v * np.exp(-abs((r_of_z_interp(single_z) - r_of_z_interp(central_z)) / (np.sqrt(2.0) * r_width)) ** power ))
        #extra_z_of_zMeas_funct = lambda single_z, central_z, v_peak_left, v_peak_ratio, r_width: redshift_of_v_funct((v_peak_left * 1.0 / (r_width * np.exp(-2.0)) * (1 if single_z < central_z else v_peak_ratio) * (r_of_z_interp(central_z) - r_of_z_interp(single_z)) * np.exp(-abs((r_of_z_interp(single_z) - r_of_z_interp(central_z)) / (np.sqrt(2.0) * r_width)) ** 2.0 )))

        #extra_z_of_zMeas_funct_vel = lambda single_z, central_z, beta_peak_left, beta_peak_ratio, r_width: redshift_of_v_funct(velocity_field_funct(single_z, central_z, beta_peak_left, beta_peak_ratio, r_width))
        extra_z_of_zMeas_funct_vel = lambda single_z, central_z, delta_s, rs, c, b_incidence: redshift_of_v_funct(velocity_field_funct(single_z, central_z, delta_s, rs, c, b_incidence))
        #fit_funct = lambda zs, central_z, betaTen, r_width, power: np.array(getErroneousRedshiftPlot(lambda single_z: extra_z_of_zMeas_funct(single_z, central_z, betaTen * speedol, r_width, power), zs, zHD = zHD, astro_arch = astro_arch, r_of_z_interp = r_of_z_interp ))
        #fit_funct = lambda zs, central_z, betaTen, r_width: np.array(getErroneousRedshiftPlot(lambda single_z: extra_z_of_zMeas_funct(single_z, central_z, betaTen * speedol, r_width, 2.0), zs, zHD = zHD, astro_arch = astro_arch, r_of_z_interp = r_of_z_interp ))
        fit_funct_vel = lambda zs, central_z, delta_s, rs, c, b_incidence: np.array(getErroneousRedshiftPlot(lambda single_z: extra_z_of_zMeas_funct_vel(single_z, central_z, delta_s, rs, c, b_incidence), zs, zHD = zHD, astro_arch = astro_arch, r_of_z_interp = r_of_z_interp ))
        #f2, axarr = plt.subplots(9,1)

        #axarr[0].plot(np.linspace(0.0, 1.0, 101), [H_of_z(z) for z in np.linspace(0.0, 1.0, 101)])
        #axarr[1].plot(np.linspace(0.0, 1.0, 101), [d_z_D_growth_over_D(z) for z in np.linspace(0.0, 1.0, 101)])
        #axarr[2].plot(np.linspace(0.0, 1.0, 101), [delta_portion_funct(z, 0.5, 300.0, 10.0, 10.0, 0.5) for z in np.linspace(0.0, 1.0, 101)])
        #axarr[3].plot(np.linspace(0.0, 1.0, 101), [r_of_z_interp(z) for z in np.linspace(0.0, 1.0, 101)])
        #axarr[4].plot(np.linspace(0.0, 1.0, 101), [sin_r_angle_from_geometry(r_of_z_interp(z), r_of_z_interp(0.5), 300.0 * 0.5) for z in np.linspace(0.0, 1.0, 101)])
        #axarr[5].plot(np.linspace(0.0, 1.0, 101), [velocity_field_funct(z, 0.5, 1.0, 300.0, 10.0, 0.5) for z in np.linspace(0.0, 1.0, 101)] )
        #axarr[6].plot(np.linspace(0.0, 1.0, 101), [redshift_of_v_funct(velocity_field_funct(z, 0.5, 1.0, 300.0, 10.0, 0.5)) for z in np.linspace(0.0, 1.0, 101)] )
        #axarr[7].plot(np.linspace(0.0, 1.0, 101), [extra_z_of_zMeas_funct_vel(z, 0.5, 1.0, 300.0, 10.0, 0.5) for z in np.linspace(0.0, 1.0, 101)] )
        #axarr[8].plot(np.linspace(0.0, 1.0, 101), fit_funct_vel(np.linspace(0.0, 1.0, 101), 0.5, 1.0, 300.0, 10.0, 0.5)  )
        #plt.show()
        max_z = max(all_zs)
        min_z = min(all_zs)
        vel_bounds = ([-0.0, 0.0, 0.4, 50.0], [1.0, 0.1, 1.0, 500.0])
        vel_bounds = [(min_z, 1.5), (-1.0, 50.0), (10.0, 500.0), (1.0, 30.0), (0.01, 3.0)]
        n_vel_fit_params = 5
        fit_central_zs = np.linspace(min_z, 0.5, 11)
        fit_delta_ss = np.linspace(0.0, 20.0, 5)
        fit_rss = np.linspace(10.0, 500.0, 5)
        fit_cs = np.linspace(5.0, 10.0, 2)
        fit_b_incidences = np.linspace(0.5, 3.0, 5)
        #vel_method =


        bounds = vel_bounds # G_bounds; vel_bounds
        n_fit_params = n_vel_fit_params# n_G_fit_params; n_vel_fit_params

        fit_A = fit_central_zs #fit_central_zs
        fit_B = fit_delta_ss #fit_deltas; fit_beta_peaks
        fit_C = fit_rss# fit_Rs; fit_beta_peak_ratios
        fit_D = fit_cs #None; fit_rWidths
        fit_E = fit_b_incidences
        fit_funct = fit_funct_vel # fit_funct_G, fit_funct_vel



        #f2, axarr10 = plt.subplots(3)
        #axarr10[0].plot(zs_to_plot, [velocity_field_funct (z_to_plot, 1.0, 0.04, 1.0, 500.0) for z_to_plot in zs_to_plot])
        #axarr10[0].plot(zs_to_plot, [velocity_field_funct (z_to_plot, 1.0, 0.06, 1.0, 500.0) for z_to_plot in zs_to_plot])
        #axarr10[0].plot(zs_to_plot, [velocity_field_funct (z_to_plot, 1.0, 0.09, 1.0, 500.0) for z_to_plot in zs_to_plot])
        #axarr10[1].plot(zs_to_plot, [extra_z_of_zMeas_funct (z_to_plot, 1.0, 0.04, 1.0, 500.0) for z_to_plot in zs_to_plot])
        #axarr10[1].plot(zs_to_plot, [extra_z_of_zMeas_funct (z_to_plot, 1.0, 0.06, 1.0, 500.0) for z_to_plot in zs_to_plot])
        #axarr10[1].plot(zs_to_plot, [extra_z_of_zMeas_funct (z_to_plot, 1.0, 0.09, 1.0, 500.0) for z_to_plot in zs_to_plot])
        #axarr10[2].plot(zs_to_plot, fit_funct (zs_to_plot, 1.0, 0.04, 1.0, 500.0) )
        #axarr10[2].plot(zs_to_plot, fit_funct (zs_to_plot, 1.0, 0.06, 1.0, 500.0) )
        #axarr10[2].plot(zs_to_plot, fit_funct (zs_to_plot, 1.0, 0.09, 1.0, 500.0) )
        #plt.show()
        #fit_funct = lambda zs, A, mean, sig: A * np.exp(-(mean - zs) ** 2.0 / (2.0 * sig ** 2.0))

        maxfev = 10000
        n_zs_to_display_fit = 1001

        if len(all_zs) > n_fit_params + 2:
            #Try to force the fitting function to sample the entire z-range
            z_bin_centers = np.linspace(min_z, max_z, n_z_bins) + (max_z - min_z) / n_z_bins
            null_dof = len(all_zs) - 1
            dof = len(all_zs) - (n_fit_params + 1)
            weighted_mean = np.sum(np.array(all_resids) * np.array(all_errs) ** -2.0) / np.sum( np.array(all_errs) ** -2.0)
            null_rchisqr = np.sum((np.array( all_resids) - weighted_mean) / np.array(all_errs)) ** 2.0 / null_dof

            #brute_minimization = []
            #brute_chisqr = np.inf
            #start = time.time()
            #for fit_z in fit_central_zs:
            #    print ('working on fit_z = ' + str(fit_z))
            #    for fit_beta in fit_betas:
            #        for fit_rwidth in fit_rWidths:
            #        #    for fit_power in fit_powers:
            #                #this_fit = [fit_z, fit_beta, fit_rwidth, fit_power]
            #                this_fit = [fit_z, fit_beta, fit_rwidth]
            #                this_fit_resids =  fit_funct(all_zs, *this_fit)
            #                this_fit_mean = c.weighted_mean(this_fit_resids - all_resids, all_errs)
            #                #print ('this_fit_mean = ' + str(this_fit_mean ))
            #                this_fit_chisqr = printingFitFunct(fit_funct, all_zs, all_resids, all_errs, np.array(this_fit))
            #                #print ('this_fit + [this_fit_chisqr] = ' + str(this_fit + [this_fit_chisqr]))
            #                #this_fit_rchisqr_argmin = np.argmin(this_fit_rchisqrs)
            #                if this_fit_chisqr < brute_chisqr:
            #                   brute_chisqr = this_fit_chisqr
            #                   brute_minimization = this_fit[:]
            #end = time.time()
            #print ('Brute minimization took ' + str(end - start) + 's')
            #print ('[brute_minimization, brute_chisqr] = ' + str([brute_minimization, brute_chisqr]))
            #print ('Now trying to refine using curve_fit...')
            #true_fit_rchisqr = brute_chisqr
            #funct_fit = brute_minimization
            #try:
            #    print ('brute_minimization = ' + str(brute_minimization))
            #    print ('bounds = ' + str(bounds))
            #    #final_minimization = scipy.optimize.curve_fit(fit_funct, all_zs, all_resids, p0 = brute_minimization, sigma=all_errs, maxfev = maxfev, bounds = bounds)
            #    final_minimization = scipy.optimize.minimize(lambda params: printingFitFunct(fit_funct, all_zs, all_resids, all_errs, params), brute_minimization, bounds = bounds, method='Nelder-Mead')
            #    #fit_res = optimize.curve_fit(lambda zs, A, mu: fit_funct(zs, A, mu, sig, 0.2, 0.0), all_zs, all_resids, p0 = [0.0, max(all_zs) / 2.0], sigma=all_errs, maxfev = maxfev)
            #    this_funct_fit = final_minimization['x']
            #except RuntimeError:
            #    print("Curve_fit failed!.  Plotting initial guess. ")
            #    final_minimization = np.array(brute_minimization)
            this_funct_fit, true_chi_sqr = doMinimization(fit_funct, all_zs, all_resids, all_errs, fit_A, fit_B, fit_C, dof, fit_params_D = fit_D, fit_params_E = fit_E, bounds = bounds)
            #plt.plot(zs_to_plot, fit_funct(zs_to_plot, *this_funct_fit) - printingFitFunct(fit_funct, all_zs, all_resids, all_errs, this_funct_fit)[2] * 0.0, c ='yellow')
            plt.plot(zs_to_plot, fit_funct(zs_to_plot, *this_funct_fit) - printingFitFunct(fit_funct, all_zs, all_resids, all_errs, this_funct_fit)[2], c ='k')
            plt.xlim(0.0, 1.0)
            rand_chi_sqrs = [np.nan for i in range(n_randomizations)]
            rand_fits = [[] for i in range(n_randomizations)]
            n_rand_done = 0
            while n_rand_done <  n_randomizations:
                start = time.time()
                if (n_rand_done % 1 == 0 and n_rand_done > 0): print ('Working on randomization ' + str(n_rand_done))
                rand_resid_indeces = list(range(len(all_resids)))
                np.random.shuffle(rand_resid_indeces)
                rand_resids = [all_resids[index] for index in rand_resid_indeces]
                rand_errs = [all_errs[index] for index in rand_resid_indeces]
                rand_fits[n_rand_done], rand_chi_sqrs[n_rand_done] = doMinimization(fit_funct, all_zs, rand_resids, rand_errs, fit_central_zs, fit_beta_peaks, fit_beta_peak_ratios, dof, fit_params_D = fit_rWidths, bounds= bounds)
                #brute_rchisqr = np.inf
                #brute_minimization = []
                #for fit_A in fit_As:
                #    for fit_mu in fit_mus:
                #        for fit_sig in fit_sigs:
                #            this_fit = [fit_A, fit_mu, fit_sig, 0.0]
                #            this_fit_resids =  fit_funct(all_zs, *this_fit)
                #            this_fit_mean = c.weighted_mean(this_fit_resids - rand_resids, rand_errs)
                #            #print ('this_fit_mean = ' + str(this_fit_mean ))
                #            this_fit_resids = this_fit_resids - this_fit_mean
                #            this_fit[-1] = -this_fit_mean
                #            this_fit_rchisqr = np.sum(np.array(((this_fit_resids - rand_resids) / np.array(rand_errs)) ** 2.0)) / dof
                #            #this_fit_rchisqr_argmin = np.argmin(this_fit_rchisqrs)
                #            if this_fit_rchisqr < brute_rchisqr:
                #               #print ('this_fit: this_fit_rchisqr = ' + str(this_fit) + ':' + str(this_fit_rchisqr))
                #               brute_rchisqr = this_fit_rchisqr
                #               brute_minimization = this_fit
                print('[true_chi_sqr, rand_chi_sqrs[n_rand_done]] = ' + str([true_chi_sqr, rand_chi_sqrs[n_rand_done]]))
                n_rand_done = n_rand_done + 1
                end = time.time()
                if (n_rand_done % 1 == 0 and n_rand_done > 0):  print (str(n_rand_done) + 'th randomization took ' + str(end - start) + 's')
                #rand_resids = np.random.normal([0.0 for resid in all_resids], all_errs)
                #print ('rand_resids = ' +str(rand_resids))

                #try:
                #    rand_fit_res = optimize.curve_fit(fit_funct, all_zs, rand_resids, p0 = initial_guess, sigma=rand_errs, maxfev = maxfev)
                #    rand_funct_fit = rand_fit_res[0]
                #    rand_fit_resids = fit_funct(all_zs, *rand_funct_fit)
                #    #null_rand_rchisqr = np.sum(np.array(((rand_resids) / np.array(rand_errs)) ** 2.0)) / null_dof
                #    #print ('null_rand_rchisqr = ' + str(null_rand_rchisqr))
                #    rand_fit_rchisqr = np.sum(np.array(((rand_fit_resids - rand_resids) / np.array(rand_errs)) ** 2.0)) / dof
                #    rand_rchiSqrs[n_rand_done] = rand_fit_rchisqr
                #    n_rand_done = n_rand_done + 1
                #except RuntimeError:
                #    n_rand_done = n_rand_done
                #    print("Curve_fit to randomized data failed!. ")
        if len(all_zs) > 0:
            #if n_randomizations < 1 or len(all_zs) <= len(initial_guess) + 2:
            label_str = field_str + ': best fit is ' + str( [cant.round_to_n(fitted_val,3) for fitted_val in this_funct_fit]) + ' with chisqr ' + str(cant.round_to_n(true_chi_sqr, 5))
            if n_randomizations > 0:
                label_str = label_str + '\n' + str(cant.round_to_n(len([rand_chi_sqr for rand_chi_sqr in rand_chi_sqrs if rand_chi_sqr< true_chi_sqr]) / n_randomizations * 100, 3)) + '% chance randomized data have better fit.' r'$\chi^2$=' + str(cant.round_to_n(true_chi_sqr, 4))
            print ('label_str = ' + str(label_str))
            ax.text(0.5, 0.8, label_str,  horizontalalignment='center', transform=ax.transAxes, fontsize = 7)

            #print ('[null_rchisqr, true_fit_rchisqr, np.mean(rand_rchiSqrs), np.std(rand_rchiSqrs)] = ' + str([null_rchisqr, true_fit_rchisqr, np.mean(rand_rchiSqrs), np.std(rand_rchiSqrs)]))
            print (str(len([rand_chi_sqr for rand_chi_sqr in rand_chi_sqrs if rand_chi_sqr < true_chi_sqr])) + ' of ' + str(n_randomizations) + ' randomizations can be better fit than the the true data.')


        #ax.set_xlim(z_range[0], z_range[1])
        #axarr[1].set_xlim(z_range[0], z_range[1])

        ax.set_ylim(res_range[0], res_range[1])
        #axarr[1].set_ylim(binned_res_range[0], binned_res_range[1])

        if separate_fields_by_plot or plot_index[0] == fig_side - 1:
            ax.set_xlabel('z')
        #axarr[1].set_ylabel('Binned mu residual')
        if separate_fields_by_plot or plot_index[1] == 0:
            ax.set_ylabel(r'$\Delta \mu$')

        #ax.set_xticks([-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1])

        plt.suptitle(r'$\Delta \mu$ of Pan Stars 1 medium-deep fields' )

        #print ('len(survey_plots) = ' + str(len(survey_plots)))
        #plt.legend(survey_plots, surveys_to_display, loc = 'upper left', prop = {'size':8})

        master_end = time.time()
        print ('Whole thing took ' + str(master_end - master_start) + 's')
        if separate_fields_by_plot:
            if save: plt.savefig(plot_dir + 'SN_residuals_v_z_PS1MD_field' + str(field) + '_bin_' + bin_scheme + str(n_bins) + '_fit_' + fit_information['funct'] + '.pdf')
            if show: plt.show()

    if not(separate_fields_by_plot):
        if save: plt.savefig(plot_dir + 'SN_residuals_v_z_PS1MD_fields_bin_' + bin_scheme + str(n_bins) + '_fit_' + fit_information['funct'] + '.pdf')
        if show: plt.show()
    plt.close('all')


    def __init__(self, data_set,n_z_bins = 10, zHD = 0, OmM = 0.3, OmL = 0.7, Om0 = 1.0, OmR = 0.0, H0 = 70.0, interp_z_params = [0.0, 100.0, 1001]):

        self.data_set = data_set
        self.n_z_bins = n_z_bins
        self.zHD = zHD
        self.OmM = OmM
        self.OmL = OmL
        self.Om0 = Om0
        self.OmR = OmR
        self.H0 = H0
        self.interp_z_params = interp_z_params 

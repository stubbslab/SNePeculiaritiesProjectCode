import cantrips as c
import numpy as np
import loadSN as lsn
import scipy.integrate as integrate
import CosmologicalParameterArchive as cpa
import matplotlib.pyplot as plt
import time
import matplotlib.gridspec as gridspec
import binData as bd
from matplotlib.colors import LogNorm
import randomSimulationFunctions as rSimF
from scipy.ndimage import gaussian_filter1d

#Hi.
#This function should take in a canonical cosmology with tunable cosmological parameters OmM, OmL, H0, Om0, and OmR
# and provide the set of the parameters that you need to compare measured and predicted autocorrelations.
# Specifically, it returns z_i, z_j, r_i, r_j, theta_i, theta_j, r_sep, r_perp, r_par, r_cent, ang_seps, xi, xi_errs.
# The returned parameters are chosen to match the predictions of Biern and Yoo 2017, https://arxiv.org/pdf/1704.07380.pdf and Gordon, Land, and Slosar 2007, https://arxiv.org/pdf/0705.1718.pdf
# The parameter descriptions are as follows, for SNe pair i_j.  There are (N_sn X (N_sn - 1))/2  = 1048 X 1047 / 2 individual pairs.  That's 548628  for the full data set
#  z_i = redshift of SNe i
#  z_i = redshift of SNe j
#  r_i = comoving distance to SNe i
#  r_i = comoving distance to SNe j
#  theta_i = the angle between the comoving vector to SNe i and the comoving vector connecting SNe j to SNe i (see Gordon 2007)
#  theta_i = the angle between the comoving vector to SNe j and the comoving vector connecting SNe j to SNe i (see Gordon 2007)
#  r_sep = the comoving separation between SNe i and SNe j
#  r_perp = the comoving separation between SNe i and SNe j perpendicular to line of sight (see Biern 2017)
#  r_par = the comoving separation between SNe i and SNe j perpendicular to line of sight (see Biern 2017)
#  r_center = the vector that the connects the observer to the center of the comoving separation between SNe i and SNe j (see Biern 2017, where it's called r_c)
#  ang_seps = the angular separation between SNe i and SNe j on the sky (in radians)
#  xi = the luminosity distance residual correlation (\delta d_{L,i} \times \delta d_{L,j}).  This is the signal.
#  xi_errs = the luminosity distance residual correlation uncertainty.  This is the uncertainty in the signal.
#Move all of Sasha's scripts into a directory (or make sure they are in the Python path).
#Then
# $ python
# Python 3.7.6 (default, Jan 26 2020, 20:40:08)
# [Clang 11.0.0 (clang-1100.0.33.12)] on darwin
# Type "help", "copyright", "credits" or "license" for more information.
# >>> import MeasureSNAutocorrelation as msna
# >>> OmM, OmL, H0, OmR, Om0 = [0.3, 0.7, 100, 0.0, 1.0]
# >>> z_i,  z_j,  r_i, r_j, theta_i,  theta_j,  r_sep, r_perp, r_par, r_cent, ang_seps, sky_coords_i, sky_coords_j, xi, xi_errs =  msna.measureAutocorrelation(OmM = OmM, OmL = OmL, H0 =  H0, OmR = OmR, Om0 = Om0)
# Note: The given value of H0 DOES NOT MATTER.  A single overall distance modulus residual is subtracted from the data so that the weighted mean of \delta \mu values is 0
#       These SNe (high and low redshift) are not useful for constraining H0 because the uncertainty in the true value of their normalized luminosity is comparatively large.
# You can enter in your preferred values of OmL, OmM, H0, and OmR

def measureAutocorrelation(normal_randomize = 0, rearrange_randomize = 0, OmM = 'default', OmL = 'default', H0 = 'default', OmR = 'default', Om0 = 'default', zHD = 0, surveys_of_interest = ['all'], surveys_to_excise = [], pull_extinctions = 0, z_include_range = [-np.inf, np.inf], remove_residuals = 1, param_to_correlate_str = 'dl'):
    cosmo_arch = cpa.CosmologicalParameterArchive(params_source = 'pantheon')
    if OmM in ['default']:
        OmM = cosmo_arch.getOmegaM()[0]
    if OmL in ['default']:
        OmL = cosmo_arch.getOmegaLambda()[0]
    if OmR in ['default']:
        OmR = cosmo_arch.getOmegaR()[0]
    if Om0 in ['default']:
        Om0 = cosmo_arch.getOmega0()[0]
    if H0 in ['default']: #Hubble constant in km/s/Mpc
        H0 = cosmo_arch.getH0()[0]
    speedOfLight = cosmo_arch.getc() # speed of light in km/s

    dl_scaling = speedOfLight / H0
    all_sn = lsn.loadSN(1, zHD = zHD, surveys_of_interest = surveys_of_interest, pull_extinctions = pull_extinctions, OmM = OmM, OmL = OmL, OmR = OmR, H0 = H0, Om0 = Om0 )
    all_sn = [sn for sn in all_sn if (sn['RA'] > 0.0 or sn['Dec'] > 0.0)]
    all_sn = [sn for sn in all_sn if not(sn['survey'] in surveys_to_excise)]
    all_sn = [sn for sn in all_sn if (sn['z'] > z_include_range[0] and sn['z'] < z_include_range[1])]
    n_sn = len(all_sn)
    print ('Number of included SN is: ' + str(n_sn))
    all_zs = [sn['z'] for sn in all_sn]
    all_ras = [sn['RA'] for sn in all_sn]
    all_decs = [sn['Dec'] for sn in all_sn]
    all_mus = [sn['mu'] for sn in all_sn]
    all_muErrs = [sn['muErr'] for sn in all_sn]
    all_muDiffs = [sn['muDiff'] for sn in all_sn]
    if normal_randomize:
        all_muDiffs = [np.random.normal(0.0, mu_err) for mu_err in all_muErrs]
    elif rearrange_randomize:
        all_muDiffs = [sn['muDiff'] for sn in all_sn]
        all_nSigs = np.array(all_muDiffs) / np.array(all_muErrs)
        rand_zs, rand_muDiffs, rand_muErrs  = rSimF.randomSortNSigma(all_zs, all_muDiffs, all_muErrs, [0.0 for mu in all_mus])
        all_zs = rand_zs
        all_muDiffs = rand_muDiffs
        all_muErrs = rand_muErrs

    all_muPreds = [sn['muPred'] for sn in all_sn]
    #If you want to measure residuals relativ
    #all_muDiffs = [all_mus[i] - all_muPreds[i] for i in range(len(all_sn))]
    all_dls = [sn['dl'] for sn in all_sn]
    all_dlErrs = [sn['dlErr'] for sn in all_sn]

    if remove_residuals:
        weights = 1.0 / (np.array(all_muErrs) ** 2.0)
        mu_offset = (np.sum(np.array(all_muDiffs) * np.array(weights)) / np.sum(weights) )
        print ('mu_offset = ' + str(mu_offset))
        all_muDiffs = [diff - mu_offset for diff in all_muDiffs]

    all_dLDiffs = [10.0 ** (muDiff / 5) - 1 for muDiff in all_muDiffs]
    all_dLDiffErrs = [all_muErrs[i] * np.log(10) / 5 * abs((all_dLDiffs[i] + 1)) for i in range(len(all_muDiffs))]

    Hubble_funct = lambda zp: np.sqrt(OmM * (1.0 + zp) ** 3.0 + OmL + OmR * (1.0 + zp) ** 4.0 + (1.0 - OmM - OmR - OmL) * (1.0 + zp) ** 2.0 )
    dl_funct = lambda z: dl_scaling * (1.0 + z) * integrate.quad(lambda zp: 1.0 / Hubble_funct(zp) , 0.0, z)[0]
    #exp_dls = [dl_funct(sn_z) for sn_z in all_zs]
    #exp_mus = [25.0 + 5.0 * np.log10(exp_dl / 1.0) for exp_dl in exp_dls]
    #delta_dLs = [(all_dls[i] - exp_dls[i]) / exp_dls[i] for i in range(n_sn)]
    #delta_dLErrs = [all_dlErrs[i] / exp_dls[i] for i in range(n_sn)]
    #f, axarr = plt.subplots(1,1, squeeze = False)
    #axarr[0, 0].scatter(all_zs, 10.0 ** (np.array(np.array(all_mus) - np.array(exp_mus)) / 5.0) - 1, marker = '.', c = 'b')
    #axarr[0, 0].errorbar(all_zs, 10.0 ** (np.array(np.array(all_mus) - np.array(exp_mus)) / 5.0) - 1, y_errs = delta_dlErrs, fmt = 'none', c = 'b')
    #param_to_correlate = [(10.0 ** (all_muDiffs[i] / 5.0) - 1.0) * (1.0 - (1.0 + all_zs[i]) ** 2.0 / (Hubble_funct(all_zs[i]) * exp_dls[i] / dl_scaling) ) ** -1.0 for i in range(len(all_muDiffs)) ]

    #These did not seem to work very well
    #param_to_correlate = [10.0 ** (all_muDiff / 5.0) - 1 for all_muDiff in all_muDiffs]
    #param_errs = [np.abs(np.log(10.0) / 5.0 * (10.0 ** (all_muDiffs[i] / 5.0) - 1) * all_muErrs[i]) for i in range(len(all_muErrs))]

    if param_to_correlate_str.lower() == 'dmu':
        param_to_correlate = all_muDiffs
        param_errs = all_muErrs
    else:
        param_to_correlate = all_dLDiffs
        param_errs = all_dLDiffErrs

    #param_to_correlate = np.array(all_dLDiffs) / np.array(all_dLDiffErrs)
    #para_errs = [1.0 for elem in param_to_correlate]


    print ('computing comoving separations....')
    print ('[speedOfLight, H0] = ' + str([speedOfLight, H0] ))
    comoving_dists = [speedOfLight / H0 * integrate.quad(lambda z_int: 1.0 / Hubble_funct(z_int) , 0 , z)[0] for z in all_zs ]

    print ('computing autocorrelations....')
    start = time.time()
    autocorrelations = c.flattenListOfLists([[param_to_correlate[i] * param_to_correlate[j] for j in range(i+1, n_sn) ] for i in range(n_sn)] )
    autocorrelation_errs = c.flattenListOfLists([[ np.sqrt((param_to_correlate[i] * param_errs[j]) ** 2.0 + (param_to_correlate[j] * param_errs[i]) ** 2.0) for j in range(i+1, n_sn) ] for i in range(n_sn)] )
    end = time.time()
    print ('Took ' + str(end -start) + 's')
    print ('computing angular separations....')
    start = time.time()
    z_i = c.flattenListOfLists([[all_zs[i]             for j in range(i+1, n_sn)] for i in range(n_sn)])
    z_j = c.flattenListOfLists([[all_zs[j]             for j in range(i+1, n_sn)] for i in range(n_sn)])
    r_i = np.array(c.flattenListOfLists([[comoving_dists[i]     for j in range(i+1, n_sn)] for i in range(n_sn)]))
    r_j = np.array(c.flattenListOfLists([[comoving_dists[j]     for j in range(i+1, n_sn)] for i in range(n_sn)]))
    ra_i = np.array(c.flattenListOfLists([[all_ras[i]   for j in range(i+1, n_sn)] for i in range(n_sn)]))
    dec_i = np.array(c.flattenListOfLists([[all_decs[i] for j in range(i+1, n_sn)] for i in range(n_sn)]))
    ra_j = np.array(c.flattenListOfLists([[all_ras[j]   for j in range(i+1, n_sn)] for i in range(n_sn)]))
    dec_j = np.array(c.flattenListOfLists([[all_decs[j] for j in range(i+1, n_sn)] for i in range(n_sn)]))
    ang_seps = np.array( [c.measureAngularSeparationOnSky([ra_i[k], dec_i[k]], [ra_j[k], dec_j[k]], return_radian = 1) for k in range(len(ra_i))] )
    end = time.time()
    print ('Took ' + str(end - start) + 's')
    #comoving_products = [[comoving_dists[i] * comoving_dists[j] for j in range(i+1, len(all_ras))] for i in range(len(all_ras))]
    comoving_seps = np.sqrt(r_i ** 2.0 + r_j ** 2.0 - 2.0 * np.cos(ang_seps) * r_i * r_j)
    comoving_centers = 0.5 * np.sqrt(r_i ** 2.0 + 2.0 * r_i * r_j * np.cos(ang_seps) + r_j ** 2.0 )
    comoving_parallels = 0.5 * abs(r_i ** 2.0 - r_j ** 2.0) / comoving_centers
    comoving_perps = np.sqrt(2.0 * (r_i **2.0 + r_j ** 2.0) - 4 * comoving_centers ** 2.0 - comoving_parallels ** 2.0)

    print ('computing comoving_seps....')
    cos_theta_i = (r_i - np.cos(ang_seps) * r_j) / (comoving_seps)
    theta_i = np.arccos(cos_theta_i)
    cos_theta_j = (-r_j + np.cos(ang_seps) * r_i) / (comoving_seps)
    theta_j = np.arccos(cos_theta_j)
    #print ('(thetajs - thetais) - np.array(ang_seps) = ' + str((thetajs - thetais) - np.array(ang_seps)))
    #comoving_seps = c.flattenListOfLists([[ comoving_dists[i] ** 2.0 + comoving_dists[j] ** 2.0 - 2.0 * np.cos(ang_seps[j - (i+1) + sum([len(comoving_dists) - (iprime + 1) for iprime in  range(i) ]) ]) * comoving_dists[i] * comoving_dists[j]
    #                                            for j in range(i+1, len(comoving_dists)) ] for i in range(len(param_to_correlate)) ] )
    #end = time.time()
    #print ('Took ' + str(end -start) + 's')

    xi = autocorrelations
    xi_errs = autocorrelation_errs
    r_sep = comoving_seps
    r_perp = comoving_perps
    r_par = comoving_parallels
    r_cent = comoving_centers
    sky_coords_i = [[ra_i[k], dec_i[k]] for k in range(len(ra_i))]
    sky_coords_j = [[ra_j[k], dec_j[k]] for k in range(len(ra_j))]
    #return [z_is, z_js, ris, rjs, theta_is, theta_js, r_s,   r_perp, r_par, r_cent, ang_seps, xi_s, xi_errs]
    return [z_i,  z_j,  r_i, r_j, theta_i,  theta_j,  r_sep, r_perp, r_par, r_cent, ang_seps, sky_coords_i, sky_coords_j, xi, xi_errs]


def plotAutocorrelation(my_xis = None, my_xi_errs = None,
                        OmM = 0.305, OmL = 0.695, H0 = 'default', OmR = 'default', zHD = 0, correlation_color_lims = [None, None],
                        surveys_of_interest = ['all'], surveys_to_excise = [], pull_extinctions = 0, n1dBins_set = [100], gauss_smoothing = [0, 0.01],
                        ylims = [-0.0005, 0.0005], figsize = [14, 8],
                        img_n_pixels = [100, 100], ang_sep_cuts = [[0.0, np.pi ]],
                        center_sky_coords = [[90, 0.0]], anisotropy_search_rad = 180.0, #The parameters to vary for looking for anisotropy on the sky
                        make_1d_plots = 1, make_2d_plots = 1, show_fig = 1, save_fig =1, save_fig_name = None, param_to_correlate_str = 'dL',
                        save_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/SNTwoPointCorrelationsProject/plots/', bin_value = 'mean',
                        z_cuts = [[-np.inf, np.inf]], mean_z_bounds = None, r_sep_bounds = None, r_perp_bounds = None, r_par_bounds = None,
                        normal_randomize = 0, rearrange_randomize = 0, ):

    z_is, z_js, r_is, r_js, theta_is, theta_js, r_seps, r_perp, r_par, r_cent, ang_seps, sky_coords_i, sky_coords_j, xis, xi_errs = measureAutocorrelation(normal_randomize = normal_randomize, rearrange_randomize = rearrange_randomize, OmM = OmM, OmL = OmL, H0 = H0, OmR = OmR, zHD = zHD, surveys_of_interest = surveys_of_interest, surveys_to_excise = surveys_to_excise, pull_extinctions = pull_extinctions, z_include_range = [-np.inf, np.inf], param_to_correlate_str = param_to_correlate_str)
    print ('[r_seps[0], r_seps[-1], np.mean(r_seps), xis[0], xis[-1], np.mean(xis)] = ' + str([r_seps[0], r_seps[-1], np.mean(r_seps), xis[0], xis[-1], np.mean(xis)]))
    if not(my_xis is None):
        xis = my_xis
    if not(my_xi_errs is None):
        xi_errs = my_xi_errs
    n_pairs = len(xis)

    #f, axarr = plt.subplots(2,1)
    #print('[np.mean(r_par), np.median(r_par), np.std(r_par), np.mean(r_perp), np.median(r_perp), np.std(r_perp)] = ' + str([np.mean(r_par), np.median(r_par), np.std(r_par), np.mean(r_perp), np.median(r_perp), np.std(r_perp)]  ))
    #axarr[0].hist(r_perp, bins = 101)
    #axarr[1].hist(r_par, bins = 101)
    #plt.show()
    mean_zs = [ (z_is[k] + z_js[k]) / 2.0  for k in range(len(z_is)) ]
    mean_rs = [ (z_is[k] + z_js[k]) / 2.0  for k in range(len(z_is)) ]

    if mean_z_bounds == None: mean_z_bounds = [min(mean_zs), max(mean_zs)]
    if r_sep_bounds == None: r_sep_bounds = [min(r_seps), max(r_seps)]
    if r_perp_bounds == None: r_perp_bounds = [min(r_perp), max(r_perp)]
    if r_par_bounds == None: r_par_bounds = [min(r_par), max(r_par)]

    fig = plt.figure(tight_layout=True, figsize = figsize )
    #gs = gridspec.GridSpec(2, 3)
    gs = gridspec.GridSpec(4, (make_1d_plots + make_2d_plots) * len(ang_sep_cuts))

    for cut_num in range(len(ang_sep_cuts)):
        n1dBins = n1dBins_set[cut_num]
        ang_sep_cut = ang_sep_cuts[cut_num]
        z_cut = z_cuts[cut_num]
        center_sky_coord = center_sky_coords[cut_num]
        offset_from_center_sky_coord_i = np.array( [c.measureAngularSeparationOnSky(sky_coord, center_sky_coord, return_radian = 0) for sky_coord in sky_coords_i] )
        #print ('offset_from_center_sky_coord_i = ' + str(offset_from_center_sky_coord_i))
        #plt.hist(offset_from_center_sky_coord_i )
        #plt.show()
        offset_from_center_sky_coord_j = np.array( [c.measureAngularSeparationOnSky(sky_coord, center_sky_coord, return_radian = 0) for sky_coord in  sky_coords_j] )
        cut_conditions = [( (ang_seps[i] <= ang_sep_cut[1] and ang_seps[i] >= ang_sep_cut[0])
                             and (z_is[i] <= z_cut[1] and z_is[i] >= z_cut[0] and z_js[i] <= z_cut[1] and z_js[i] >= z_cut[0])
                             and ( (offset_from_center_sky_coord_i[i] < anisotropy_search_rad) and (offset_from_center_sky_coord_j[i] > 180.0 - anisotropy_search_rad ) )
                          ) for i in range(n_pairs)]
        n_cut_pairs = len([condition for condition in cut_conditions if condition])
        print ('n_cut_pairs = ' + str(n_cut_pairs))
        cut_z_is = [z_is[i]         for i in range(n_pairs) if cut_conditions[i] ]
        cut_z_js = [z_js[i]         for i in range(n_pairs) if cut_conditions[i] ]
        cut_mean_zs = [mean_zs[i]   for i in range(n_pairs) if cut_conditions[i]]
        cut_r_is = [r_is[i]         for i in range(n_pairs) if cut_conditions[i] ]
        cut_r_js = [r_js[i]         for i in range(n_pairs) if cut_conditions[i] ]
        cut_mean_rs = [mean_rs[i]   for i in range(n_pairs) if cut_conditions[i] ]
        cut_theta_is = [theta_is[i] for i in range(n_pairs) if cut_conditions[i] ]
        cut_theta_js = [theta_js[i] for i in range(n_pairs) if cut_conditions[i] ]
        cut_r_seps = [r_seps[i]     for i in range(n_pairs) if cut_conditions[i] ]
        cut_r_perp = [r_perp[i]     for i in range(n_pairs) if cut_conditions[i] ]
        cut_r_par = [r_par[i]       for i in range(n_pairs) if cut_conditions[i] ]
        cut_r_cent = [r_cent[i]     for i in range(n_pairs) if cut_conditions[i] ]
        cut_xis = [xis[i]           for i in range(n_pairs) if cut_conditions[i] ]
        cut_xi_errs = [xi_errs[i]   for i in range(n_pairs) if cut_conditions[i] ]
        cut_ang_seps = [ang_seps[i] for i in range(n_pairs) if cut_conditions[i] ]

        sorted_mean_zs_on_mean_zs, sorted_r_seps_on_mean_zs, sorted_r_perp_on_mean_zs, sorted_r_par_on_mean_zs, sorted_xis_on_mean_zs, sorted_xierrs_on_mean_zs = c.safeSortOneListByAnother(cut_mean_zs, [cut_mean_zs, cut_r_seps, cut_r_perp, cut_r_par, cut_xis, cut_xi_errs])
        sorted_mean_zs_on_r_seps, sorted_r_seps_on_r_seps, sorted_xis_on_r_seps, sorted_xierrs_on_r_seps = c.safeSortOneListByAnother(cut_r_seps, [cut_mean_zs, cut_r_seps, cut_xis, cut_xi_errs])
        sorted_mean_zs_on_r_perps, sorted_r_perps_on_r_perps, sorted_xis_on_r_perps, sorted_xierrs_on_r_perps = c.safeSortOneListByAnother(cut_r_perp, [cut_mean_zs, cut_r_perp, cut_xis, cut_xi_errs])
        sorted_mean_zs_on_r_pars, sorted_r_pars_on_r_pars, sorted_xis_on_r_pars, sorted_xierrs_on_r_pars = c.safeSortOneListByAnother(cut_r_par, [cut_mean_zs, cut_r_par, cut_xis, cut_xi_errs])

        print('[sorted_r_seps_on_r_seps[0], sorted_r_seps_on_r_seps[-1], np.mean(sorted_r_seps_on_r_seps), sorted_xis_on_r_seps[0], sorted_xis_on_r_seps[-1], np.mean(sorted_xis_on_r_seps), sorted_xierrs_on_r_seps[0], sorted_xierrs_on_r_seps[-1], np.mean(sorted_xis_on_r_seps)] = ' + str([sorted_r_seps_on_r_seps[0], sorted_r_seps_on_r_seps[-1], np.mean(sorted_r_seps_on_r_seps), sorted_xis_on_r_seps[0], sorted_xis_on_r_seps[-1], np.mean(sorted_xis_on_r_seps), sorted_xierrs_on_r_seps[0], sorted_xierrs_on_r_seps[-1], np.mean(sorted_xis_on_r_seps)] ))
        #print ('[mean_z_bounds, r_sep_bounds, r_perp_bounds, r_par_bounds] = ' + str([mean_z_bounds, r_sep_bounds, r_perp_bounds, r_par_bounds]))
        #print ('[sorted_r_perps_on_r_perps[0], sorted_r_perps_on_r_perps[-1],  sorted_r_pars_on_r_pars[0], sorted_r_pars_on_r_pars[-1]] = ' + str([sorted_r_perps_on_r_perps[0], sorted_r_perps_on_r_perps[-1],  sorted_r_pars_on_r_pars[0], sorted_r_pars_on_r_pars[-1]]))
        #print ('[np.mean(sorted_xis_on_mean_zs), np.std(sorted_xis_on_mean_zs), np.mean(sorted_xis_on_r_seps), np.std(sorted_xis_on_r_seps), np.mean(sorted_xis_on_r_perps), np.std(sorted_xis_on_r_perps), np.mean(sorted_xis_on_r_pars), np.std(sorted_xis_on_r_pars)] = ' + str([np.mean(sorted_xis_on_mean_zs), np.std(sorted_xis_on_mean_zs), np.mean(sorted_xis_on_r_seps), np.std(sorted_xis_on_r_seps), np.mean(sorted_xis_on_r_perps), np.std(sorted_xis_on_r_perps), np.mean(sorted_xis_on_r_pars), np.std(sorted_xis_on_r_pars)]))
        if make_1d_plots:
            if gauss_smoothing[0]:
                binned_zs = np.linspace(mean_z_bounds[0], mean_z_bounds[1], n1dBins)
                binned_xis_on_zs = c.convolveGaussian(sorted_mean_zs_on_mean_zs, sorted_xis_on_mean_zs, binned_zs, (mean_z_bounds[1] - mean_z_bounds[0]) * gauss_smoothing[1], y_errs = sorted_xierrs_on_mean_zs)
                binned_r_seps = np.linspace(r_sep_bounds[0], r_sep_bounds[1], n1dBins)
                binned_xis_on_r_seps = c.convolveGaussian(sorted_r_seps_on_r_seps, sorted_xis_on_r_seps, np.linspace(r_sep_bounds[0], r_sep_bounds[1], n1dBins), (r_sep_bounds[1] - r_sep_bounds[0]) * gauss_smoothing[1], y_errs = sorted_xierrs_on_r_seps)
                binned_r_perps = np.linspace(r_perp_bounds[0], r_perp_bounds[1], n1dBins)
                binned_xis_on_r_perps = c.convolveGaussian(sorted_r_perps_on_r_perps, sorted_xis_on_r_perps, np.linspace(r_perp_bounds[0], r_perp_bounds[1], n1dBins), (r_perp_bounds[1] - r_perp_bounds[0]) * gauss_smoothing[1], y_errs = sorted_xierrs_on_r_perps)
                binned_r_pars = np.linspace(r_par_bounds[0], r_par_bounds[1], n1dBins)
                binned_xis_on_r_pars = c.convolveGaussian(sorted_r_pars_on_r_pars, sorted_xis_on_r_pars, np.linspace(r_par_bounds[0], r_par_bounds[1], n1dBins), (r_par_bounds[1] - r_par_bounds[0]) * gauss_smoothing[1], y_errs = sorted_xierrs_on_r_pars)
            else:
                binned_zs, binned_xis_on_zs = bd.binData(sorted_mean_zs_on_mean_zs, sorted_xis_on_mean_zs, n_bins = n1dBins, y_errs = sorted_xierrs_on_mean_zs, bin_boundaries = mean_z_bounds, computed_value = bin_value)
                binned_r_seps, binned_xis_on_r_seps = bd.binData( sorted_r_seps_on_r_seps, sorted_xis_on_r_seps, n_bins = n1dBins, y_errs = sorted_xierrs_on_r_seps, bin_boundaries = r_sep_bounds, computed_value = bin_value)
                binned_r_perps, binned_xis_on_r_perps = bd.binData( sorted_r_perps_on_r_perps, sorted_xis_on_r_perps, n_bins = n1dBins, y_errs = sorted_xierrs_on_r_perps, bin_boundaries = r_perp_bounds, computed_value = bin_value)
                binned_r_pars, binned_xis_on_r_pars = bd.binData( sorted_r_pars_on_r_pars, sorted_xis_on_r_pars, n_bins = n1dBins, y_errs = sorted_xierrs_on_r_pars, bin_boundaries = r_par_bounds, computed_value = bin_value)

                print ('[np.mean(binned_xis_on_zs[0]), np.std(binned_xis_on_zs[0]), np.mean(binned_xis_on_r_seps[0]), np.std(binned_xis_on_r_seps[0]), np.mean(binned_xis_on_r_perps[0]), np.std(binned_xis_on_r_perps[0]), np.mean(binned_xis_on_r_pars[0]), np.std(binned_xis_on_r_pars[0])] = ' + str([np.mean(binned_xis_on_zs[0]), np.std(binned_xis_on_zs[0]), np.mean(binned_xis_on_r_seps[0]), np.std(binned_xis_on_r_seps[0]), np.mean(binned_xis_on_r_perps[0]), np.std(binned_xis_on_r_perps[0]), np.mean(binned_xis_on_r_pars[0]), np.std(binned_xis_on_r_pars[0])]))
            #print ('binned_xis_on_zs = ' + str(binned_xis_on_zs))
        #f, axarr = plt.subplots(2,2)
        #axarr[0,0].hist(binned_xis_on_zs[0], bins = 101)
        #axarr[1,0].hist(binned_xis_on_r_seps[0], bins = 101)
        #axarr[0,1].hist(binned_xis_on_r_perps[0], bins = 101)
        #axarr[1,1].hist(binned_xis_on_r_pars[0], bins = 101)
        #plt.show()

        #comoving_bins, z_bins, binned_corrs, binned_corr_errs = c.binDataOfTwoVariables(sorted_r_seps_on_mean_zs[0:len(sorted_r_seps_on_mean_zs)//2], sorted_mean_zs_on_mean_zs[0:len(sorted_r_seps_on_mean_zs)//2], sorted_xis_on_mean_zs[0:len(sorted_r_seps_on_mean_zs)//2], *img_n_pixels)
        if make_2d_plots:
            comoving_bins, z_bins, binned_corrs, avg_binned_corrs, avg_binned_corr_errs, n_corrs_in_bin = c.binDataOfTwoVariables(sorted_r_seps_on_mean_zs, sorted_mean_zs_on_mean_zs, sorted_xis_on_mean_zs, *img_n_pixels, x_bin_boundaries = r_sep_bounds, y_bin_boundaries = mean_z_bounds )
        else:
            comoving_bins = np.linspace(*r_sep_bounds, img_n_pixels[0])
            z_bins = np.linspace(*mean_z_bounds, img_n_pixels[1])
        #print ('np.shape(binned_corrs) = ' + str(np.shape(binned_corrs)))
        #n_corrs_in_bin = np.sum(np.zeros(np.shape(binned_corrs)) + 1.0, axis = 0)
        #print ('np.shape(n_corrs_in_bin) = ' + str(np.shape(n_corrs_in_bin)))

        if param_to_correlate_str.lower()  == 'dmu':
            correlation_str = r'$\delta \mu_i \times \delta \mu_j$'
        else:
            correlation_str = r'$\delta d_{L,i} \times \delta d_{L,j}$'
        if make_2d_plots:
            img_ax = fig.add_subplot(gs[0:2, 2 * cut_num + 1])
            img = img_ax.imshow(avg_binned_corrs, vmin = correlation_color_lims[0], vmax = correlation_color_lims[1])
            img_colorbar = fig.colorbar(img, ax = img_ax)
            img_colorbar.set_label(correlation_str + ' in bin', rotation=270, labelpad=20)
    #img = img_ax.imshow(np.array(avg_binned_corrs) / np.array(avg_binned_corr_errs)) #, extent = [comoving_bins[0], comoving_bins[-1], comoving_bins[0], comoving_bins[-1]] )

        if make_2d_plots:
            n_corrs_ax = fig.add_subplot(gs[2:, 2 * cut_num + 1])
            n_corrs_img = n_corrs_ax.imshow(n_corrs_in_bin, norm=LogNorm())
            n_corrs_colorbar = fig.colorbar(n_corrs_img , ax = n_corrs_ax)
            n_corrs_colorbar.set_label(r'$N$ SNe pairs in bin', rotation=270, labelpad=20)

        comoving_ticks = [int(tick) for tick in np.linspace(1, len(comoving_bins)- 2, 7)]
        comoving_labels = [c.round_to_n(comoving_bins[tick] / 1000.0, 2) for tick in comoving_ticks]

        mean_z_ticks = [int(tick) for tick in np.linspace(1, len(z_bins) - 2, 7)]
        mean_z_labels = [c.round_to_n(z_bins[tick], 2) for tick in mean_z_ticks]

        comoving_perp_ticks = [c.round_to_n(tick, 2) for tick in np.linspace(r_perp_bounds[0] + (r_perp_bounds[1] - r_perp_bounds[0]) * 0.01, r_perp_bounds[1] - (r_perp_bounds[1] - r_perp_bounds[0]) * 0.01, 7)]
        comoving_perp_labels = [c.round_to_n(tick / 1000.0, 2) for tick in comoving_perp_ticks]

        comoving_par_ticks = [c.round_to_n(tick, 2) for tick in np.linspace(r_par_bounds[0] + (r_par_bounds[1] - r_par_bounds[0]) * 0.01, r_par_bounds[1] - (r_par_bounds[1] - r_par_bounds[0]) * 0.01, 7)]
        comoving_par_labels = [c.round_to_n(tick / 1000.0, 2) for tick in comoving_par_ticks]
        if make_2d_plots:
            img_ax.set_xticks(comoving_ticks)
            img_ax.set_yticks(mean_z_ticks)
            img_ax.set_xticklabels(comoving_labels)
            img_ax.set_yticklabels(mean_z_labels)
            img_ax.set_xlabel('Comoving separation between SNe (Gpc)')
            img_ax.set_ylabel('Mean redshift of SNe pairs')
            img_ax.set_xlim(-1, len (comoving_bins) - 1)
            img_ax.set_ylim(-1, len (z_bins) - 1)
            n_corrs_ax.set_xticks(comoving_ticks)
            n_corrs_ax.set_yticks(mean_z_ticks)
            n_corrs_ax.set_xticklabels(comoving_labels)
            n_corrs_ax.set_yticklabels(mean_z_labels)
            n_corrs_ax.set_xlabel('Comoving separation between SNe (Gpc)')
            n_corrs_ax.set_ylabel('Mean redshift of SNe pairs')
            n_corrs_ax.set_xlim(-1, len (comoving_bins) - 1)
            n_corrs_ax.set_ylim(-1, len (z_bins) - 1)

        comoving_scat_ax = fig.add_subplot(gs[0, 2 * cut_num])
        mean_z_scat_ax = fig.add_subplot(gs[1, 2 * cut_num])
        rperp_scat_ax = fig.add_subplot(gs[2, 2 * cut_num])
        rpar_scat_ax = fig.add_subplot(gs[3, 2 * cut_num])

        if make_1d_plots:
            print ('[np.shape(binned_zs), np.shape(binned_xis_on_zs[0])] = ' + str([np.shape(binned_zs), np.shape(binned_xis_on_zs[0])]))
            mean_z_scat_ax.scatter(binned_zs, binned_xis_on_zs[0], marker = '.', c = 'k')
            mean_z_scat_ax.errorbar(binned_zs, binned_xis_on_zs[0], fmt = 'none', yerr = binned_xis_on_zs[1], color = 'k')
            mean_z_x_lims = [mean_z_bounds[0] - (mean_z_bounds[1] - mean_z_bounds[0]) * 0.01, mean_z_bounds[1] + (mean_z_bounds[1] - mean_z_bounds[0]) * 0.01]
            mean_z_scat_ax.set_xlim(mean_z_x_lims)
            mean_z_scat_ax.plot(mean_z_x_lims, np.array(mean_z_x_lims) * 0.0, c = 'r')
            if not(ylims == None): mean_z_scat_ax.set_ylim(ylims)
            mean_z_scat_ax.set_xticks([float(tick_label) for tick_label in mean_z_labels])
            mean_z_scat_ax.set_xticklabels(mean_z_labels)
            mean_z_scat_ax.set_xlabel('Mean redshift of SNe pairs')
            mean_z_scat_ax.set_ylabel(correlation_str + ' in bin')

            comoving_scat_ax.scatter(binned_r_seps, binned_xis_on_r_seps[0], marker = '.', c = 'k')
            comoving_scat_ax.errorbar(binned_r_seps, binned_xis_on_r_seps[0], fmt = 'none', yerr = binned_xis_on_r_seps[1], color = 'k')
            comoving_x_lims = [r_sep_bounds[0] - (r_sep_bounds[1] - r_sep_bounds[0]) * 0.01, r_sep_bounds[1] + (r_sep_bounds[1] - r_sep_bounds[0]) * 0.01]
            comoving_scat_ax.set_xlim(comoving_x_lims)
            comoving_scat_ax.plot(comoving_x_lims, np.array(comoving_x_lims) * 0.0, c = 'r')
            if not(ylims == None): comoving_scat_ax.set_ylim(ylims)
            comoving_scat_ax.set_xticks([float(tick_label) * 1000.0 for tick_label in comoving_labels])
            comoving_scat_ax.set_xticklabels(comoving_labels)
            comoving_scat_ax.set_xlabel('Comoving separation between SNe (Gpc)')
            comoving_scat_ax.set_ylabel(correlation_str + ' in bin')

            rpar_scat_ax.scatter(binned_r_pars, binned_xis_on_r_pars[0], marker = '.', c = 'k')
            rpar_scat_ax.errorbar(binned_r_pars, binned_xis_on_r_pars[0], fmt = 'none', yerr = binned_xis_on_r_pars[1], color = 'k')
            rpar_x_lims = [r_par_bounds[0] - (r_par_bounds[1] - r_par_bounds[0]) * 0.01, r_par_bounds[1] + (r_par_bounds[1] - r_par_bounds[0]) * 0.01]
            rpar_scat_ax.set_xlim(rpar_x_lims)
            rpar_scat_ax.plot(rpar_x_lims, np.array(rpar_x_lims) * 0.0, c = 'r')
            if not(ylims == None): rpar_scat_ax.set_ylim(ylims)
            rpar_scat_ax.set_xticks(comoving_par_ticks)
            rpar_scat_ax.set_xticklabels(comoving_par_labels)
            rpar_scat_ax.set_xlabel('Comoving separation along line of sight (Gpc)')
            rpar_scat_ax.set_ylabel(correlation_str + ' in bin')

            rperp_scat_ax.scatter(binned_r_perps, binned_xis_on_r_perps[0], marker = '.', c = 'k')
            rperp_scat_ax.errorbar(binned_r_perps, binned_xis_on_r_perps[0], fmt = 'none', yerr = binned_xis_on_r_perps[1], color = 'k')
            rperp_x_lims = [r_perp_bounds[0] - (r_perp_bounds[1] - r_perp_bounds[0]) * 0.01, r_perp_bounds[1] + (r_perp_bounds[1] - r_perp_bounds[0]) * 0.01]
            rperp_scat_ax.set_xlim(rperp_x_lims)
            rperp_scat_ax.plot(rperp_x_lims, np.array(rperp_x_lims) * 0.0, c = 'r')
            if not(ylims == None): rperp_scat_ax.set_ylim(ylims)
            rperp_scat_ax.set_xticks(comoving_perp_ticks)
            rperp_scat_ax.set_xticklabels(comoving_perp_labels)
            rperp_scat_ax.set_xlabel('Comoving separation perpendicular to line of sight (Gpc)')
            rperp_scat_ax.set_ylabel(correlation_str + ' in bin')

    plt.tight_layout()

    if save_fig_name == None:
        save_fig_prefix = 'DeltaMuCorrelations_' + ''
        excised_save_str = '_'.join(surveys_to_excise)
        included_save_str = '_'.join(surveys_of_interest)
        zHD_str = ('zHD' if zHD else 'zCMB')
        ang_inclusion_str = str(c.round_to_n(ang_sep_cut[0], 3)) + 'To' + str(c.round_to_n(ang_sep_cut[1], 3))
        bin_save_strs = str(img_n_pixels[0]) + 'x' + str(img_n_pixels[0])
        save_fig_name = save_fig_prefix + 'included' + included_save_str + '_excluded' + excised_save_str + '_zHD' + zHD_str + '_angInclusion' + ang_inclusion_str + '_imgBins' + bin_save_strs + '.pdf'
    if save_fig:
        plt.savefig(save_dir + save_fig_name)
    if show_fig:
        plt.show()
    plt.close('all')

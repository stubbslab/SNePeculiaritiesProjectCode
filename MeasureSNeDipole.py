import cantrips as c
import numpy as np
import loadSN as lsn
import CosmologicalParameterArchive as cpa
import AstronomicalParameterArchive as apa
import scipy
from matplotlib.colors import LogNorm
from logList import logList
import samplingFunctions as sf
import VisitedDsphMCMCPoint as vdp
import matplotlib.pyplot as plt
import math
import randomSimulationFunctions as rSimF
import random
import itertools
import time

def readIn_snData(normal_randomize = 0, rearrange_randomize = 0,
                  OmM = 'default', OmL = 'default', H0 = 'default', OmR = 'default', Om0 = 'default', zHD = 1, surveys_of_interest = ['all'], surveys_to_excise = [], pull_extinctions = 0, z_include_range = [-np.inf, np.inf], remove_residuals = 1):
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
        all_nSigs = np.array(all_muDiffs) / np.array(all_muErrs)
        ordered_indeces = list(range(0, n_sn))
        random_z_indeces = ordered_indeces[:]
        random.shuffle(random_z_indeces)
        random_sky_indeces = ordered_indeces[:]
        random.shuffle(random_sky_indeces)
        rand_zs = [all_zs[index] for index in random_z_indeces]
        rand_ras = [all_ras[index] for index in random_sky_indeces]
        rand_decs = [all_decs[index] for index in random_sky_indeces]
        #rand_zs, rand_muDiffs, rand_muErrs  = rSimF.randomSortNSigma(all_zs, all_muDiffs, all_muErrs, [0.0 for mu in all_mus])
        all_zs = rand_zs
        all_ras = rand_ras
        all_decs = rand_decs

    all_muPreds = [sn['muPred'] for sn in all_sn]

    if remove_residuals:
        weights = 1.0 / (np.array(all_muErrs) ** 2.0)
        mu_offset = (np.sum(np.array(all_muDiffs) * np.array(weights)) / np.sum(weights) )
        print ('mu_offset = ' + str(mu_offset))
        all_muDiffs = [diff - mu_offset for diff in all_muDiffs]

    all_dLDiffs = [10.0 ** (muDiff / 5) - 1 for muDiff in all_muDiffs]
    all_dLDiffErrs = [all_muErrs[i] * np.log(10) / 5 * abs((all_dLDiffs[i] + 1)) for i in range(len(all_muDiffs))]

    param_to_correlate = all_dLDiffs
    param_errs = all_dLDiffErrs

    autocorrelations = c.flattenListOfLists([[param_to_correlate[i] * param_to_correlate[j] for j in range(i+1, n_sn) ] for i in range(n_sn)] )
    autocorrelation_errs = c.flattenListOfLists([[ np.sqrt((param_to_correlate[i] * param_errs[j]) ** 2.0 + (param_to_correlate[j] * param_errs[i]) ** 2.0) for j in range(i+1, n_sn) ] for i in range(n_sn)] )
    mean_zs = c.flattenListOfLists([[(all_zs[i] + all_zs[j]) / 2.0 for j in range(i+1, n_sn) ] for i in range(n_sn)] )

    Hubble_funct = lambda zp: np.sqrt(OmM * (1.0 + zp) ** 3.0 + OmL + OmR * (1.0 + zp) ** 4.0 + (1.0 - OmM - OmR - OmL) * (1.0 + zp) ** 2.0 )
    dl_funct = lambda z: dl_scaling * (1.0 + z) * integrate.quad(lambda zp: 1.0 / Hubble_funct(zp) , 0.0, z)[0]

    return [all_ras, all_decs, all_zs, all_muDiffs, all_muErrs, autocorrelations, autocorrelation_errs, mean_zs]

def makeResidualPlotOnSky(OmM = 'default', OmL = 'default', H0 = 'default', OmR = 'default', Om0 = 'default', zHD = 1,
                           surveys_of_interest = ['all'], surveys_to_excise = [], z_include_range = [-np.inf, np.inf], remove_residuals = 1,
                           figsize = [10, 8]):
    all_ras, all_decs, all_zs, all_muDiffs, all_muErrs, all_xis, all_xiErrs = readIn_snData(OmM = OmM, OmL = OmL, H0 = H0, OmR = OmR, Om0 = Om0, zHD = zHD, surveys_of_interest = surveys_of_interest, surveys_to_excise = surveys_to_excise, z_include_range = z_include_range, remove_residuals = remove_residuals)

    c.plotStarsOnSky(all_ras, all_decs, color = all_muDiffs, figsize = figsize)

    return 1


class SNDipoleFitter:

    def getInitParams(self, start_index = 0):
        #init_dipole, init_monopole, init_fit_shape, init_z_scaling, init_sky_vec = [0.0, 0.0, 1.0, 1.0, [0.0, 0.0, 0.1]]
        init_dipole, init_monopole, init_fit_shape, init_scaling =  [0.0, 0.1, 1.0, 0.1]
        init_params = [ [-0.05, -0.05, 0.1, 0.1],
                        [0.05, -0.05, 0.1, 0.1],
                        [-0.05, 0.05, 0.1, 0.1],
                        [0.05, 0.05, 0.1, 0.1],
                        [-0.05, -0.05, 0.3, 0.3],
                        [0.05, -0.05, 0.3, 0.3],
                        [-0.05, 0.05, 0.3, 0.3],
                        [0.05, 0.05, 0.3, 0.3], ]
                        #[-0.1, -0.1, 0.5, 0.08],
                        #[0.1, -0.1, 2.0, 0.08],
                        #[-0.1, 0.1, 2.0, 0.08],
                        #[0.1, 0.1, 2.0, 0.08],
                        #[-0.1, -0.1, 2.0, 0.16],
                        #[0.1, -0.1, 2.0, 0.16],
                        #[-0.1, 0.1, 2.0, 0.16],
                        #[0.1, 0.1, 2.0,  0.16] ]
        init_dipole, init_monopole, init_mu_dip, init_mu_mono =  init_params[start_index]
        #return [ init_dipole, init_monopole, init_fit_shape, init_z_scaling, init_sky_vec]
        return [ init_dipole, init_monopole, init_mu_dip, init_mu_mono ]

    def varyParams1(self):
        #o_dipole_amp, o_monopole_amp, o_shape_param, o_z_scaling, o_sky_unit_vector = self.current_params
        o_dipole_amp, o_monopole_amp, o_shape_param, o_scaling = self.current_params
        n_dipole_amp = sf.constrainedGauss(o_dipole_amp, *self.dipole_samp_params)
        n_monopole_amp = sf.constrainedGauss(o_monopole_amp, *self.monopole_samp_params)
        n_shape_param = sf.constrainedGauss(o_shape_param, *self.shape_param_samp_params)
        n_scaling = sf.constraintedGauss(o_scaling, *self.z_scaling_samp_params) #[0.1, [0.1, 2.0]]
        #n_sky_unit_vector = sf.varyUnitVector(o_sky_unit_vector, *self.unit_vector_samp_params)
        #n_params = [n_dipole_amp, n_monopole_amp, n_shape_param, n_z_scaling, n_sky_unit_vector]
        n_params = [n_dipole_amp, n_monopole_amp, n_shape_param, n_scaling]

        return n_params

    def varyParams2(self, jump_multiplier = 1):
        #o_dipole_amp, o_monopole_amp, o_shape_param, o_z_scaling, o_sky_unit_vector = self.current_params
        o_dipole_amp, o_monopole_amp, o_mu_dip, o_mu_mono = self.current_params
        n_dipole_amp = sf.constrainedGauss(o_dipole_amp, self.dipole_samp_params[0] * jump_multiplier, self.dipole_samp_params[1])
        n_monopole_amp = sf.constrainedGauss(o_monopole_amp, self.monopole_samp_params[0] * jump_multiplier, self.monopole_samp_params[1])
        n_mu_dip = sf.constrainedGauss(o_mu_dip, self.z_mu_dip_samp_params[0] * jump_multiplier, self.z_mu_dip_samp_params[1])
        n_mu_mono = sf.constrainedGauss(o_mu_mono, self.z_mu_mono_samp_params[0] * jump_multiplier, self.z_mu_mono_samp_params[1])
        #n_z_scaling_rat = sf.elSurrogateGauss(o_z_scaling_rat_param, *self.z_scaling_rat_samp_params) #[0.1, [0.1, 2.0]]
        #n_sky_unit_vector = sf.varyUnitVector(o_sky_unit_vector, *self.unit_vector_samp_params)
        #n_params = [n_dipole_amp, n_monopole_amp, n_shape_param, n_z_scaling, n_sky_unit_vector]
        n_params = [n_dipole_amp, n_monopole_amp, n_mu_dip, n_mu_mono]

        return n_params

    #Some information about skewnorm from here: https://en.wikipedia.org/wiki/Skew_normal_distribution
    def radial_funct(self, zs_in, amp, mu, sig): #radial_funct(self, zs_in, amp, shape_param, scaling):
        #delta = shape_param / np.sqrt(1.0 + shape_param ** 2.0)
        #skewness = (4.0 - np.pi) / 2.0 * (delta * np.sqrt(2.0 / np.pi)) ** 3.0 / (1.0 - 2.0 * delta ** 2.0 / np.pi) ** 2.0
        #mu_z = np.sqrt(2.0 / np.pi) * delta

        #Does not seem to work well
        #approx_peak_val = mu_z - skewness * np.sqrt(1.0 - mu_z ** 2.0) / 2.0 - np.sign(shape_param) / 2.0 * np.exp(-2.0 * np.pi / np.abs(shape_param))
        #min_res = scipy.optimize.minimize_scalar(lambda xs: -scipy.stats.skewnorm.pdf(xs, shape_param, 0.0, scaling), max(zs_in))
        #max_curve_val = -min_res['fun']
        #return amp * scipy.stats.skewnorm.pdf(zs_in, shape_param, 0.0, scaling) / max_curve_val
        return amp * np.exp(-(zs_in - mu) ** 2.0 / (2.0 * sig ** 2.0))


    def radial_funct2(self, zs_in, amp, z_scale1, power1, z_scale2):
        #z_scale2 = z_scale1 * z_scale_rat
        term1 = (zs_in / z_scale1) ** power1
        term2 = np.exp(-zs_in / z_scale2)
        peak_loc = power1 * z_scale2
        peak_val = (peak_loc / z_scale1) ** power1 * np.exp(-peak_loc / z_scale2)
        return amp * term1  * term2  / peak_val

    def ang_funct(self, phi, d):
        return d * np.cos(phi)

    def full_fit_funct(self, zs_in, phi, m, d, mu_dip, mu_mono, n_dip_sig_left = 3.0, n_mono_sig_left = 3.0):

        sig_dip = (mu_dip - 0.0) / (n_dip_sig_left)
        sig_mono = (mu_mono - 0.0) / (n_mono_sig_left)
        dip_part = self.ang_funct(phi, d)
        mono_part = m
        radial_part_dip = self.radial_funct(zs_in, 1.0, mu_dip, sig_dip)
        radial_part_mono = self.radial_funct(zs_in, 1.0, mu_mono, sig_mono)
        pred_residuals = dip_part * radial_part_dip + mono_part * radial_part_mono
        #if self.remove_residuals:
        #    weights = 1.0 / (np.array(self.all_muErrs) ** 2.0)
        #    mu_offset = (np.sum(np.array(pred_residuals) * np.array(weights)) / np.sum(weights) )
        #    pred_residuals = [residual - mu_offset for residual in pred_residuals]
        return pred_residuals

    def full_fit_funct2(self, zs_in, phi, m, d, z_power, z_scale2):
        ang_part = self.ang_funct(phi, m, d)
        radial_part = self.radial_funct2(zs_in, 1.0, 1.0, z_power, z_scale2)
        pred_residuals = ang_part * radial_part
        if self.remove_residuals:
            weights = 1.0 / (np.array(self.all_muErrs) ** 2.0)
            mu_offset = (np.sum(np.array(pred_residuals) * np.array(weights)) / np.sum(weights) )
            pred_residuals = [residual - mu_offset for residual in pred_residuals]
        return pred_residuals


    #def calcChiSqrOfParams(self, dipole_amp, monopole_amp, shape_param, z_scaling, sky_unit_vector):
    def calcChiSqrOfParams(self, angle_seps, dipole_amp, monopole_amp, mu_dip, mu_mono):

        #deg_to_rad = self.astro_arch.getDegToRad()
        #target_sky_coord = [c.goodArctan(sky_unit_vector[0], sky_unit_vector[1]) / deg_to_rad, (np.pi/2-np.arccos(sky_unit_vector[2])) / deg_to_rad]

        #angle_seps = [ang_sep * deg_to_rad for ang_sep in [c.measureAngularSeparationOnSky(target_sky_coord, [self.all_ras[i], self.all_decs[i]], return_radian = 0) for i in range(self.n_sn)]]
        #val_to_minimize = lambda m, d, shape_param, z_scaling: np.sum(((self.full_fit_funct(np.array(all_zs), np.array(angle_seps), m, d, shape_param, z_scaling) - np.array(all_muDiffs)) / np.array(all_muErrs)) ** 2.0)

        pred_mu_diffs = self.full_fit_funct(np.array(self.all_zs), np.array(angle_seps), monopole_amp, dipole_amp, mu_dip, mu_mono)




        chi_sqr = np.sum(((pred_mu_diffs - np.array(self.all_muDiffs)) / np.array(self.all_muErrs)) ** 2.0)
        r_chi_sqr = chi_sqr / (self.n_sn - self.n_fit_params - 1)
        #print ('[monopole_amp, dipole_amp, z_power, z_scale2, chi_sqr] = ' + str([monopole_amp, dipole_amp, z_power, z_scale2, chi_sqr]))
        #plt.scatter(self.all_zs, pred_mu_diffs, marker = '.')
        #plt.scatter(self.all_zs, self.all_muDiffs, marker = '.', c  = 'r' )
        #plt.pause(1.0)
        return chi_sqr

    #def calcChiSqrOfParams(self, dipole_amp, monopole_amp, shape_param, z_scaling, sky_unit_vector):
    def calcChiSqrOfParamsOnCorrs(self, angle_seps, dipole_amp, monopole_amp, z_power, z_scale2):

        #deg_to_rad = self.astro_arch.getDegToRad()
        #target_sky_coord = [c.goodArctan(sky_unit_vector[0], sky_unit_vector[1]) / deg_to_rad, (np.pi/2-np.arccos(sky_unit_vector[2])) / deg_to_rad]

        #angle_seps = [ang_sep * deg_to_rad for ang_sep in [c.measureAngularSeparationOnSky(target_sky_coord, [self.all_ras[i], self.all_decs[i]], return_radian = 0) for i in range(self.n_sn)]]
        #val_to_minimize = lambda m, d, shape_param, z_scaling: np.sum(((self.full_fit_funct(np.array(all_zs), np.array(angle_seps), m, d, shape_param, z_scaling) - np.array(all_muDiffs)) / np.array(all_muErrs)) ** 2.0)

        pred_mu_diffs = self.full_fit_funct2(np.array(self.all_zs), np.array(angle_seps), monopole_amp, dipole_amp, z_power, z_scale2)

        pred_dL_diffs = [10.0 ** (muDiff / 5) - 1 for muDiff in self.pred_mu_diffs]

        pred_xis = c.flattenListOfLists([[pred_dL_diffs[i] * pred_dL_diffs[j] for j in range(i+1, self.n_sn) ] for i in range(self.n_sn)] )

        chi_sqr = np.sum(((pred_xis - np.array(self.all_xis)) / np.array(self.all_xiErrs)) ** 2.0)
        r_chi_sqr = chi_sqr / (len(self.all_xis) - self.n_fit_params - 1)
        return chi_sqr

    #def compressed_fit_funct(self, angle_seps, param_vec):
    #    print('param_vec = ' + str(param_vec))
    #    return calcChiSqrOfParams(angle_seps, *param_vec)

    def measureSNDipoleViaGrid(self, spot_index ):
            start = time.time()
            print ('Working on spot_index = ' + str(spot_index) + ' of ' + str(self.n_sky_spots))
            spot_on_sky = self.spots_on_sky[spot_index]
            deg_to_rad = self.astro_arch.getDegToRad()
            target_sky_coord = [c.goodArctan(spot_on_sky[0], spot_on_sky[1]) / deg_to_rad, (np.pi/2-np.arccos(spot_on_sky[2])) / deg_to_rad]
            #print ('target_sky_coord = ' + str(target_sky_coord))
            current_angle_seps = [ang_sep * deg_to_rad for ang_sep in [c.measureAngularSeparationOnSky(target_sky_coord, [self.all_ras[i], self.all_decs[i]], return_radian = 0) for i in range(self.n_sn)]]
            mono_steps = np.arange(self.monopole_samp_params[1][0], self.monopole_samp_params[1][1] + self.monopole_samp_params[0] / 2.0, self.monopole_samp_params[0])
            dip_steps = np.arange(self.dipole_samp_params[1][0], self.dipole_samp_params[1][1] + self.dipole_samp_params[0] / 2.0, self.dipole_samp_params[0])
            power_steps = np.arange(self.z_power_samp_params[1][0], self.z_power_samp_params[1][1] + self.z_power_samp_params[0] / 2.0, self.z_power_samp_params[0])
            scaling_steps = np.arange(self.z_scaling_samp_params[1][0], self.z_scaling_samp_params[1][1] + self.z_scaling_samp_params[0] / 2.0, self.z_scaling_samp_params[0])

            #print ('[mono_steps, dip_steps, power_steps, scaling_steps] = ' + str([mono_steps, dip_steps, power_steps, scaling_steps]))


            #vec_fit_funct = np.vectorize(self.compressed_fit_funct, excluded=['angle_seps'])
            outer_product_of_params = list(itertools.product(dip_steps, mono_steps, power_steps, scaling_steps))

            #print ('outer_product_of_params = ' + str(outer_product_of_params))
            chi_sqrs_array = [self.calcChiSqrOfParams(current_angle_seps, *param_vec) for param_vec in outer_product_of_params]
            end = time.time()
            print ('Took ' + str(end - start) + 's')

            return [outer_product_of_params, chi_sqrs_array]

    def measureBestFitDipoleViaGrid(self):
        self.best_fit_params = [[] for i in range(self.n_sky_spots)]
        self.best_fit_chi_sqrs = [0.0 for i in range(self.n_sky_spots)]

        for spot_index in range(self.n_sky_spots):
            new_param_array, new_chi_sqrs_array = self.measureSNDipoleViaGrid( spot_index )
            best_fit_index = np.argmin(new_chi_sqrs_array)
            self.best_fit_chi_sqrs[spot_index] = new_chi_sqrs_array[best_fit_index ]
            self.best_fit_params[spot_index] = new_param_array[best_fit_index ]
            print ('[self.best_fit_params[spot_index], self.best_fit_chi_sqrs[spot_index]] = ' + str([self.best_fit_params[spot_index], self.best_fit_chi_sqrs[spot_index]]))

        return 1


    def measureSNDipoleViaMCMCOnCorrs(self, n_iters, init_index, long_jump_num = 50):

        for spot_index in range(self.n_sky_spots):
            print ('Working on spot_index = ' + str(spot_index) + ' of ' + str(self.n_sky_spots))
            spot_on_sky = self.spots_on_sky[spot_index]
            deg_to_rad = self.astro_arch.getDegToRad()
            target_sky_coord = [c.goodArctan(spot_on_sky[0], spot_on_sky[1]) / deg_to_rad, (np.pi/2-np.arccos(spot_on_sky[2])) / deg_to_rad]
            #print ('target_sky_coord = ' + str(target_sky_coord))
            current_angle_seps = [ang_sep * deg_to_rad for ang_sep in [c.measureAngularSeparationOnSky(target_sky_coord, [self.all_ras[i], self.all_decs[i]], return_radian = 0) for i in range(self.n_sn)]]
            init_params = self.getInitParams(init_index)
            self.current_params = init_params
            self.current_chi_sqr = self.calcChiSqrOfParams(current_angle_seps, *init_params)
            visited_points = [vdp.VisitedDsphMCMCPoint(self.current_params, self.current_chi_sqr, 1)]

            for i in range(n_iters):
                jump_multiplier = (100 if  i % long_jump_num == long_jump_num -1 else 1)
                new_params = self.varyParams2(jump_multiplier = jump_multiplier)
                #print ('new_params = ' + str(new_params))
                new_chi_sqr = self.calcChiSqrOfParamsOnCorrs(current_angle_seps, *new_params )
                #supposedly a good way to use chi squares with likelihoods is to look at the exponential of the chi square value
                step_prob = np.random.random()
                #print ('[self.current_params, self.current_chi_sqr] = ' + str([self.current_params, self.current_chi_sqr]))
                #print ('[new_params, new_chi_sqr] = ' + str([new_params, new_chi_sqr]))
                if i % 5000 == 5000 - 1:
                    print ('Just checked params ' + str(new_params) + ' with [new_chi_sqr, self.current_chi_sqr, np.exp(-new_chi_sqr + self.current_chi_sqr)] ' + str([new_chi_sqr, self.current_chi_sqr, np.exp(-new_chi_sqr + self.current_chi_sqr)]) + ' on iteration ' + str(i) + ' of ' + str(n_iters))
                if np.exp(-new_chi_sqr + self.current_chi_sqr ) > step_prob:
                    #self.current_params = [n_dipole_amp, n_monopole_amp, n_shape_param, n_z_scaling, n_sky_unit_vector]
                    self.current_params = new_params[:]
                    self.current_chi_sqr = new_chi_sqr
                    visited_points = visited_points + [vdp.VisitedDsphMCMCPoint(new_params[:], new_chi_sqr, 1)]
                else:
                    visited_points[-1].n_visits = visited_points[-1].n_visits + 1

            self.visited_points_dict[spot_index] = self.visited_points_dict[spot_index] + [visited_points[:]]

        return 1


    def measureSNDipoleViaMCMC(self, n_iters, init_index, long_jump_num = 50):

        for spot_index in range(self.n_sky_spots):
            print ('Working on spot_index = ' + str(spot_index) + ' of ' + str(self.n_sky_spots))
            spot_on_sky = self.spots_on_sky[spot_index]
            deg_to_rad = self.astro_arch.getDegToRad()
            target_sky_coord = [c.goodArctan(spot_on_sky[0], spot_on_sky[1]) / deg_to_rad, (np.pi/2-np.arccos(spot_on_sky[2])) / deg_to_rad]
            #print ('target_sky_coord = ' + str(target_sky_coord))
            current_angle_seps = [ang_sep * deg_to_rad for ang_sep in [c.measureAngularSeparationOnSky(target_sky_coord, [self.all_ras[i], self.all_decs[i]], return_radian = 0) for i in range(self.n_sn)]]
            init_params = self.getInitParams(init_index)
            self.current_params = init_params
            self.current_chi_sqr = self.calcChiSqrOfParams(current_angle_seps, *init_params)
            visited_points = [vdp.VisitedDsphMCMCPoint(self.current_params, self.current_chi_sqr, 1)]

            for i in range(n_iters):
                jump_multiplier = (10 if  i % long_jump_num == long_jump_num -1 else 1)
                new_params = self.varyParams2(jump_multiplier = jump_multiplier)
                #print ('new_params = ' + str(new_params))
                new_chi_sqr = self.calcChiSqrOfParams(current_angle_seps, *new_params )
                #supposedly a good way to use chi squares with likelihoods is to look at the exponential of the chi square value
                step_prob = np.random.random()
                #print ('[self.current_params, self.current_chi_sqr] = ' + str([self.current_params, self.current_chi_sqr]))
                #print ('[new_params, new_chi_sqr] = ' + str([new_params, new_chi_sqr]))
                if i % 5000 == 5000 - 1:
                    print ('Just checked params ' + str(new_params) + ' with [new_chi_sqr, self.current_chi_sqr, np.exp(-new_chi_sqr + self.current_chi_sqr)] ' + str([new_chi_sqr, self.current_chi_sqr, np.exp(-new_chi_sqr + self.current_chi_sqr)]) + ' on iteration ' + str(i) + ' of ' + str(n_iters))
                if np.exp(-new_chi_sqr + self.current_chi_sqr ) > step_prob:
                    #self.current_params = [n_dipole_amp, n_monopole_amp, n_shape_param, n_z_scaling, n_sky_unit_vector]
                    self.current_params = new_params[:]
                    self.current_chi_sqr = new_chi_sqr
                    visited_points = visited_points + [vdp.VisitedDsphMCMCPoint(new_params[:], new_chi_sqr, 1)]
                else:
                    visited_points[-1].n_visits = visited_points[-1].n_visits + 1

            self.visited_points_dict[spot_index] = self.visited_points_dict[spot_index] + [visited_points[:]]
        return 1

    def zScalingProxyFunct(self, z_scaling):
        z_scaling = np.array(z_scaling)
        new_z_scaling = np.where(z_scaling > 1.0, z_scaling, 2.0 - 1.0 / z_scaling )
        return new_z_scaling

    def zScalingInvProxyFunct(self, z_scaling_proxy):
        z_scaling_proxy = np.array(z_scaling_proxy)
        z_scaling = np.where(z_scaling_proxy > 1.0, z_scaling_proxy,  1.0 / (2.0 - z_scaling_proxy ) )
        return z_scaling

    def zPowerProxyFunct(self, z_scaling):
        return self.zScalingProxyFunct(z_scaling)

    def __init__(self, n_sky_spots = 100, normal_randomize = 0, rearrange_randomize = 0, OmM = 'default', OmL = 'default', H0 = 'default', OmR = 'default', Om0 = 'default', zHD = 1, surveys_of_interest = ['all'], surveys_to_excise = [], z_include_range = [-np.inf, np.inf], remove_residuals = 1, added_monopole_params = [0.0, 1.0, 0.1], added_dipole_params = [0.0, 0.0 * 180.0, 1.0 * 180.0,  0.1, 1.0],):


        self.all_ras, self.all_decs, self.all_zs, self.all_muDiffs, self.all_muErrs, self.all_xis, self.all_xiErrs, self.mean_zs  = readIn_snData(normal_randomize = normal_randomize, rearrange_randomize =  rearrange_randomize, OmM = OmM, OmL = OmL, H0 = H0, OmR = OmR, Om0 = Om0, zHD = zHD, surveys_of_interest = surveys_of_interest, surveys_to_excise = surveys_to_excise, z_include_range = z_include_range, remove_residuals = remove_residuals)
        self.n_sn = len(self.all_ras)
        self.n_fit_params = 6 #dipole, monopole, fit_shape, z_scaling, target_RA, target_Dec

        self.randomize = (normal_randomize or rearrange_randomize)
        self.remove_residuals = remove_residuals
        if added_dipole_params[0] != 0.0:
            artificial_ang_seps = [c.measureAngularSeparationOnSky([added_dipole_params[1], added_dipole_params[2]], [self.all_ras[i], self.all_decs[i]], return_radian = 0) for i in range(len(all_ras))]
            self.all_muDiffs = self.all_muDiffs + self.ang_funct(artificial_ang_seps, 0.0, added_dipole_params[0]) * self.radial_funct(np.array(self.all_zs), added_dipole_params[3], added_dipole_params[4])


        self.all_muDiffs = self.all_muDiffs + self.radial_funct(np.array(self.all_zs), added_monopole_params[0], added_monopole_params[1], added_monopole_params[2])
        self.null_chi_sqr = np.sum((np.array(self.all_muDiffs) / np.array(self.all_muErrs)) ** 2.0)
        all_dLDiffs = [10.0 ** (muDiff / 5) - 1 for muDiff in self.all_muDiffs]
        all_dLDiffErrs = [self.all_muErrs[i] * np.log(10) / 5 * abs((all_dLDiffs[i] + 1)) for i in range(len(self.all_muDiffs))]

        param_to_correlate = all_dLDiffs
        param_errs = all_dLDiffErrs
        self.all_xis = c.flattenListOfLists([[param_to_correlate[i] * param_to_correlate[j] for j in range(i+1, self.n_sn) ] for i in range(self.n_sn)] )
        self.all_xiErrs = c.flattenListOfLists([[ np.sqrt((param_to_correlate[i] * param_errs[j]) ** 2.0 + (param_to_correlate[j] * param_errs[i]) ** 2.0) for j in range(i+1, self.n_sn) ] for i in range(self.n_sn)] )

        #self.all_muDiffs = np.array([0.0 for mu in self.all_muDiffs]) +  self.radial_funct2(np.array(self.all_zs), added_monopole_params[0], 1.0, added_monopole_params[1], added_monopole_params[2])
        print ('self.calcChiSqrOfParams([0.0 for mu in self.all_muDiffs], 0.0, added_monopole_params[0], added_monopole_params[1], added_monopole_params[2]) ' + str(self.calcChiSqrOfParams([0.0 for mu in self.all_muDiffs], 0.0, added_monopole_params[0], added_monopole_params[1], added_monopole_params[2]) ))

        self.unit_vector_samp_params = [0.02]
        self.dipole_samp_params = [0.01, [-0.2, 0.2]] # [0.01, [-0.0, 0.0]]
        self.monopole_samp_params = [0.01, [-0.2, 0.2]]
        self.shape_param_samp_params = [0.5, [0.0, 5.0]]
        #self.z_scaling_samp_params = [1.0, [0.0005, 0.2]]
        self.z_mu_dip_samp_params = [0.025, [0.01, 0.4]] #[0.5, [0.2, 5.0]]
        self.z_mu_mono_samp_params = [0.01, [0.01, 0.4]] #Minimum z is ~0.01, so make that our minimum exponential scaling
        self.z_power_samp_params = [0.1, [0.01, 10.0]]
        self.z_scaling_samp_params = [0.001, [0.001, 0.05]]
        self.z_scaling_rat_samp_params = [0.1, [0.1, 10.0]]

        self.param_order_indeces = {'dip':0, 'mono':1, 'mu_dip':2, 'mu_mono':3}
        self.default_param_vals = {'dip':0.0, 'mono':0.0, 'mu_dip':0.2, 'mu_mono':0.2}

        self.astro_arch = apa.AstronomicalParameterArchive()

        #self.spots_on_sky = [spot for spot in c.getPointsOnSphere(n_sky_spots * 2) if spot[2] > 0.0]
        self.spots_on_sky = c.getPointsOnSphere(n_sky_spots )
        self.n_sky_spots = len(self.spots_on_sky)
        self.visited_points_dict = {spot_index:[] for spot_index in range(self.n_sky_spots)}

        self.gaussianPosteriorFits = {spot_index:{} for spot_index in range(self.n_sky_spots)}


    def getParamDisplayVals(self):
        n_bins_array = {'mono':100, 'dip':100, 'shape_param':100, 'power':100, 'scaling':100, 'scaling_rat':100, 'mu_dip':100, 'mu_mono':100}
        param_display_vals = {'mono':[*self.monopole_samp_params[1], n_bins_array['mono']],
                              'dip':[*self.dipole_samp_params[1],n_bins_array['dip']],
                              'shape_param':[*self.shape_param_samp_params[1], n_bins_array['shape_param']],
                              #'power':[*self.zPowerProxyFunct(self.z_power_samp_params[1]), n_bins_array['power']],
                              'power':[*self.z_power_samp_params[1], n_bins_array['power']],
                              'scaling':[*self.z_scaling_samp_params[1], n_bins_array['scaling']],
                              #'scaling_rat':[*self.zScalingProxyFunct(self.z_scaling_rat_samp_params[1]), n_bins_array['scaling_rat']]
                              'scaling_rat':[*self.z_scaling_rat_samp_params[1], n_bins_array['scaling_rat']],
                              'mu_dip':[*self.z_mu_dip_samp_params[1], n_bins_array['mu_dip']],
                              'mu_mono':[*self.z_mu_mono_samp_params[1], n_bins_array['mu_mono']],
                              }
        return [n_bins_array, param_display_vals]

    def getDisplayUnits(self):
        display_labels_dict = {'mono':'mag', 'dip':'mag', 'shape_param':None, 'power':None, 'scaling':None, 'scaling_rat':None, 'mu_dip':None, 'mu_mono':None}
        return display_labels_dict

    def getDisplayLabels(self):
        display_labels_dict = {'mono':'mag', 'dip':'mag', 'shape_param':None, 'power':None, 'scaling':None, 'scaling_rat':None, 'mu_dip':None, 'mu_mono':None}
        return display_labels_dict

    def getDisplayLabels(self):
        disp_labels = {'mono':r'$m$ (mag)', 'dip':r'$d$ (mag)', 'shape_param':r'$\alpha$', 'power':r'$\gamma$', 'scaling':r'$z_0$', 'scaling_rat':r'$z_2/z_1$', 'mu_dip':r'$\mu_{D}$', 'mu_mono':r'$\mu_{M}$'}
        return disp_labels

    def showChiSqrOnSky(self, min_via_mcmc = 1, plot_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/SNIsotropyProject/plots/', file_name = 'ChiSquareOfDipoleOnSky.pdf', cbar_label = r'Best $\chi^2$ of dipole fit', show_fig = 1, save_fig = 0, funct_params = ['dip', 'mono', 'mu_dip', 'mu_mono']):
        deg_to_rad = self.astro_arch.getDegToRad()
        target_sky_coords = [[c.goodArctan(spot_on_sky[0], spot_on_sky[1]) / deg_to_rad, (np.pi/2-np.arccos(spot_on_sky[2])) / deg_to_rad] for spot_on_sky in self.spots_on_sky]
        gaussian_fits = [self.gaussianPosteriorFits[sky_coord_index] for sky_coord_index in range(len(target_sky_coords))]
        if min_via_mcmc:
            fit_params = [[(gaussian_fit[funct_param][1] if funct_param in gaussian_fit.keys() else self.default_param_vals[funct_param]) for funct_param in funct_params] for gaussian_fit in gaussian_fits]
            #print ('fit_params = ' + str(fit_params))

            chi_sqrs = [0.0 for target_sky_coord in target_sky_coords]
            for i in range(len(target_sky_coords)):
                target_sky_coord = target_sky_coords[i]
                print ('target_sky_coord = ' + str(target_sky_coord))
                angle_seps = [c.measureAngularSeparationOnSky(target_sky_coord, [self.all_ras[i], self.all_decs[i]], return_radian = 0) * deg_to_rad for i in range(self.n_sn)]
                print("[fit_params[i][funct_params.index('dip')], fit_params[i][funct_params.index('mono')], self.zScalingInvProxyFunct(fit_params[i][funct_params.index('mu_dip')]), fit_params[i][funct_params.index('mu_mono')]] = " + str([fit_params[i][funct_params.index('dip')], fit_params[i][funct_params.index('mono')], self.zScalingInvProxyFunct(fit_params[i][funct_params.index('mu_dip')]), fit_params[i][funct_params.index('mu_mono')]]))
                best_fit_params = [fit_params[i][funct_params.index(param_label)] for param_label in funct_params ]
                chi_sqrs[i] = self.calcChiSqrOfParams(angle_seps,  *best_fit_params)
                print ('chi_sqrs[i] = ' + str(chi_sqrs[i]))
        else:
            chi_sqrs = [self.best_fit_chi_sqrs[spot_index] for spot_index in range(self.n_sky_spots)]

        c.plotStarsOnSky([target_sky_coord[0] for target_sky_coord in target_sky_coords], [target_sky_coord[1] for target_sky_coord in target_sky_coords], color = [chi_sqr if chi_sqr < self.null_chi_sqr else self.null_chi_sqr for chi_sqr in chi_sqrs], add_colorbar = 1, cbar_label = cbar_label,
                         save_fig = save_fig, show_fig = show_fig, plot_dir = plot_dir, file_name = file_name)
        print ('chi_sqrs = ' + str(chi_sqrs))

        return chi_sqrs



    def showMCMCResults(self, spot_index_to_disp, vars_to_disp = ['dip', 'mono', 'mu_dip', 'mu_mono'],
                        params_to_fit = {}, fit_funct = 'poly', fitting_funct = lambda mcmc_param, amp, mu, sig, shift, slope: amp * np.exp(-(np.array(mcmc_param) - mu) ** 2.0 / (2.0 * sig ** 2.0)) + np.array(mcmc_param) * slope + shift,
                        functions_to_display = [], show = 1, save_fig = 0, figpadding = 6.0,
                        results_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/randall/MCMCOutputs/',
                        plot_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/SNIsotropyProject/plots/',
                        file_name = '', n_x_ticks = 5, n_y_ticks = 5, n_fit_bin_buffer = 2,
                        n_bins = [], n_ignore = 0, theta_shift = 0.0, n_var_per_plot = 'both', fig_size_unit = 2.0, n_levels = 10, exp_level_powers = [-3.0, -2.0, -1.0],
                        existing_subplots_array = None, return_axarr = 0, fit_line_color = 'black', n_test_lines = 100.0,
                        n_fitted_points = 100, n_xy_points_to_fit = 1000, smallest_max_val_to_be_on_curve = 50.0, fancyTex = 0):

        plt.rc('font', family='serif')
        if fancyTex:
            plt.rc('text', usetex=True)
        #plt.rc('xtick', labelsize='x-small')
        #plt.rc('ytick', labelsize='x-small')
        #var_cyclicity = getVarCyclicity()
        #var_display_scalings = getVarDisplayScalings()
        #bin_types = getBinTypes()
        #var_ticks = getVarTicks()
        #var_ticklabels = getVarTickLabels()

        spot_on_sky = self.spots_on_sky[spot_index_to_disp]
        deg_to_rad = self.astro_arch.getDegToRad()
        target_sky_coord = [c.goodArctan(spot_on_sky[0], spot_on_sky[1]) / deg_to_rad, (np.pi/2-np.arccos(spot_on_sky[2])) / deg_to_rad]

        angle_seps = [ang_sep * deg_to_rad for ang_sep in [c.measureAngularSeparationOnSky(target_sky_coord, [self.all_ras[i], self.all_decs[i]], return_radian = 0) for i in range(self.n_sn)]]

        disp_units = self.getDisplayUnits()
        display_labels = self.getDisplayLabels()

        vars_to_read_in = vars_to_disp

        n_bins_array , param_display_vals = self.getParamDisplayVals()

        #measured_arrays, n_visits_array, likelihood_array = readInMCMCResults(results_files, functions_to_compute = functions_to_display, vars_to_load = vars_to_read_in,
        #                                    results_dir = results_dir, n_ignore = n_ignore )
        visited_points_to_display = self.visited_points_dict[spot_index_to_disp]
        n_chains = len(visited_points_to_display )
        measured_arrays = []
        n_visits_array = []
        chi_sqr_array = []
        for chain_num in range(n_chains):
            total_chain_visits = 0
            visited_points = visited_points_to_display[chain_num]
            single_chain_n_visits = [point.n_visits for point in visited_points]
            for n_visited_index in range(len(single_chain_n_visits)):
                n_visited =  single_chain_n_visits[n_visited_index]
                total_chain_visits = total_chain_visits + n_visited
                if total_chain_visits > n_ignore:
                    cutoff_visited_index = n_visited_index
                    break
            points_to_include = visited_points[cutoff_visited_index:]
            measured_arrays = measured_arrays + [point.parameter_storer for point in points_to_include ]
            n_visits_array = n_visits_array + [point.n_visits for point in points_to_include ]
            chi_sqr_array = chi_sqr_array + [point.log_likelihood for point in points_to_include ]

        print ('[len(measured_arrays), len(n_visits_array), len(chi_sqr_array)] = ' + str([len(measured_arrays), len(n_visits_array), len(chi_sqr_array)] ))

        #measured_arrays, n_visits_array, chi_sqr_array = [[point.parameter_storer for point in visited_points_to_display], [point.n_visits for point in visited_points_to_display], [point.log_likelihood for point in visited_points_to_display]]
        measured_arrays = [[measure_array[i] for measure_array in measured_arrays] for i in range(len(measured_arrays[0]))]
        #print ('measured_arrays = ' + str(measured_arrays))

        total_visits = sum(n_visits_array)

        for var in vars_to_disp:
            var_min = param_display_vals[var][0]
            var_max = param_display_vals[var][1]
            var_nbins = param_display_vals[var][2]
            var_step = (var_max - var_min) / var_nbins
            var_bins = [[var_min + var_step*i, var_min + var_step *(i+1)] for i in range(var_nbins)]

            var_bin_centers = [(var_bins[i][1] + var_bins[i][0])/ 2.0 for i in range(len(var_bins)) ]
            print ('[var, var_min, var_max, var_nbins, var_bin_centers] = ' + str([var, var_min, var_max, var_nbins, var_bin_centers] ))
            param_display_vals[var] = param_display_vals[var] + [var_step, var_bins, var_bin_centers]


        #plt.rcParams.update({'font.size': 6})

        if n_var_per_plot == 2 or n_var_per_plot in ['both','Both','BOTH']:

            n_x_contour_plots = len(vars_to_disp) - 1
            n_y_contour_plots = len(vars_to_disp) - 1
            if n_var_per_plot == 2:
                total_n_x_plots = n_x_contour_plots
                total_n_y_plots = n_y_contour_plots
            else:
                total_n_x_plots = n_x_contour_plots + 1
                total_n_y_plots = n_y_contour_plots + 1
            if existing_subplots_array is None:
                f, axarr = plt.subplots(total_n_x_plots, total_n_y_plots,
                                        figsize = (fig_size_unit * total_n_x_plots, fig_size_unit * total_n_y_plots),
                                        squeeze = False) #, sharex = True, sharey = True)
                plt.tight_layout(pad = figpadding )
                #for i in range(np.shape(axarr)[0]):
                #    for j in range(np.shape(axarr)[1]):
                #        axarr[j][i].set_xticks([])
                #        axarr[j][i].set_yticks([])
                #        axarr[j][i].tick_params(axis = 'x', direction = 'in')
                #        axarr[j][i].tick_params(axis = 'y', direction = 'in')
                #plt.show()

            total_n_visits = 0
            for i in range(len(vars_to_disp) - 1):
                #for j in [elem + 1 for elem in range(i+1)]:
                for j in range(i, len(vars_to_disp) - 1):
                    #print 'i = ' + str(i)
                    #print 'j = ' + str(j)
                    if n_var_per_plot is 2:
                        axarr_x_index = i
                        #axarr_y_index = j-1
                        axarr_y_index = j
                    else:
                        axarr_x_index = i + 1
                        axarr_y_index = j
                    var_x_index = i
                    var_y_index = j+1
                    print ('[axarr_y_index, axarr_x_index] = ' + str([axarr_y_index, axarr_x_index]))
                    x_var = vars_to_disp[var_x_index]
                    y_var = vars_to_disp[var_y_index]

                    print ('[x_var, y_var] = ' + str([x_var, y_var]))
                    x_bin_centers = param_display_vals[x_var][5]
                    y_bin_centers = param_display_vals[y_var][5]

                    x_mesh, y_mesh = np.meshgrid(x_bin_centers, y_bin_centers)

                    binned_n_visits = np.zeros(np.shape(x_mesh))

                    x_bins = param_display_vals[x_var][4]
                    y_bins = param_display_vals[y_var][4]
                    #if y_var == 'scaling':
                    #    fig2, axarr2 = plt.subplots(1,1)
                    #    print ('measured_arrays[var_y_index] = ' + str(measured_arrays[var_y_index]))
                    #    axarr2.hist(measured_arrays[var_y_index], bins = 101)
                    #    plt.show()
                    for k in range(len(n_visits_array)):
                        measured_x = measured_arrays[self.param_order_indeces[x_var]][k]
                        #if x_var == 'power':
                        #    measured_x = self.zPowerProxyFunct(measured_x)
                        measured_y = measured_arrays[self.param_order_indeces[y_var]][k]
                        #if y_var == 'power':
                        #    measured_y = self.zPowerProxyFunct(measured_y)
                        n_visits = n_visits_array[k]
                        x_bin_index = -1
                        y_bin_index = -1

                            #if x_var == 'scaling_rat':
                            #    print ('[measured_x, x_bins] = ' + str([measured_x, x_bins] ))
                            #if y_var == 'scaling_rat':
                            #    print ('[measured_y, y_bins] = ' + str([measured_y, y_bins] ))
                        for l in range(len(x_bins)):
                            if measured_x >=x_bins[l][0] and measured_x < x_bins[l][1]:
                                x_bin_index = l
                        for l in range(len(y_bins)):
                            if measured_y >=y_bins[l][0] and measured_y < y_bins[l][1]:
                                y_bin_index = l
                        total_n_visits = total_n_visits + n_visits
                        #print ('total_n_visits = ' + str(total_n_visits))
                        if total_n_visits > n_ignore:
                            binned_n_visits[y_bin_index,x_bin_index] = binned_n_visits[y_bin_index,x_bin_index] + n_visits

                    max_visits = np.max(np.abs(binned_n_visits + 1.0))
                    log_levels = logList(1.0, max_visits + 1.0 , n_levels)
                    exp_levels = [max_visits * np.exp(power) for power in exp_level_powers]
                    levels = log_levels

                    xlabel = display_labels[x_var]
                    ylabel = display_labels[y_var]
                    #if bin_types[x_var] == 'log':
                    #    axarr[axarr_y_index][axarr_x_index].set_xscale('log')
                    #    axarr[axarr_y_index][axarr_x_index].tick_params(top = False, bottom = False, which = 'minor')
                    #    if (n_var_per_plot in ['both','Both','BOTH']):
                    #        axarr[-1][axarr_x_index].set_xscale('log')
                    #        axarr[-1][axarr_x_index].tick_params(top = False, bottom = False, which = 'minor')
                    #if bin_types[y_var] == 'log':
                    #    axarr[axarr_y_index][axarr_x_index].set_yscale('log')
                    #    axarr[axarr_y_index][axarr_x_index].tick_params(right = False, left= False, which = 'minor')
                    #    if (n_var_per_plot in ['both','Both','BOTH']):
                    #        axarr[axarr_y_index][0].set_yscale('log')
                    #        axarr[axarr_y_index][0].tick_params(right = False, left= False, which = 'minor')
                    print ('[axarr_y_index, axarr_x_index] ' + str([axarr_y_index, axarr_x_index]))
                    CS = axarr[axarr_y_index][axarr_x_index].contour(x_mesh, y_mesh, binned_n_visits + 1.0, levels = levels, norm = LogNorm())
                    #CS = axarr[axarr_y_index][axarr_x_index].imshow(binned_n_visits + 1.0, norm = LogNorm(), extent = [x_bin_centers[0], x_bin_centers[-1], y_bin_centers[0], y_bin_centers[-1]], origin = 'lower')

                    if (x_var in params_to_fit and y_var in params_to_fit):
                        peak_xs, peak_ys = fitMCMCResults(results_files, [x_var, y_var], measured_arrays = measured_arrays, n_visits_array = n_visits_array, param_ranges_to_fit = [ params_to_fit[x_var],params_to_fit[y_var]],
                                                          results_dir = results_dir, n_ignore = n_ignore, theta_shift = theta_shift, n_fitted_points = n_fitted_points, n_xy_points_to_fit = n_xy_points_to_fit, smallest_max_val_for_fit = smallest_max_val_to_be_on_curve)
                        print ('peak_xs = ' + str(peak_xs))
                        print ('peak_ys = ' + str(peak_ys))
                        peak_interpolator = interpolate.interp1d( *safeSortOneListByAnother(peak_xs, [peak_xs, peak_ys]), kind = 'linear')
                        #axarr[0][1].contour(x_mesh, y_mesh, binned_n_visits_interp((x_mesh, y_mesh)) + 1.0, 20, levels = levels, norm = LogNorm())

                        axarr[axarr_y_index][axarr_x_index].scatter(np.array(peak_xs) * x_scaling, np.array(peak_ys) * y_scaling, c = 'black', s = 2.0 * (scaled_x_bin_centers[1] - scaled_x_bin_centers[0]) )
                        axarr[axarr_y_index][axarr_x_index].plot([x for x in scaled_x_bin_centers if (x > min(peak_xs) and x < max(peak_xs))], [peak_interpolator(x / x_scaling) * y_scaling for x in scaled_x_bin_centers if (x > min(peak_xs) and x < max(peak_xs))], c = fit_line_color)

                    #axarr[i][j-1].colorbar(CS,format = '%.6f')
                    #xticks = var_ticks[x_var]
                    #xticklabels = var_ticklabels[x_var]
                    #yticks = var_ticks[y_var]
                    #yticklabels = var_ticklabels[y_var]
                    #axarr[axarr_y_index][axarr_x_index].set_xticks(xticks)
                    #axarr[axarr_y_index][axarr_x_index].set_yticks(yticks)
                    if (axarr_y_index == len(vars_to_disp) - 2 and n_var_per_plot is 2):
                        axarr[axarr_y_index][axarr_x_index].set_xlabel(xlabel, fontsize = 3.0 * (fig_size_unit) * n_y_contour_plots)
                    #    axarr[axarr_y_index][axarr_x_index].set_xticklabels(xticklabels, fontsize = 5.0 * fig_size_unit * n_y_contour_plots/ (n_y_ticks), rotation = 0)
                    #else:
                    #    axarr[axarr_y_index][axarr_x_index].set_xticklabels([])
                    if (axarr_x_index == 0 and n_var_per_plot is 2) :
                        axarr[axarr_y_index][axarr_x_index].set_ylabel(ylabel, fontsize = 2.0 * (fig_size_unit) * n_y_contour_plots)
                    #    axarr[axarr_y_index][axarr_x_index].set_yticklabels(yticklabels, fontsize = 10.0 * fig_size_unit * n_y_contour_plots/ (n_y_ticks), rotation = 0)
                    #else:
                    #    axarr[axarr_y_index][axarr_x_index].set_yticklabels([])

                    if n_var_per_plot is 'both':
                    #    axarr[-1][axarr_x_index].set_xticks(xticks)
                        axarr[-1][axarr_x_index].set_xlabel(xlabel, fontsize = 2.0 * (fig_size_unit) * n_y_contour_plots)
                    #    axarr[-1][axarr_x_index].set_xticklabels(xticklabels, fontsize = 10.0 * fig_size_unit * n_y_contour_plots/ (n_y_ticks), rotation = 0)
                    #    axarr[axarr_y_index][0].set_yticks(yticks)
                        axarr[axarr_y_index][0].set_ylabel(ylabel, fontsize = 2.0 * (fig_size_unit) * n_y_contour_plots)

                    #print ('[x_bin_centers, y_bin_centers] = ' + str([x_bin_centers, y_bin_centers] ))
                    axarr[axarr_y_index][axarr_x_index].set_xlim([min(x_bin_centers), max(x_bin_centers)])
                    axarr[axarr_y_index][axarr_x_index].set_ylim([min(y_bin_centers), max(y_bin_centers)])
                    if n_var_per_plot is 'both':
                        axarr[-1][axarr_x_index].set_xlim([min(x_bin_centers), max(x_bin_centers)])
                        axarr[axarr_y_index][0].set_ylim([min(y_bin_centers), max(y_bin_centers)])

            f.subplots_adjust(hspace = 0.0, wspace = 0.0)

        if n_var_per_plot == 1 or n_var_per_plot in ['both','Both','BOTH']:

            if n_var_per_plot == 1:
                n_x_plots = int( math.ceil( math.sqrt(len(vars_to_disp)) ) )
                n_y_plots = int( math.ceil( len(vars_to_disp) / float(n_x_plots) ) )
                f, axarr = plt.subplots(n_x_plots, n_y_plots, figsize = (fig_size_unit * n_x_plots, fig_size_unit * n_y_plots), squeeze = False)

            max_n_visits = 0
            best_fit_params = {'dip':0.0, 'mono':0.0, 'power':1.0, 'scaling':0.03}
            for i in range(len(vars_to_disp)):
                var = vars_to_disp[i]
                xlabel = display_labels[var]
                #x_bin_centers = x_bin_centers + [param_display_vals[var][5]]
                x_bin_centers = param_display_vals[var][5]
                binned_n_visits = np.zeros((len(x_bin_centers)))
                for k in range(len(n_visits_array)):
                    measured_x = measured_arrays[self.param_order_indeces[var]][k]
                    #if var == 'power':
                    #    measured_x = self.zPowerProxyFunct(measured_x)
                    x_bins = param_display_vals[var][4]
                    n_visits = n_visits_array[k]
                    x_bin_index = -1
                    for j in range(len(x_bins)):
                        if measured_x >=x_bins[j][0] and measured_x < x_bins[j][1]:
                            x_bin_index = j

                    binned_n_visits[x_bin_index] = binned_n_visits[x_bin_index] + n_visits

                x_index = i + 1
                y_index = i - 1
                ['dip','mono', 'power', 'scaling']
                fit_res = [0.0, 0.0, 1.0, 0.03, 1.0]
                if n_fit_bin_buffer == 0:
                    x_bin_centers_to_fit = x_bin_centers[:]
                    binned_n_visits_to_fit = binned_n_visits[:]
                else:
                    x_bin_centers_to_fit = x_bin_centers[n_fit_bin_buffer:-n_fit_bin_buffer]
                    binned_n_visits_to_fit = binned_n_visits[n_fit_bin_buffer:-n_fit_bin_buffer]

                max_index = np.argmax(binned_n_visits_to_fit)
                oneDFitGuess = [binned_n_visits_to_fit[max_index], x_bin_centers_to_fit[max_index], np.abs(x_bin_centers_to_fit[max_index] - x_bin_centers_to_fit[np.argmin(np.abs(np.array(binned_n_visits_to_fit) - binned_n_visits_to_fit[max_index] * 0.5))]), 0.0, 0.0]
                lower_bounds = [-np.inf, x_bin_centers_to_fit[0], x_bin_centers_to_fit[1] - x_bin_centers_to_fit[0], -np.inf, -np.inf]
                upper_bounds = [np.inf, x_bin_centers_to_fit[-1], (x_bin_centers_to_fit[-1] - x_bin_centers_to_fit[0])/ 2.0, np.inf, np.inf]
                for bound_index in range(len(lower_bounds)):
                    if upper_bounds[bound_index] <= lower_bounds[bound_index]:
                        lower_bounds[bound_index] = -np.inf
                        upper_bounds[bound_index] = np.inf
                    if oneDFitGuess[bound_index] < lower_bounds[bound_index]: oneDFitGuess[bound_index] = lower_bounds[bound_index]
                    if oneDFitGuess[bound_index] > upper_bounds[bound_index]: oneDFitGuess[bound_index] = upper_bounds[bound_index]
                try:
                    fit_res = scipy.optimize.curve_fit(fitting_funct, x_bin_centers_to_fit, binned_n_visits_to_fit, p0 = oneDFitGuess, bounds = (lower_bounds, upper_bounds) )[0]
                except RuntimeError:
                    print ('Failed to fit binning of single variable ' + var + '. Using null fit.')

                print ('fit_res = ' + str(fit_res))
                #If the edge bins is larger than the peak bin, we should just center the gaussian fit at that edge
                #if max([binned_n_visits[0], binned_n_visits[-1]]) > fit_res[0]: fit_res[1] = x_bin_centers[np.argmax([binned_n_visits[0], binned_n_visits[-1]])]
                self.gaussianPosteriorFits[spot_index_to_disp][var] = fit_res[:]

                #fig2 = plt.figure(2)
                #plt.plot(x_bin_centers_to_fit, binned_n_visits_to_fit)
                #plt.plot(x_bin_centers, fitting_funct(x_bin_centers, *oneDFitGuess), c = 'r')
                #plt.plot(x_bin_centers, fitting_funct(x_bin_centers, *fit_res), c = 'g')
                #plt.show()

                local_max_n_visits = max(binned_n_visits)
                max_n_visits = max(max_n_visits, local_max_n_visits)
                n_visits_disp_max = max_n_visits * 1.01
                if x_index < len(vars_to_disp):
                    axarr[-1][x_index].plot(x_bin_centers, binned_n_visits)
                    axarr[-1][x_index].plot(x_bin_centers, fitting_funct(x_bin_centers, *fit_res), c = 'r')
                    axarr[-1][x_index].text(x_bin_centers[1], local_max_n_visits * 1.2, r'$(\mu, \sigma)$ = ' + '\n' + str([c.round_to_n(fit_res[1], 2), c.round_to_n(fit_res[2], 2)]))
                if y_index >=  0:
                    axarr[y_index][0].plot(binned_n_visits, x_bin_centers)
                    axarr[y_index][0].plot(fitting_funct(x_bin_centers, *fit_res), x_bin_centers, c = 'r')
                    axarr[y_index][0].text(local_max_n_visits * 1.2, x_bin_centers[1], r'$(\mu, \sigma)$ = ' + '\n' + str([c.round_to_n(fit_res[1], 2), c.round_to_n(fit_res[2], 2)]))
                best_fit_params[var] = fit_res[1]
                if best_fit_params[var] < min(x_bin_centers): best_fit_params[var] = min(x_bin_centers)
                if best_fit_params[var] > max(x_bin_centers): best_fit_params[var] = max(x_bin_centers)
            print ("[best_fit_params[var_to_disp] for var_to_disp in vars_to_disp] = " + str([best_fit_params[var_to_disp] for var_to_disp in vars_to_disp]))
            best_fit_chi_sqr = self.calcChiSqrOfParams(angle_seps, best_fit_params['dip'], best_fit_params['mono'], best_fit_params['mu_dip'], best_fit_params['mu_mono'])
            null_chi_sqr = self.null_chi_sqr #self.calcChiSqrOfParams(angle_seps, 0.0, 0.0, best_fit_params['power'], best_fit_params['scaling'])
            f.suptitle(r'best $\chi^2$ = ' + str(c.round_to_n(best_fit_chi_sqr, 5)) + r'; null $\chi^2 = $' + str(c.round_to_n(null_chi_sqr,5)))

            if n_var_per_plot is 1:
                f.text(0.04, 0.5, 'Number of Chain Steps', va='center', rotation='vertical')

            else :
                axarr[-1][0].set_ylim([0.0, n_visits_disp_max])
                axarr[-1][0].set_xlim([0.0, n_visits_disp_max])
                bin_fractions = [1.0/6.0, 2.0/6.0, 3.0/6.0, 4.0/6.0, 5.0/6.0]
                bin_fractions = [1.0/6.0, 3.0/6.0, 5.0/6.0 ]
                axarr[-1][0].set_xticks([int(n_visits_disp_max * i) for i in [1.0/6.0, 3.0/6.0, 5.0/6.0 ]])
                axarr[-1][0].set_yticks([int(n_visits_disp_max * i) for i in [1.0/6.0, 3.0/6.0, 5.0/6.0 ]])
                binned_xticks = [int(n_visits_disp_max * i ) for i in [1.0/6.0, 3.0/6.0, 5.0/6.0 ]]
                binned_yticks = [int(n_visits_disp_max * i ) for i in [1.0/6.0, 3.0/6.0, 5.0/6.0 ]]
                binned_xticklabels = [int(n_visits_disp_max * i / 1000.0) for i in [1.0/6.0, 3.0/6.0, 5.0/6.0 ]]
                binned_yticklabels = [int(n_visits_disp_max * i / 1000.0) for i in [1.0/6.0, 3.0/6.0, 5.0/6.0 ]]
                axarr[-1][0].set_xticklabels(binned_xticklabels, fontsize = 10.0 * fig_size_unit * n_y_contour_plots / (n_y_ticks), rotation = 0)
                axarr[-1][0].set_yticklabels(binned_yticklabels, fontsize = 10.0 * fig_size_unit * n_x_contour_plots / (n_x_ticks), rotation = 0)
                bin_label = r'$N_{\mathrm{bin}}$' + ' ' + r'$(10^3)$'
                axarr[-1][0].set_xlabel(bin_label, fontsize = 2.0 * (fig_size_unit) * n_y_contour_plots)
                axarr[-1][0].set_ylabel(bin_label, fontsize = 2.0 * (fig_size_unit) * n_y_contour_plots)
                axarr[axarr_y_index][0].yaxis.set_label_coords(-0.4, 0.5)


                for i in range(1, total_n_x_plots):
                    axarr[-1][i].set_ylim([0.0, n_visits_disp_max])
                    axarr[-1][i].set_yticks(binned_yticks)
                    axarr[-1][i].set_yticklabels([])
                    axarr[i-1][0].set_xlim([0.0, n_visits_disp_max])
                    axarr[i-1][0].set_xticks(binned_xticks)
                    axarr[i-1][0].set_xticklabels([])

        print ('save_fig = '  + str(save_fig) )
        if save_fig:
            if file_name == '': file_name = 'MCMC_results_' + ('randData' if self.randomize else 'trueData') + '_' + 'dipRA' + str(c.round_to_n(target_sky_coord[0], 3)) + 'Dec' + str(c.round_to_n(target_sky_coord[1], 3)) + '_vars_' + '_'.join([str(elem) for elem in vars_to_disp]) + '_' + str(n_var_per_plot) + '_per_plot' + '.png'
            print ('saving figure to ' + plot_dir + file_name)
            plt.savefig(plot_dir + file_name)

        if show:
            plt.show()
        else:
            plt.close('all')

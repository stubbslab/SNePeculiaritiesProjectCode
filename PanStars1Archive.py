#Define an archive specific to the Pan Stars-1 medium field Sn surveys.
# This can just help me to expedite analysis that want to focus just on those fields.
# There are 10 fields, which are (arbtrarily) numbers 0-9.  They can be specifically requested.

import AstronomicalParameterArchive as apa
from SNDataArchive import SNDataArchive
from RawSNDataStorer import RawSNDataStorer
import cantrips as c
import matplotlib.pyplot as plt
import CosmologicalParameterArchive as cpa
import scipy.integrate as integrate
import numpy as np
from PyAstronomy import pyasl
import time
from loadSN import loadSN
import scipy
import healpy as hp
from astropy import units as u


#Here's the basic use:
# $ python
# import PanStars1Archive as psa
# import cantrips as c
# import numpy as np
# import matplotlib.pyplot as plt
# sdssdir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/SNIsotropyProject/SDSSGalaxies/'
# fulldata_fastRead = c.readInColumnsToList('SDSS_fullCoverage_SDSSGals_pzAll.csv', sdssdir, delimiter = ',', n_ignore = 2, all_np_readable = 1)
# gal_plotter = psa.PanStars1Archive(full_sdss_gal_data_file = 'SDSS_PSSuperField3_SDSSGals_Allpz.csv', preloaded_sdss_gals = fulldata_fastRead)
# angular_seps = np.linspace(0.5, 10.0, 21)
# densities_by_ang_radius = []
# f, axarr = plt.subplots(7,3, figsize = (14, 10), sharex = True)
# for ang_sep in angular_seps:
#     print ('Working on ang_sep = ' + str(ang_sep))
#     new_rel_gal_densities = gal_plotter.getGalRelativeDensitiesInAngularRadiusAroundCoordInRedshiftBins([36.0, -4.0], ang_sep, max_bin_z = 0.4, n_z_bins = 41, )
#     densities_by_ang_radius = densities_by_ang_radius + [new_rel_gal_densities]
#
#for i in range(len(angular_seps)):
#    axarr[i//3, i % 3].scatter(gal_plotter.z_bin_edges[0:-1], densities_by_ang_radius[i])
#
#plt.show()

class PanStars1Archive:

    def loadInterp(self, interp_file):
        loaded_interp = np.load(self.SDSSDir + interp_file, allow_pickle = True).item()

        return loaded_interp

    #def loadRedshiftToProperDistInterp(self, interp_file = 'PlanckRedshifToComovingInterpForSDSS1.npy'):
    #    self.comoving_interp = np.load(self.SDSSDir + interp_file, allow_pickle = True).item()
    #
    #    return 1

    #We need to fill
    def sampleGalDensity(self, min_z, max_z, grid_comoving_sep, search_comov_rad_Mpc, sky_sample_frac = 1.0, verbose = 1):
        min_r, max_r = [self.z_to_comoving_interp(min_z), self.z_to_comoving_interp(max_z)]
        max_r = max_r + grid_comoving_sep
        print ('[min_r, max_r] = ' + str([min_r, max_r]))
        deg_to_rad = self.astro_arch.getDegToRad( )
        #Now set up a grid of points
        comoving_points_on_grid = (-1 * np.arange(min_r, max_r, grid_comoving_sep)).tolist() + np.arange(min_r, max_r, grid_comoving_sep).tolist()
        print ('comoving_points_on_grid')
        comoving_interp_x_mesh, comoving_interp_y_mesh, comoving_interp_z_mesh = np.meshgrid(comoving_points_on_grid, comoving_points_on_grid, comoving_points_on_grid)
        comoving_interp_rs = (comoving_interp_x_mesh ** 2.0 + comoving_interp_y_mesh ** 2.0 + comoving_interp_z_mesh ** 2.0) ** 0.5
        interp_phis = c.goodArctan(comoving_interp_x_mesh, comoving_interp_y_mesh) / deg_to_rad
        interp_thetas = np.arccos(comoving_interp_z_mesh / comoving_interp_rs) / deg_to_rad

        interp_redshifts = c.flattenListOfLists(c.flattenListOfLists(self.comoving_to_z_interp(comoving_interp_rs)))
        interp_RAs = c.flattenListOfLists(c.flattenListOfLists( interp_phis))
        interp_Decs = c.flattenListOfLists(c.flattenListOfLists( -interp_thetas + 90.0))
        interp_thetas = c.flattenListOfLists(c.flattenListOfLists(interp_thetas))

        #comoving_interp_rs = c.flattenListOfLists(c.flattenListOfLists(comoving_interp_rs))

        #interp_points = [[comoving_interp_rs[i], interp_RAs[i], interp_thetas[i]] for i in range(len(interp_redshifts)) if interp_redshifts[i] < np.inf]
        interp_points = [[interp_redshifts[i], interp_RAs[i], interp_Decs[i]] for i in range(len(interp_redshifts)) ]
        print ('len(interp_points) = ' + str(len(interp_points)))
        interp_points = [[interp_redshifts[i], interp_RAs[i], interp_Decs[i]] for i in range(len(interp_redshifts)) if (interp_redshifts[i] < max_z )]
        print ('len(interp_points) = ' + str(len(interp_points)))

        interp_healpix_indeces = self.hp.lonlat_to_healpix( [interp_point[1] for interp_point in interp_points] * u.deg, [interp_point[2] for interp_point in interp_points] * u.deg )
        points_in_footprint = np.in1d(interp_healpix_indeces, self.covered_pixels)
        interp_points = [interp_points[i] for i in range(len(interp_points)) if points_in_footprint[i] ]
        self.interp_points = interp_points
        print ('len(interp_points) = ' + str(len(interp_points)))

        gal_counts_at_sample_coords = [0 for interp_point in interp_points]
        for i in range(len(interp_points)):
            interp_point = interp_points[i]
            print ('Working on point number ' + str(i) + ' of ' + str(len(interp_points))  +  '... ')
            gal_counts_at_sample_coords[i] = len(self.countNGalsAroundCoord(*interp_point, search_comov_rad_Mpc), verbose = verbose)
        self.gal_counts_at_sample_coords = gal_counts_at_sample_coords

        sky_n_gals_interp = scipy.interpolate.LinearNDInterpolator(np.array(self.interp_points), gal_counts_at_sample_coords)
        #Now we need to normalize
        healpix_pixel_lonlats = [[self.hp.healpix_to_lonlat(i)[0].value, self.hp.healpix_to_lonlat(i)[1].value] for i in range(self.hp.npix)]
        interp_r_shells = np.arange(min_r, max_r, grid_comoving_sep)
        interp_z_shells = [float(self.comoving_to_z_interp(r)) for r in interp_r_shells]
        n_gals_on_shells = [[float(sky_n_gals_interp(z, pix[0] / deg_to_rad, pix[1] / deg_to_rad)) for pix in healpix_pixel_lonlats] for z in interp_z_shells]
        n_gals_on_shells = [[n_gal for n_gal in n_gals_on_shell if not(np.isnan(n_gal))] for n_gals_on_shell in n_gals_on_shells]
        shell_normalizations = [1.0 / (np.sum(n_gals_on_shell) / len(n_gals_on_shell)) if len(n_gals_on_shell) > 0 else 0.0 for n_gals_on_shell in n_gals_on_shells ]
        print ('shell_normalizations = ' + str(shell_normalizations))
        shell_normalizations_interp = scipy.interpolate.interp1d(interp_z_shells, shell_normalizations)
        self.norm_gal_counts_at_sample_coords = [self.gal_counts_at_sample_coords[i] * shell_normalizations_interp(self.interp_points[i][0]) for i in range(len(self.gal_counts_at_sample_coords))]
        self.norm_sky_n_gals_interp = scipy.interpolate.LinearNDInterpolator(np.array(self.interp_points), self.norm_gal_counts_at_sample_coords, fill_value = 0.0 )

        return gal_counts_at_sample_coords


    def autoFindRelativeGalDensityAroundCoordToSN(self,  targ_z, targ_RA, targ_Dec, target_SN, inner_annulus_ratio = 1, outer_annulus_ratio = 2, start_mpc_step = 2.0, mpc_step = 2.0):
         S_to_N = 0.0
         current_inner_sphere_rad = start_mpc_step = mpc_step
         while S_to_N < target_SN:
             current_inner_sphere_rad = current_inner_sphere_rad + mpc_step
             current_gal_n_dens, current_gal_n_dens_err = self.measureRelativeGalDensityAroundCoord(targ_z, targ_RA, targ_Dec, current_inner_sphere_rad, current_inner_sphere_rad * inner_annulus_ratio, current_inner_sphere_rad * outer_annulus_ratio)
             S_to_N = current_gal_n_dens / current_gal_n_dens_err
             print ('For inner ratio of ' + str(current_inner_sphere_rad) + ', S/N = ' + str(S_to_N) )

         return current_gal_n_dens, current_gal_n_dens_err, current_inner_sphere_rad


    def measureRelativeGalDensityAroundCoord(self,  targ_z, targ_RA, targ_Dec, inner_sphere_comoving_rad, annulus_inner_comoving_rad, annulus_outer_comoving_rad, verbose = 1):

        #local_gal_indeces = self.countNGalsInAnnulusAroundCoord(targ_z, targ_RA, targ_Dec, 0.0, inner_sphere_comoving_rad)
        #background_gal_indeces = self.countNGalsInAnnulusAroundCoord(targ_z, targ_RA, targ_Dec, annulus_inner_comoving_rad, annulus_outer_comoving_rad)
        #inner_ang_size = np.arcsin(inner_sphere_comoving_rad / np.sqrt(inner_sphere_comoving_rad ** 2.0 + targ_r ** 2.0))
        #outer_ang_size = np.arcsin(inner_sphere_comoving_rad / np.sqrt(inner_sphere_comoving_rad ** 2.0 + targ_r ** 2.0))

        local_gal_indeces = self.countNGalsInAnnulusAroundCoord(targ_z, targ_RA, targ_Dec, 0.0, inner_sphere_comoving_rad, verbose = verbose)
        background_gal_indeces = self.countNGalsInAnnulusAroundCoord(targ_z, targ_RA, targ_Dec, annulus_inner_comoving_rad, annulus_outer_comoving_rad, verbose = verbose)
        local_n_gal, background_n_gal = [len(local_gal_indeces), len(background_gal_indeces)]
        local_n_gal_err, background_n_gal_err = [np.sqrt(local_n_gal), np.sqrt(background_n_gal)]

        inner_vol = self.getCoveredVolumeInAnnulusAroundCoord(targ_z, targ_RA, targ_Dec, 0.0, inner_sphere_comoving_rad)
        outer_vol = self.getCoveredVolumeInAnnulusAroundCoord(targ_z, targ_RA, targ_Dec, annulus_inner_comoving_rad, annulus_outer_comoving_rad)
        vols_from_geometry = [4.0 * np.pi / 3.0 * inner_sphere_comoving_rad ** 3.0, 4.0 * np.pi / 3.0 * (annulus_outer_comoving_rad ** 3.0 - annulus_inner_comoving_rad ** 3.0)]
        #print ('vols_from_geometry = ' + str(vols_from_geometry))
        #print ('vols_from_pixelation = ' + str([inner_vol, outer_vol]))
        comoving_vol_ratio = inner_vol / outer_vol
        #print ('[local_n_gal, background_n_gal, inner_vol, outer_vol] = ' + str([local_n_gal, background_n_gal, inner_vol, outer_vol]))

        if len(background_gal_indeces) > 0 and (inner_vol > 0 and outer_vol > 0):
            gal_n_density_ratio = local_n_gal / background_n_gal / comoving_vol_ratio
            if local_n_gal > 0:
                gal_n_density_ratio_err = np.sqrt((local_n_gal_err / background_n_gal) ** 2.0 +  (local_n_gal * background_n_gal_err / (background_n_gal ** 2.0)) ** 2.0) / comoving_vol_ratio
            else:
                gal_n_density_ratio_err = 1.0
        else:
            gal_n_density_ratio, gal_n_density_ratio_err = [1,1]
        #print ('[len(local_gal_indeces), len(background_gal_indeces), comoving_vol_ratio, gal_n_density_ratio] = ' + str([len(local_gal_indeces), len(background_gal_indeces), comoving_vol_ratio, gal_n_density_ratio]))
        #print('local n gal over/under density = ' + str(gal_n_density_ratio) + r'+/-' + str(gal_n_density_ratio_err))
        return [gal_n_density_ratio, gal_n_density_ratio_err]

    def computeIntersectionsOfLineOfSightsWithSphere(self, sphere_center_r, sphere_rad, ang_sep):
        term1 = sphere_center_r * np.cos(ang_sep)
        term2 = (sphere_center_r *  np.cos(ang_sep)) ** 2.0 - (sphere_center_r ** 2.0 - sphere_rad ** 2.0)
        if term2 >= 0.0:
            intersecting_rs_by_pix = [term1 - np.sqrt(term2), term1 + np.sqrt(term2)]
        else:
            intersecting_rs_by_pix = [sphere_center_r, sphere_center_r]
        return intersecting_rs_by_pix

    def getCoveredVolumeInAnnulusAroundCoord(self, targ_z, targ_RA, targ_Dec, annulus_inner_comoving_rad, annulus_outer_comoving_rad):
        deg_to_rad = self.astro_arch.getDegToRad()
        targ_r = self.z_to_comoving_interp(targ_z)
        outer_ang_size_rad = np.arcsin(annulus_outer_comoving_rad / np.sqrt(annulus_outer_comoving_rad ** 2.0 + targ_r ** 2.0))
        targ_RA, targ_Dec = c.moduloSkyCoords(targ_RA, targ_Dec)

        outer_ipix = hp.query_disc(nside=self.n_healpix_sides, vec=hp.ang2vec((90.0 - targ_Dec) * deg_to_rad, (targ_RA) * deg_to_rad), radius=outer_ang_size_rad, inclusive = True)
        #print ('[outer_ang_size_rad, self.n_healpix_sides, hp.ang2vec((90.0 - targ_Dec) * deg_to_rad, (targ_RA) * deg_to_rad), outer_ipix] = ' + str([outer_ang_size_rad, self.n_healpix_sides, hp.ang2vec((90.0 - targ_Dec) * deg_to_rad, (targ_RA) * deg_to_rad), outer_ipix]))
        pixels_included = np.in1d(outer_ipix, self.covered_pixels)
        #print ('pixels_included = ' + str(pixels_included))
        covered_pix = [outer_ipix[i] for i in range(len(outer_ipix)) if pixels_included[i]]
        covered_pix_centers = [hp.pix2ang(nside=self.n_healpix_sides, ipix = pix) for pix in covered_pix]
        covered_pix_centers = [[center[1] / deg_to_rad, -(center[0] - np.pi / 2.0) / deg_to_rad] for center in covered_pix_centers]
        covered_pix_ang_seps = [c.measureAngularSeparationOnSky([targ_RA, targ_Dec], covered_pix_center, return_radian = 1) for covered_pix_center in covered_pix_centers]
        inner_intersecting_rs_by_pix = [self.computeIntersectionsOfLineOfSightsWithSphere(targ_r, annulus_inner_comoving_rad, covered_pix_ang_sep) for covered_pix_ang_sep in covered_pix_ang_seps]
        outer_intersecting_rs_by_pix = [self.computeIntersectionsOfLineOfSightsWithSphere(targ_r, annulus_outer_comoving_rad, covered_pix_ang_sep) for covered_pix_ang_sep in covered_pix_ang_seps]
        outer_volume_by_pix = [self.healpix_area / (4.0 * np.pi) * np.pi * 4.0 / 3.0 * (outer_rs[1] ** 3.0 - outer_rs[0] ** 3.0) for outer_rs in outer_intersecting_rs_by_pix]
        inner_volume_by_pix = [self.healpix_area / (4.0 * np.pi) * np.pi * 4.0 / 3.0 * (inner_rs[1] ** 3.0 - inner_rs[0] ** 3.0) for inner_rs in inner_intersecting_rs_by_pix]
        total_volume = np.sum(outer_volume_by_pix) - np.sum(inner_volume_by_pix)
        #print ('[inner_volume_by_pix, outer_volume_by_pix, total_volume] = ' + str([inner_volume_by_pix, outer_volume_by_pix, total_volume]))

        return total_volume

    #I need to check which galaxies are in the pixelated volumes, which approximate but are not strictly the same as a true spherical radius
    def countNGalsInAnnulusAroundCoord(self, targ_z, targ_RA, targ_Dec, annulus_inner_comoving_rad, annulus_outer_comoving_rad, verbose = 1):
        start1 = time.time()
        deg_to_rad = self.astro_arch.getDegToRad()
        targ_r = float(self.z_to_comoving_interp(targ_z))
        #print ('targ_r = ' + str(targ_r))
        #First due a course cut so that we do not need to compute properties for every galaxy
        min_z = self.comoving_to_z_interp(max(targ_r - annulus_outer_comoving_rad, 0.0))
        max_z = self.comoving_to_z_interp(targ_r + annulus_outer_comoving_rad)
        all_zs = np.array(self.sdss_all_gals[self.redshift_index])
        min_z_bound_index = 0
        min_r, max_r = [self.z_to_comoving_interp(min_z), self.z_to_comoving_interp(max_z)]

        first_pass_z_cuts = []
        rough_z_bound_indeces = [[self.sdss_z_indeces[z] for z in list(self.sdss_z_indeces.keys()) if z < min_z],
                                 [self.sdss_z_indeces[z] for z in list(self.sdss_z_indeces.keys()) if z > max_z]]
        #print ('rough_z_bound_indeces = ' + str(rough_z_bound_indeces))
        rough_z_bound_indeces = [rough_z_bound_indeces[0][-1] if len(rough_z_bound_indeces[0]) > 0 else 0,
                                 rough_z_bound_indeces[1][0] if len(rough_z_bound_indeces[1]) > 0 else self.n_sdss_gals - 1]
        #print ('rough_z_bound_indeces = ' + str(rough_z_bound_indeces))

        #z_bound_indeces = [np.where(all_zs > min_z, 1, 0).tolist().index(1), np.where(all_zs > max_z, 1, 0).tolist().index(1) - 1]
        #print ('old z_bound_indeces = ' + str(z_bound_indeces))
        zs_above_min = np.where(all_zs[rough_z_bound_indeces[0]:rough_z_bound_indeces[1]] > min_z, 1, 0).tolist()
        zs_below_max = np.where(all_zs[rough_z_bound_indeces[0]:rough_z_bound_indeces[1]] > max_z, 1, 0).tolist()
        if 1 in zs_above_min and 1 in zs_below_max:
            z_bound_indeces = [zs_above_min.index(1) + rough_z_bound_indeces[0], zs_below_max.index(1) + rough_z_bound_indeces[0] - 1]
        else:
            z_bound_indeces = [0,0]
        #print ('new z_bound_indeces = ' + str(z_bound_indeces))
        if z_bound_indeces[1] < 0: z_bound_indeces[1] = 0

        end1 = time.time()
        #print ('end1 - start1 = ' + str(end1 - start1) + 's')

        start2 = time.time()
        min_ang_sep = np.arcsin(annulus_outer_comoving_rad / np.sqrt(targ_r ** 2.0 + annulus_outer_comoving_rad ** 2.0))
        possible_gal_indeces = list(range(*z_bound_indeces))
        #print ('possible_gal_indeces = ' + str(possible_gal_indeces))
        dec_ang_seps = [abs(dec - targ_Dec) * deg_to_rad for dec in np.array(self.sdss_all_gals[self.dec_index])[z_bound_indeces[0]:z_bound_indeces[1]] ]
        possible_gal_indeces = [possible_gal_indeces[i] for i in range(len(dec_ang_seps)) if dec_ang_seps[i] < min_ang_sep]
        #print ('possible_gal_indeces = ' + str(possible_gal_indeces))
        possible_zs = [self.sdss_all_gals[self.redshift_index][index] for index in possible_gal_indeces]
        possible_ras = [self.sdss_all_gals[self.ra_index][index] for index in possible_gal_indeces]
        possible_decs = [self.sdss_all_gals[self.dec_index][index] for index in possible_gal_indeces]
        #print ('possible_ras = ' + str(possible_ras))
        possible_ang_seps = [c.measureAngularSeparationOnSky([targ_RA, targ_Dec], [possible_ras[i], possible_decs[i]], return_radian = 1) for i in range(len(possible_ras))]

        #print ('min_ang_sep = ' + str(min_ang_sep))
        #print ('possible_ang_seps = ' + str(possible_ang_seps))
        renarrowed_possible_gal_indeces = [i for i in range(len(possible_gal_indeces)) if min_ang_sep > possible_ang_seps[i]]
        possible_gal_indeces = [possible_gal_indeces[index] for index in renarrowed_possible_gal_indeces]
        possible_zs = [possible_zs[index]  for index in renarrowed_possible_gal_indeces]
        possible_ras = [possible_ras[index]  for index in renarrowed_possible_gal_indeces]
        possible_decs = [possible_decs[index]  for index in renarrowed_possible_gal_indeces]
        possible_ang_seps = [possible_ang_seps[index]  for index in renarrowed_possible_gal_indeces]

        #print ('len(renarrowed_possible_gal_indeces) = ' + str(len(renarrowed_possible_gal_indeces) ))
        possible_rs = [float(self.z_to_comoving_interp(z)) for z in possible_zs]
        possible_gal_pix = hp.ang2pix(nside=self.n_healpix_sides, theta = (90.0 - np.array(possible_decs)) * deg_to_rad, phi = np.array(possible_ras) * deg_to_rad)
        #print ('possible_gal_pix = ' + str(possible_gal_pix))
        #print ('[possible_gal_indeces, possible_zs, possible_rs, possible_ras, possible_decs, possible_ang_seps] = ' + str([possible_gal_indeces, possible_zs, possible_rs, possible_ras, possible_decs, possible_ang_seps]))


        inner_rs_by_pix = [self.computeIntersectionsOfLineOfSightsWithSphere(targ_r, annulus_inner_comoving_rad, ang_sep) for ang_sep in possible_ang_seps]
        outer_rs_by_pix = [self.computeIntersectionsOfLineOfSightsWithSphere(targ_r, annulus_outer_comoving_rad, ang_sep) for ang_sep in possible_ang_seps]

        #comoving_seps = [np.sqrt(targ_r ** 2.0 + possible_rs[i] ** 2.0 - 2.0 * targ_r * possible_rs[i] * np.cos(possible_ang_seps[i])) for i in range(len(possible_ang_seps))]
        #print ('comoving_seps = ' + str(comoving_seps))
        #print ('[[possible_rs[i], inner_rs_by_pix[i], outer_rs_by_pix[i]] for i in range(len(possible_rs))] = ' + str([[possible_rs[i], inner_rs_by_pix[i], outer_rs_by_pix[i]] for i in range(len(possible_rs))]))
        final_possible_gal_indeces = [possible_gal_indeces[i] for i in range(len(possible_rs)) if (possible_rs[i] > outer_rs_by_pix[i][0] and possible_rs[i] <= inner_rs_by_pix[i][0]) or (possible_rs[i] < outer_rs_by_pix[i][1] and possible_rs[i] >= inner_rs_by_pix[i][1])]
        #if verbose: print ('Found ' + str(len(final_possible_gal_indeces)) + ' gals around [redshift, RA, Dec]: ' + str([targ_z, targ_RA, targ_Dec]) + ' with inner/outer Mpc comoving annuli: ' + str([annulus_inner_comoving_rad, annulus_outer_comoving_rad]))
        end2 = time.time()
        #print ('end2 - start2 ' + str(end2 - start2) + 's')

        return final_possible_gal_indeces


    #From all SDSS gals w/ photo-zs, measure the number of galaxies in an angular radius within an angular diameter (given in degrees) fixed line of sight comoving bin
    def getGalRelativeDensitiesInAngularRadiusAroundCoordInRedshiftBins(self, target_coord, max_angular_sep, print_frequency = 10000, min_bin_z = 0.0, max_bin_z = 1.0, n_z_bins = 10, verbose = 1):
        sdss_full_angular_size = 1.0 #I need to actually calculating the SDSS photo-z footprint.  But to look at relative densities, this should be fine
        self.getGalsInAngularRadiusAroundCoordInRedshiftBins(target_coord, max_angular_sep, print_frequency = print_frequency, min_bin_z = min_bin_z, max_bin_z = max_bin_z, n_z_bins = n_z_bins)
        deg_to_rad = self.astro_arch.getDegToRad()
        los_angular_area = 2.0 * np.pi * (1.0 - np.cos(max_angular_sep * deg_to_rad))
        los_sdss_angular_area = los_angular_area #I need to figure out how to calculate the intersection of the SDSS survey with a given line of sight angular opening
        gal_rel_densities_in_redshift_bins = [-1.0 for i in range(n_z_bins)]
        for i in range(n_z_bins):
            n_sdss_gals_in_bin = self.gal_index_ranges_in_z_bin[i][1] - self.gal_index_ranges_in_z_bin[i][0] + 1
            n_sdss_gals_in_bin_in_los = len(self.gal_indeces_in_ang_rad_by_z_bin[i])
            gal_rel_densities_in_redshift_bin = (n_sdss_gals_in_bin_in_los  / los_sdss_angular_area) / (n_sdss_gals_in_bin / sdss_full_angular_size) # both volumes have an r^2 element, which cancel out
            gal_rel_densities_in_redshift_bins[i] = gal_rel_densities_in_redshift_bin

        return gal_rel_densities_in_redshift_bins


    #From all SDSS gals w/ photo-zs, measure the number of galaxies in a comoving radius within a fixed line of sight comoving bin
    def getGalsInAngularRadiusAroundCoordInRedshiftBins(self, target_coord, max_angular_sep, print_frequency = 10000, min_bin_z = 0.0, max_bin_z = 1.0, n_z_bins = 10, ):
        self.getGalsInRedshifts(min_bin_z = min_bin_z, max_bin_z = max_bin_z, n_z_bins = n_z_bins, show_gal_slices = 0)
        self.gal_indeces_in_ang_rad_by_z_bin = [ [] for bins in self.gal_index_ranges_in_z_bin ]
        if self.target_coord != target_coord:
            self.target_coord = target_coord

        for i in range(len(self.gal_index_ranges_in_z_bin)):
            #angular_seps = self.angular_seps_by_z_bin[i]
            print ('Calculating angular separation for galaxie in z-slice number ' + str(i) + ' of ' + str(len(self.gal_index_ranges_in_z_bin)))
            z_bins = self.gal_index_ranges_in_z_bin[i]
            angular_seps = [c.measureAngularSeparationOnSky(self.target_coord, [self.sdss_all_gals[self.ra_index][j], self.sdss_all_gals[self.dec_index][j]], return_radian = 0) for j in range(z_bins[0], z_bins[1])]
            #print ('angular_seps = ' + str(angular_seps))
            in_radius = np.where(np.array(angular_seps) <= max_angular_sep, 1, 0)
            indeces_in_radius = [index for index in range(z_bins[0], z_bins[1]) if in_radius[index - z_bins[0]]]
            self.gal_indeces_in_ang_rad_by_z_bin[i] = indeces_in_radius

        return 1

    def getGalsInRedshifts(self, min_bin_z = 0.0, max_bin_z = 1.0, n_z_bins = 10, show_gal_slices = 0, display_pause_time = 1.0, figsize = (9, 4)):
        self.gal_index_ranges_in_z_bin = [[0,0] for i in range(n_z_bins)]
        gal_zs = self.sdss_all_gals[self.redshift_index]
        max_z = gal_zs[-1]
        self.z_bin_edges = np.linspace(min_bin_z, max_bin_z, n_z_bins + 1)
        for i in range(n_z_bins):
            upper_z_edge = self.z_bin_edges[i+1]
            for j in range(self.gal_index_ranges_in_z_bin[i][1], len(gal_zs)):
                if gal_zs[j] > upper_z_edge:
                    print ('[j, gal_zs[j], upper_z_edge] = ' + str([j, gal_zs[j], upper_z_edge]))
                    self.gal_index_ranges_in_z_bin[i][1] = j-1
                    if i < n_z_bins - 1:
                        self.gal_index_ranges_in_z_bin[i+1][0] = j
                    break
            if j == len(gal_zs) - 1:
                self.gal_index_ranges_in_z_bin[i][1] = j
            j = 0
        print ('self.gal_index_ranges_in_z_bin = ' + str(self.gal_index_ranges_in_z_bin))
        if show_gal_slices:
            for i in range(len(self.gal_index_ranges_in_z_bin)):
                f = plt.figure(figsize=figsize)
                gal_indeces = self.gal_index_ranges_in_z_bin[i]
                ras = self.sdss_all_gals[self.ra_index][gal_indeces[0]:gal_indeces[1]]
                decs = self.sdss_all_gals[self.dec_index][gal_indeces[0]:gal_indeces[1]]
                c.plotStarsOnSky(ras, decs, show_fig = 0, save_fig = 0, fig = f, fig_indeces = 111, marker = '.')
                plt.draw()
                plt.pause(display_pause_time)
                plt.close('all')


    def generateSDSSFullGalDensHistogram(self, show = 1, save = 0, n_bins = 100, xlabel = 'Shell central photo-z ', ylabel = r'$\log_{10}$ (density of gals in shell $\times$ Mpc$^{3}$)', plt_title = 'SDSS photo-zs in shell volumes (DR-12)',
                                         label_size = 14, title_size = 18, color = 'k',
                                         save_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/SNIsotropyProject/plots/', plt_name = 'SDSS_all_photo_zs_density_in_z.pdf'):
        all_zs = self.sdss_all_gals[self.redshift_index]
        #all_rs = self.z_to_comoving_interp(all_zs)

        min_r, max_r = [self.z_to_comoving_interp(all_zs[0]), self.z_to_comoving_interp(all_zs[-1])]
        r_bin_edges = np.linspace(min_r, max_r, n_bins + 1)
        r_bin_centers = [(r_bin_edges[i] + r_bin_edges[i-1]) / 2.0 for i in range(1, len(r_bin_edges))]
        r_bin_volumes = [(r_bin_edges[i] ** 3.0 - r_bin_edges[i-1] ** 3.0) * 4.0 * np.pi / 3.0 for i in range(1, len(r_bin_edges)) ]
        z_bin_edges = self.comoving_to_z_interp(r_bin_edges)
        z_bin_centers = self.comoving_to_z_interp(r_bin_centers)
        #print ('z_bin_edges = ' + str(z_bin_edges))
        #print ('all_zs[-1] = ' + str(all_zs[-1]))

        n_zs = [0 for num in range(len(z_bin_edges) - 1)]
        current_bin_edge_index = 0
        for i in range(len(all_zs)):
            z = all_zs[i]
            #if i % 100000 == 0 : print ('On i = ' + str(i) + ' of ' + str(len(all_zs)) + ' z = ' + str(z))
            if z > z_bin_edges[current_bin_edge_index]:
                #print ('[z, current_bin_edge_index, z_bin_edges[current_bin_edge_index]] = ' + str([z, current_bin_edge_index, z_bin_edges[current_bin_edge_index]] ))
                current_bin_edge_index = current_bin_edge_index + 1
            n_zs[current_bin_edge_index - 1]= n_zs[current_bin_edge_index - 1] + 1
        print ('n_zs = ' + str(n_zs))

        plt.scatter(z_bin_centers, np.log10(np.array(n_zs) / np.array(r_bin_volumes) ), c = color)
        plt.xlabel(xlabel, fontsize = label_size)
        plt.ylabel(ylabel, fontsize = label_size)
        plt.title(plt_title, fontsize = title_size)
        plt.tight_layout()
        if save:
            plt.savefig(save_dir + plt_name)
        if show:
            plt.show()
        plt.close('all')
        return 1


    def generateSDSSFullZGalHistogram(self, show = 1, save = 0, n_bins = 100, xlabel = 'Photo-z bin', ylabel = '# gals in bin', plt_title = 'SDSS photo-zs (DR-12)',
                                      label_size = 14, title_size = 18,
                                     save_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/SNIsotropyProject/plots/', plt_name = 'SDSS_all_photo_zs_hist.pdf'):
        all_zs = self.sdss_all_gals[self.redshift_index]

        hist_res = plt.hist(all_zs, bins= n_bins, fill = None)
        plt.xlabel(xlabel, fontsize = label_size)
        plt.ylabel(ylabel, fontsize = label_size)
        plt.title(plt_title, fontsize = title_size)
        plt.tight_layout()
        if save:
            plt.savefig(save_dir + plt_name, bbox_inches='tight')
        if show:
            plt.show()
        plt.close('all')
        return 1

    def generateSDSSNGalHistograms(self, fields = [0, 2, 3, 4, 5, 6, 7, 8, 9], hist_cropped_field = 1, hist_full_field = 1, plot_array_rows = 5, hist_colors = ['k','r'],
                               n_bins = 50, figsize = (12,10)):
        ra_index, dec_index, redshift_index = [self.ra_index, self.dec_index, self.redshift_index]
        sdss_loads = {field:c.readInColumnsToList(self.SDSSPhotozs[field], self.SDSSDir, n_ignore = 2, delimiter= ',') for field in fields}
        f, axarr = plt.subplots(plot_array_rows, (len(fields) // plot_array_rows) + 1, sharex = True, sharey = True, figsize = figsize, squeeze = False)
        plt.subplots_adjust(hspace=0.0)
        plt.suptitle('SDSS Photo-zs in PS1MD Fields (white = 7 deg. diam. circle around field; orange = just field)')
        for i in range(len(fields)):
            ax = axarr[plot_array_rows - i % plot_array_rows - 1, i // plot_array_rows]
            field = fields[i]
            ras = [float(ra) for ra in sdss_loads[field][ra_index]]
            decs = [float(dec) for dec in sdss_loads[field][dec_index]]
            zs = [float(z) for z in sdss_loads[field][redshift_index]]
            #rs = self.z_to_comoving_interp(zs)
            field_bounds = self.fields[field]
            field_zs = [zs[i] for i in range(len(zs)) if (ras[i] > field_bounds[0] and ras[i] < field_bounds[1] and decs[i] > field_bounds[2] and decs[i] < field_bounds[3]) ]
            #field_rs = self.z_to_comoving_interp(field_zs)
            full_binned_zs, sub_binned_zs = [[], []]
            if hist_full_field:
                #full_n_gals_in_bin, bin_r_edges = np.histogram(rs, bins = n_bins)
                #bin_r_edges = self.z_to_comoving_interp(bin_z_edges)
                #central_rs = [(bin_r_edges[i] + bin_r_edges[i-1]) / 2.0 for i in range(1, len(bin_r_edges))]
                #central_rs = self.z_to_comoving_interp(central_zs)
                full_hist_output = ax.hist(zs, edgecolor = hist_colors[0], bins = n_bins, fill = False)
                #angular_area_sq_deg = 0.1
                #bin_volumes = [angular_area_sq_deg * deg_to_rad ** 2.0 * central_rs[i] ** 2.0 * (bin_r_edges[i+1] - bin_r_edges[i]) for i in range(len(central_rs))]
                #full_n_gal_per_mpc = full_n_gals_in_bin / bin_volumes
                #ax.scatter(central_rs, full_n_gals_in_bin / bin_volumes, c = hist_colors[0])
                #ax.scatter(central_zs, full_n_gals_in_bin / 1.0, c = hist_colors[0])
            if hist_cropped_field:
                #angular_area_sq_deg = 0.01
                #sub_n_gals_in_bin, bin_r_edges = np.histogram(field_rs, bins = n_bins)
                #bin_r_edges = self.z_to_comoving_interp(bin_z_edges)
                #central_rs = [(bin_r_edges[i] + bin_r_edges[i-1]) / 2.0 for i in range(1, len(bin_r_edges))]
                #central_rs = self.z_to_comoving_interp(central_zs)
                #bin_volumes = [angular_area_sq_deg * deg_to_rad ** 2.0 * central_rs[i] ** 2.0 * (bin_r_edges[i+1] - bin_r_edges[i]) for i in range(len(central_rs))]
                sub_hist_output = ax.hist(field_zs, edgecolor = hist_colors[0], bins = n_bins, fill = hist_colors[1])
                #sub_n_gal_per_mpc = sub_n_gals_in_bin / bin_volumes
                #ax.scatter(central_rs, sub_n_gals_in_bin / bin_volumes, c = hist_colors[1])
                #ax.scatter(central_zs, sub_n_gals_in_bin / 1.0, c = hist_colors[1])
            #n_gal_dens_ratio = sub_n_gal_per_mpc / full_n_gal_per_mpc
            #ax.scatter(central_rs, n_gal_dens_ratio, c = hist_colors[0])

            max_hist_val = np.max(full_hist_output[0].tolist() + sub_hist_output[0].tolist())
            if i % plot_array_rows == 0: ax.set_xlabel('SDSS photometric redshift')
            if i // plot_array_rows == 0: ax.set_ylabel(r'# Gals in bin in $\Delta \Omega$')
            text_coords = [0.5, max_hist_val * 0.75]
            ax.text(*text_coords, 'PS1MD Field ' + str(field))

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.show()
        return 1

    def generateSDSSRelativeGalDensityHistograms(self, fields = [0, 2, 3, 4, 5, 6, 7, 8, 9], hist_cropped_field = 1, hist_full_field = 1, plot_array_rows = 5, hist_colors = ['k','r'],
                                                  n_bins = 50, figsize = (12,10)):
        ra_index, dec_index, redshift_index = [self.ra_index, self.dec_index, self.redshift_index]
        deg_to_rad = self.astro_arch.getDegToRad()
        sdss_loads = {field:c.readInColumnsToList(self.SDSSPhotozs[field], self.SDSSDir, n_ignore = 2, delimiter= ',') for field in fields}
        f, axarr = plt.subplots(plot_array_rows, (len(fields) // plot_array_rows) + 1, sharex = True, sharey = True, figsize = figsize, squeeze = False)
        plt.subplots_adjust(hspace=0.0)
        plt.suptitle('SDSS Photo-zs in PS1MD Fields (white = 7 deg. diam. circle around field; orange = just field)')
        for i in range(len(fields)):
            ax = axarr[plot_array_rows - i % plot_array_rows - 1, i // plot_array_rows]
            field = fields[i]
            ras = [float(ra) for ra in sdss_loads[field][ra_index]]
            decs = [float(dec) for dec in sdss_loads[field][dec_index]]
            zs = [float(z) for z in sdss_loads[field][redshift_index]]
            rs = self.z_to_comoving_interp(zs)
            field_bounds = self.fields[field]
            field_zs = [zs[i] for i in range(len(zs)) if (ras[i] > field_bounds[0] and ras[i] < field_bounds[1] and decs[i] > field_bounds[2] and decs[i] < field_bounds[3]) ]
            field_rs = self.z_to_comoving_interp(field_zs)
            full_binned_zs, sub_binned_zs = [[], []]
            if hist_full_field:
                full_n_gals_in_bin, bin_r_edges = np.histogram(rs, bins = n_bins)
                #bin_r_edges = self.z_to_comoving_interp(bin_z_edges)
                central_rs = [(bin_r_edges[i] + bin_r_edges[i-1]) / 2.0 for i in range(1, len(bin_r_edges))]
                #central_rs = self.z_to_comoving_interp(central_zs)
                #full_hist_output = ax.hist(zs, edgecolor = hist_colors[0], bins = n_bins, fill = False)
                angular_area_sq_deg = 0.1
                bin_volumes = [angular_area_sq_deg * deg_to_rad ** 2.0 * central_rs[i] ** 2.0 * (bin_r_edges[i+1] - bin_r_edges[i]) for i in range(len(central_rs))]
                full_n_gal_per_mpc = full_n_gals_in_bin / bin_volumes
                #ax.scatter(central_rs, full_n_gals_in_bin / bin_volumes, c = hist_colors[0])
                #ax.scatter(central_zs, full_n_gals_in_bin / 1.0, c = hist_colors[0])
            if hist_cropped_field:
                angular_area_sq_deg = 0.01
                sub_n_gals_in_bin, bin_r_edges = np.histogram(field_rs, bins = n_bins)
                #bin_r_edges = self.z_to_comoving_interp(bin_z_edges)
                central_rs = [(bin_r_edges[i] + bin_r_edges[i-1]) / 2.0 for i in range(1, len(bin_r_edges))]
                #central_rs = self.z_to_comoving_interp(central_zs)
                bin_volumes = [angular_area_sq_deg * deg_to_rad ** 2.0 * central_rs[i] ** 2.0 * (bin_r_edges[i+1] - bin_r_edges[i]) for i in range(len(central_rs))]
                #        edgecolor = hist_colors[0], bins = n_bins, fill = hist_colors[1])
                sub_n_gal_per_mpc = sub_n_gals_in_bin / bin_volumes
                #ax.scatter(central_rs, sub_n_gals_in_bin / bin_volumes, c = hist_colors[1])
                #ax.scatter(central_zs, sub_n_gals_in_bin / 1.0, c = hist_colors[1])
            n_gal_dens_ratio = sub_n_gal_per_mpc / full_n_gal_per_mpc
            ax.scatter(central_rs, n_gal_dens_ratio, c = hist_colors[0])

            max_hist_val = np.max(full_n_gals_in_bin[0].tolist() + sub_n_gals_in_bin[0].tolist())
            if i % plot_array_rows == 0: ax.set_xlabel('SDSS comoving-r (Mpc) from photometric redshift')
            if i // plot_array_rows == 0: ax.set_ylabel(r'# Gals Mpc$^{-3}$')
            text_coords = [0.5, max_hist_val * 0.75]
            ax.text(*text_coords, 'PS1MD Field ' + str(field))

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.show()
        return 1

    def measureSDSSGalDensity(self, target_fields,
                              selection_rad = 'full', smoothing_rad_Mpc = 50, ra_index = 1, dec_index = 2, redshift_index = 11, field = 0, title_strs = [], figsize = [15, 7.5]):


        target_positions = [[(self.fields[field][1] + self.fields[field][0]) / 2, (self.fields[field][3] + self.fields[field][2]) / 2] for field in target_fields]
        cosmo_arch = cpa.CosmologicalParameterArchive()
        OmM = cosmo_arch.getOmegaM()[0]
        OmL = cosmo_arch.getOmegaLambda()[0]
        OmR = cosmo_arch.getOmegaR()[0]
        H_scaled = lambda z: np.sqrt((1.0 + z) ** 3.0 * OmM + (1.0 + z) ** 4.0 * OmR + OmL)
        c_light = cosmo_arch.getc() # in km/s
        H0 = cosmo_arch.getH0()[0]

        print ('Loading full galaxy archive...')
        start = time.time()
        sdss_load = c.readInColumnsToList(self.SDSSPhotozs[field],self.SDSSDir, n_ignore = 2, delimiter= ',')
        n_all_gals = len(sdss_load[0])
        end = time.time()
        print ('Took ' + str(end - start) + 's')
        print ('Computing parameters for ' + str(n_all_gals) + ' loaded galaxies...')
        start = time.time()
        all_ras = [float(ra) for ra in sdss_load[ra_index]]
        all_decs = [float(dec) for dec in sdss_load[dec_index]]
        all_zs = [float(z) for z in sdss_load[redshift_index]]
        max_z = max(all_zs)
        del sdss_load
        all_a_s = [1.0 / (1.0 + z) for z in all_zs]
        #all_comoving_dists = [c_light / H0 * integrate.quad(H_scaled, 0, z)[0]  for z in all_zs] # in units of Mpc
        all_comoving_dists = self.z_to_comoving_interp(all_zs)
        all_proper_dists = [ all_a_s[i] * all_comoving_dists[i] for i in range(len(all_zs)) ]
        end = time.time()
        print ('Took ' + str(end - start) + 's')
        print ('Selecting subgroup of galaxies in viewing radius...')
        start = time.time()
        if selection_rad in ['full','Full','FULL','all','All','ALL']:
            ras = [all_ras]
            decs = [all_decs]
            zs = [all_zs]
        else:
            target_ras = [[] for target in target_positions]
            target_decs = [[] for target in target_positions]
            target_zs = [[] for target in target_positions]
            target_a_s = [[] for target in target_positions]
            target_comoving = [[] for target in target_positions]
            target_proper = [[] for target in target_positions]
            for i in range(len(target_positions)):
                target_pos = target_positions[i]
                print ('target_pos = ' + str(target_pos))
                valid_indeces = [i for i in range(len(all_zs)) if pyasl.getAngDist(target_pos[0], target_pos[1], all_ras[i], all_decs[i]) < selection_rad]
                target_ras[i], target_decs[i], target_zs[i] = [ [all_ras[i] for i in valid_indeces], [all_decs[i] for i in valid_indeces], [all_zs[i] for i in valid_indeces] ]
                target_a_s[i] = [1.0 / (1.0 + z) for z in target_zs[i]]
                target_comoving[i] = self.z_to_comoving_interp(target_zs[i])
                target_proper[i] = [ target_a_s[i][j] * target_comoving[i][j] for j in range(len(target_zs[i])) ]
        #a_s = [1.0 / (1.0 + z) for z in zs]
        target_n_selected_gals = [len(zs) for zs in target_zs]
        #print ('[target_n_selected_gals, n_all_gals] = ' + str([target_n_selected_gals, n_all_gals]))
        #comoving_dists = [c_light / H0 * integrate.quad(H_scaled, 0, z)[0]  for z in zs] # in units of Mpc
        #comoving_dists = self.comoving_interp(zs)
        #proper_dists = [ a_s[i] * comoving_dists[i] for i in range(len(zs)) ]

        n_gal_starts = np.arange(min(all_comoving_dists), max(all_comoving_dists), smoothing_rad_Mpc)
        n_gal_centers = n_gal_starts + smoothing_rad_Mpc / 2.0
        end = time.time()
        print ('Took ' + str(end - start) + 's')
        print ('n_gal_starts = ' + str(n_gal_starts))
        comoving_bins = [[[] for bin_start in n_gal_starts] for target_position in target_positions]
        print ('Binning subselection galaxies...')
        start = time.time()
        for k in range(len(target_positions)):
            comoving_dists = target_comoving[k]
            for j in range(len(n_gal_starts)-1):
                gals_in_this_bin = (np.array(comoving_dists) > n_gal_starts[j]) * (np.array(comoving_dists) < n_gal_starts[j+1])
                comoving_bins[k][j] = [i for i in range(len(gals_in_this_bin)) if gals_in_this_bin[i] ]

        end = time.time()
        print ('Took ' + str(end - start) + 's')
        print ('Binning all galaxies...')
        start = time.time()
        all_comoving_bins = [[] for bin_start in n_gal_starts]
        for j in range(len(n_gal_starts)-1):
            gals_in_this_bin = (np.array(all_comoving_dists) > n_gal_starts[j]) * (np.array(all_comoving_dists) < n_gal_starts[j+1])
            all_comoving_bins[j] = [i for i in range(len(gals_in_this_bin)) if gals_in_this_bin[i] ]
            print ('For j = ' + str(j) + ': len(all_comoving_bins[j]) = ' + str(len(all_comoving_bins[j])))

        end = time.time()
        print ('Took ' + str(end - start) + 's')
        print ('Plotting...')
        n_gals_in_shells = []
        #plt.plot(n_gal_starts, [( len(all_comoving_bins[i]) - len(comoving_bins[i]) ) if len(all_comoving_bins[i]) > 0 else 0.0 for i in range(len(comoving_bins))], c = 'r')
        density_fluctuations = [[(( len(comoving_bins[i][j]) / target_n_selected_gals[i] - len(all_comoving_bins[j]) / n_all_gals) / (len(all_comoving_bins[j]) / n_all_gals))
                                            if len(all_comoving_bins[j]) > 0 else 0.0 for j in range(len(comoving_bins[i]))] for i in range(len(target_fields))]
        comoving_to_density_interps = [ scipy.interpolate.interp1d(n_gal_centers, density_fluctuations[i], bounds_error = False, fill_value = 0) for i in range(len(target_fields))]
        comoving_to_density_interps = [ scipy.interpolate.interp1d(n_gal_centers, 2.0 * (np.random.random(len(n_gal_starts)) - 0.5 ), bounds_error = False, fill_value = 0) for i in range(len(target_fields))]
        comoving_to_density_interps = [ scipy.interpolate.interp1d(n_gal_centers, [0 if (start < 500 or start > 1000) else -1 for start in n_gal_starts], bounds_error = False, fill_value = 0) for i in range(len(target_fields))]

        #comoving_to_density_interps = [ interpolate.interp1d(comoving_bins[i], density_fluctuations[i], bounds_error = False, fill_value = 0) for i in range(len(target_fields))]
        print ('n_gal_starts = ' + str(n_gal_starts))


        k = 0
        j = 0
        i = 0
        #print ('[comoving_bins[i], density_fluctuations[i]] = ' + str([comoving_bins[i], density_fluctuations[i]]))
        #print ('(np.array(n_gal_starts) + smoothing_rad_Mpc / 2.0) = ' +str((np.array(n_gal_starts) + smoothing_rad_Mpc / 2.0)))
        #print ('(comoving_to_density_interps[i]) = ' + str((comoving_to_density_interps[i])))
        #print ('(comoving_to_density_interps[i])(300.0) = ' + str((comoving_to_density_interps[i])(300.0) ))
        #print ('(comoving_to_density_interps[i])(np.array(n_gal_starts) + smoothing_rad_Mpc / 2.0) = ' + str((comoving_to_density_interps[i])(np.array(n_gal_starts) + smoothing_rad_Mpc / 2.0)))
        #print ('(1.0 + self.comoving_to_z_interp(np.array(n_gal_starts) + smoothing_rad_Mpc / 2.0)) = ' + str((1.0 + self.comoving_to_z_interp(np.array(n_gal_starts) + smoothing_rad_Mpc / 2.0))))
        #print ('((n_gal_starts[k] + smoothing_rad_Mpc / 2.0) / (n_gal_starts[j] + smoothing_rad_Mpc / 2.0) ) = ' + str(((n_gal_starts[k] + smoothing_rad_Mpc / 2.0) / (n_gal_starts[j] + smoothing_rad_Mpc / 2.0) )))
        #print ('((n_gal_starts[j] + smoothing_rad_Mpc / 2.0) - (n_gal_starts[k] + smoothing_rad_Mpc / 2.0) ) = ' + str(((n_gal_starts[j] + smoothing_rad_Mpc / 2.0) - (n_gal_starts[k] + smoothing_rad_Mpc / 2.0) )))

        convergences = [ [ (3.0 * (H0 / c_light) ** 2.0 * OmM / 2.0 * smoothing_rad_Mpc) *
                                                np.sum([comoving_to_density_interps[i](n_gal_centers[k]) * (1.0 + self.comoving_to_z_interp(n_gal_centers[k] )) * ((n_gal_centers[k]) / (n_gal_centers[j]) ) * ((n_gal_centers[j]) - (n_gal_centers[k]) ) for k in range(0,j+1)])
                                                  for j in range(len(n_gal_centers)) ] for i in range(len(target_fields)) ]

        j = len(n_gal_starts) - 1
        i = 0
        f, axarr = plt.subplots(6, 2)
        for index in range(0, 6):
            j = len(n_gal_starts) - 20 + index
            axarr[index][0].plot(n_gal_centers[0:j+1], [comoving_to_density_interps[i](n_gal_centers[k]) * (1.0 + self.comoving_to_z_interp(n_gal_centers[k] )) * ((n_gal_centers[k]) / (n_gal_centers[j]) ) * ((n_gal_centers[j]) - (n_gal_centers[k]) ) for k in range(0,j+1)])
            axarr[index][1].plot(n_gal_centers[0:j+1], [(1.0 + self.comoving_to_z_interp(n_gal_centers[k] )) * ((n_gal_centers[k]) / (n_gal_centers[j]) ) * ((n_gal_centers[j]) - (n_gal_centers[k]) ) for k in range(0,j+1)])
            print ('np.sum([comoving_to_density_interps[i](n_gal_centers[k]) * (1.0 + self.comoving_to_z_interp(n_gal_centers[k] )) * ((n_gal_centers[k]) / (n_gal_centers[j]) ) * ((n_gal_centers[j]) - (n_gal_centers[k]) ) for k in range(0,j+1)]) = ' + str(np.sum([comoving_to_density_interps[i](n_gal_centers[k]) * (1.0 + self.comoving_to_z_interp(n_gal_centers[k] )) * ((n_gal_centers[k]) / (n_gal_centers[j]) ) * ((n_gal_centers[j]) - (n_gal_centers[k]) ) for k in range(0,j+1)])))
        plt.show()


        convergences_interps = [ scipy.interpolate.interp1d(np.array(n_gal_centers), convergences[i], bounds_error = False, fill_value = convergences[i][-1]) for i in range(len(target_fields)) ]


        all_sns = loadSN(1, ['all'], pull_extinctions = 0)
        surveys_to_display = [sn['survey'] for sn in all_sns]
        sn_by_field = {}
        for field_num in target_fields:
            field = self.fields[field_num]
            sn_by_field[field_num] = [sn for sn in all_sns if ( sn['RA']>field[0] and sn['RA'] < field[1] and sn['Dec'] > field[2] and sn['Dec'] < field[3] )]

        f, axarr = plt.subplots(len(target_positions), 4, figsize = figsize)
        for i in range(len(target_positions)):
            field = target_fields[i]
            #comoving_bins = target_comoving[i]
            n_selected_gals = target_n_selected_gals[i]
            print ('n_gal_starts = ' + str(n_gal_centers))
            axarr[i][0].plot(np.array(n_gal_centers), [len(one_bin) / n_all_gals for one_bin in all_comoving_bins], c = 'b')
            axarr[i][0].plot(np.array(n_gal_centers), [len(one_bin) / n_selected_gals for one_bin in comoving_bins[i]], c = 'r')
            #axarr[0].plot(n_gal_starts, [( len(all_comoving_bins[i]) - len(comoving_bins[i]) ) if len(all_comoving_bins[i]) > 0 else 0.0 for i in range(len(comoving_bins))], c = 'purple')
            axarr[i][1].plot(np.array(n_gal_centers), comoving_to_density_interps[i](np.array(n_gal_centers)), c = 'k')  # we use total number of galaxies as stand in for sky area subtended
            axarr[i][1].plot(np.array(n_gal_centers), comoving_to_density_interps[i](np.array(n_gal_centers)), c = 'green')
            axarr[i][1].set_ylim(-1.0, 1.0)
            axarr[i][0].set_ylabel('Frac. of galaxies' + ' for f' + str(field) + ')')
            axarr[i][1].set_ylabel(r'$\delta$' + ' for f' + str(field))
            axarr[i][2].set_ylabel(r'$\Delta \mu$' + ' for f' + str(field) + ')')
            if i == len(target_positions) - 1:
                axarr[i][0].set_xlabel('Comoving (Mpc)')
                axarr[i][1].set_xlabel('Comoving (Mpc)')
                axarr[i][2].set_xlabel('Comoving (Mpc)')
            if len(title_strs) > 0:
                axarr[0][i].set_title(title_strs[i])
            sn_in_field = sn_by_field[field]
            sn_zs = [sn['z'] for sn in sn_in_field]
            sn_comoving = self.z_to_comoving_interp(sn_zs)
            sn_muDiffs = [sn['muDiff'] for sn in sn_in_field]
            sn_muErrs = [sn['muErr'] for sn in sn_in_field]
            sn_colors = [sn['color'] for sn in sn_in_field]
            print ('sn_comoving = ' + str(sn_comoving))
            axarr[i][2].scatter(sn_comoving, sn_muDiffs, c = sn_colors)
            axarr[i][2].errorbar(sn_comoving, sn_muDiffs, yerr = sn_muErrs, fmt = 'none', ecolor = sn_colors)
            x_lims = [0.0, max(sn_comoving) * 1.05]
            axarr[i][3].plot(n_gal_centers, (3.0 * (H0 / c_light) ** 2.0 * OmM / 2.0 * smoothing_rad_Mpc) * comoving_to_density_interps[i](np.array(n_gal_centers)) * (1.0 + self.comoving_to_z_interp(np.array(n_gal_centers)))  )
            axarr[i][2].plot(np.linspace(*x_lims, 1001), (convergences_interps[i])(np.linspace(*x_lims, 1001)), c = 'k')
            axarr[i][0].set_xlim(x_lims)
            axarr[i][1].set_xlim(x_lims)
            axarr[i][2].set_xlim(x_lims)
        end = time.time()
        print ('Took ' + str(end - start) + 's')
        plt.tight_layout()
        plt.show()

        return 1

    #def generateFullSDSSGalDensityHistogram(n_bins = 101):
    #    all_zs = self.sdss_all_gals[self.redshift_index]
    #    all_rs = self.z_to_comoving_interp(all_zs)


    def checkIfSkyCoordsInFootprint(self, RAs, Decs):
        deg_to_rad = self.astro_arch.getDegToRad()
        phis = np.array(RAs) * deg_to_rad
        thetas = (90.0 - np.array(Decs)) * deg_to_rad
        #print ('[thetas, phis] = ' + str([thetas, phis]))
        target_pixels = hp.ang2pix(self.n_healpix_sides, thetas, phis)
        #print ('[target_pixels, self.covered_pixels] = ' + str([target_pixels, self.covered_pixels]))
        points_in_footprint = np.in1d(target_pixels, self.covered_pixels)
        return points_in_footprint

    def setSDSSFootprintHEALPixMap(self, n_healpix_sides):
        deg_to_rad = self.astro_arch.getDegToRad()
        self.n_healpix_sides = n_healpix_sides
        self.n_healpix_pix = hp.nside2npix(self.n_healpix_sides)
        self.healpix_area = hp.pixelfunc.nside2pixarea(self.n_healpix_sides, degrees=False)
        #self.hp = HEALPix(nside=n_healpix_sides)
        gal_RAs = np.array(self.sdss_all_gals[self.ra_index])
        gal_Decs = np.array(self.sdss_all_gals[self.dec_index])
        print ('[max(gal_Decs * deg_to_rad), max(gal_RAs * deg_to_rad)] = ' + str([max(gal_Decs * deg_to_rad), max(gal_RAs * deg_to_rad)] ))
        self.covered_pixels = list(set(hp.ang2pix(self.n_healpix_sides, (90.0 - gal_Decs) * deg_to_rad  , gal_RAs * deg_to_rad)))
        #self.covered_pixels = self.hp.lonlat_to_healpix( )
        #print ('self.covered_pixels = ' + str(self.covered_pixels))

        return 1

    def showSDSSHealpixFootprint(self):
        pix_coloration = np.array([int(bool_val) for bool_val in np.in1d(np.array(list(range(self.n_healpix_pix))), self.covered_pixels)])
        print ('pix_coloration = ' + str(pix_coloration))
        print ('len(pix_coloration) = ' + str(len(pix_coloration)))

        #fig = plt.figure()
        m = np.array(list(range(self.n_healpix_pix)))
        print ('len(m) = ' + str(len(m)))
        res = hp.visufunc.mollview(pix_coloration, title="SDSS Footprint using HEALPix", cbar = False, fig = 111, hold = 1)
        hp.graticule()
        print ('res = ' + str(res))
        #c.plotStarsOnSky([], [], fig = fig)
        plt.show()

        return 1


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
        self.z_to_comoving_interp = scipy.interpolate.interp1d([-1000.0] + self.interp_zs, [0.0] + [self.r_of_z(z) for z in self.interp_zs], fill_value = (0.0, self.r_of_z(1000)), bounds_error = False)
        self.comoving_to_z_interp = scipy.interpolate.interp1d([0.0] + [self.r_of_z(z) for z in self.interp_zs], [-1000.0] + self.interp_zs, fill_value = (0.0, self.r_of_z(1000)), bounds_error = False)
        return 1

    def __init__(self, full_sdss_gal_data_file = 'SDSS_fullCoverage_SDSSGals_pzAll.csv', planck_cosmology = 0, H0 = 70.0, OmM = 0.3, OmL = 0.7, OmR = 0.0,
                 interp_z_params = [0.0, 100.0, 1001], ra_index = 1, dec_index = 2, redshift_index = 11, preloaded_sdss_gals = None, n_healpix_sides = 16, n_quick_index_zs = 100 ):
        self.ps_survey = 'PS1MD'
        self.astro_arch = apa.AstronomicalParameterArchive()
        self.sdss_full_gal_data_file = full_sdss_gal_data_file

        self.ra_index = ra_index
        self.dec_index = dec_index
        self.redshift_index = redshift_index
        self.target_coord = [np.nan, np.nan]

        #These are the on-sky regions that contain a given set of observations.
        #  They are not actual regions of observation; just regions in which all observations
        #  from a given field are contained.
        self.fields = {0:[34.5,37.5,-6.0,-2.5], #in: SDSS
                       1:[51.5,55.0,-29.5,-26.0],#in:
                       2:[128.0,133.0,43.0,45.5],#in: SDSS
                       3:[148.5,152.0,0.5,4.0],#in: SDSS
                       4:[159.5,164.0,56.5,59.5],#in: SDSS
                       5:[183.0,187.0,45.5,48.5],#in: SDSS
                       6:[211.5,216.5,51.0,54.5],#in: SDSS
                       7:[241.0,246.0,53.5,56.5],#in: SDSS
                       8:[332.0,336.0,-1.5,2.0],#in: SDSS
                       9:[350.0,354.0,-2.0,1.0] #in: SDSS
        }
        self.field_centers = {field: [(self.fields[field][1] + self.fields[field][0]) / 2.0, (self.fields[field][3] + self.fields[field][2]) / 2.0 ]  for field in self.fields.keys()}
        self.field_colors = {0:'blue', 1:'red', 2:'yellow', 3:'orange',
                             4:'cyan', 5:'limegreen', 6:'Peru', 7:'Salmon',
                             8:'Purple', 9:'Thistle'}

        self.SDSSDir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/SNIsotropyProject/SDSSGalaxies/'
        if preloaded_sdss_gals is None:
            self.sdss_all_gals = c.readInColumnsToList(self.sdss_full_gal_data_file, self.SDSSDir, n_ignore = 2, delimiter= ',', all_np_readable = 1)
        else:
            self.sdss_all_gals = preloaded_sdss_gals
        self.n_sdss_gals = len(self.sdss_all_gals[redshift_index])
        quick_index_indeces = [int(i * self.n_sdss_gals / n_quick_index_zs) for i in range(n_quick_index_zs)]
        quick_index_zs = [self.sdss_all_gals[self.redshift_index][index] for index in quick_index_indeces ]
        print ('quick_index_indeces = ' + str(quick_index_indeces))
        self.sdss_z_indeces = {z: np.where(self.sdss_all_gals[self.redshift_index] > z, 1, 0).tolist().index(1) for z in quick_index_zs }
        print ('self.sdss_z_indeces = ' + str(self.sdss_z_indeces))
        self.SDSSPhotozs = {key: 'SDSS_PSSuperField' + str(key) + '_SDSSGals_Allpz.csv' for key in [0, 2, 3, 4, 5, 6, 7, 8, 9]}
        #self.SDSSPhotozs[0] = 'SDSS_fullCoverage_SDSSGals_pzAll.csv'
        #self.loadRedshiftToProperDistInterp(interp_file = 'PlanckRedshifToComovingInterpForSDSS1.npy')
        if planck_cosmology:
            self.z_to_comoving_interp = self.loadInterp(interp_file = 'PlanckRedshifToComovingInterpForSDSS1.npy')
            self.comoving_to_z_interp = self.loadInterp(interp_file = 'PlanckComovingToRedshiftInterpForSDSS1.npy')
        else:
            self.H0 = H0
            self.OmM = OmM
            self.OmL = OmL
            self.OmR = OmR
            self.interp_z_params = interp_z_params
            self.initialize_r_of_z_interp(self.interp_z_params)

        self.setSDSSFootprintHEALPixMap(n_healpix_sides)

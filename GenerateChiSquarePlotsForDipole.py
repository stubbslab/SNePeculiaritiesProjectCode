import MeasureSNeDipole as msnd
import cantrips as c
import AstronomicalParameterArchive as apa
import numpy as np
import time

if __name__=="__main__":
    n_sky_points = 200 # 200 #200
    n_randomizations = 5
    n_init_param_sets = 8 #8
    n_MCMC_steps = 10000
    n_ignore = 1000
    astro_arch = apa.AstronomicalParameterArchive()
    deg_to_rad = astro_arch.getDegToRad()

    do_real = 1

    if do_real:
        start = time.time()
        d_fitter_real = msnd.SNDipoleFitter(n_sky_spots = n_sky_points, added_monopole_params = [0.00, 1.0, 0.15], rearrange_randomize = 0)

        true_n_sky_points = d_fitter_real.n_sky_spots
        print ('true_n_sky_points = ' + str(true_n_sky_points))

        for init_indeces in range(n_init_param_sets):
            print ('Working on init_indeces ' + str(init_indeces))
            d_fitter_real.measureSNDipoleViaMCMC(n_MCMC_steps, init_index = init_indeces)

        for i in range(true_n_sky_points):
            spot_on_sky = d_fitter_real.spots_on_sky[i]
            target_sky_coord = [c.goodArctan(spot_on_sky[0], spot_on_sky[1]) / deg_to_rad, (np.pi/2-np.arccos(spot_on_sky[2])) / deg_to_rad]
            file_name = 'MCMC_results_realData' + '_sepMonopoleAndDipoleFuncts' + '_dipRA' + str(c.round_to_n(target_sky_coord[0], 3)) + 'Dec' + str(c.round_to_n(target_sky_coord[1], 3)) + '.png'
            d_fitter_real.showMCMCResults(i, n_levels = 20, n_ignore = n_ignore, save_fig = 1, show = 0, file_name = file_name, vars_to_disp = ['dip', 'mono', 'mu_dip', 'mu_mono'] )

        d_fitter_real.showChiSqrOnSky(save_fig = 1, show_fig = 0, file_name = 'ChiSquareOfDipoleOnSky_realData.pdf',funct_params = ['dip', 'mono', 'mu_dip', 'mu_mono'])

        end = time.time()
        print ('<<<< TOOK ' + str(end - start) + 's FOR FULL REAL ITERATION >>>> ')

    for rand_index in range(n_randomizations):
        start = time.time()
        print ('!!!Working on rand_index ' + str(rand_index) + ' of ' + str(n_randomizations) + '!!!!')
        d_fitter_fake = msnd.SNDipoleFitter(n_sky_spots = n_sky_points, added_monopole_params = [0.00, 1.0, 0.15], rearrange_randomize = 1)

        true_n_sky_points = d_fitter_fake.n_sky_spots

        for init_indeces in range(8):
            d_fitter_fake.measureSNDipoleViaMCMC(n_MCMC_steps, init_index = init_indeces)

        for i in range(true_n_sky_points):
             spot_on_sky = d_fitter_fake.spots_on_sky[i]
             target_sky_coord = [c.goodArctan(spot_on_sky[0], spot_on_sky[1]) / deg_to_rad, (np.pi/2-np.arccos(spot_on_sky[2])) / deg_to_rad]
             file_name = 'MCMC_results_RandomRearrangeNumber' + str(rand_index) + '_dipRA' + str(c.round_to_n(target_sky_coord[0], 3)) + 'Dec' + str(c.round_to_n(target_sky_coord[1], 3)) + '.png'
             d_fitter_fake.showMCMCResults(i, n_levels = 20, n_ignore = 0, save_fig = 1, show = 0, file_name = file_name, vars_to_disp = ['dip', 'mono', 'mu_dip', 'mu_mono'] )

        d_fitter_fake.showChiSqrOnSky(save_fig = 1, show_fig = 0, file_name = 'ChiSquareOfDipoleOnSky_fakeData' + str(rand_index) + '.pdf',funct_params = ['dip', 'mono', 'mu_dip', 'mu_mono'])

        end = time.time()
        print ('<<<< TOOK ' + str(end - start) + 's FOR FULL RAND ITERATION >>>> ')

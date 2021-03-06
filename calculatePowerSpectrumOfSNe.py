import matplotlib.pyplot as plt
import numpy as np
import math
from loadSN import loadSN
import cantrips as c
import matplotlib.pyplot as plt
from CosmologicalParameterArchive import CosmologicalParameterArchive  
import scipy.integrate as integrate
import binData as bd 

if __name__ == "__main__":
    
    all_sn = loadSN(1, ['all'], pull_extinctions = 0)
    print ('len(all_sn) = ' + str(len(all_sn))) 
    all_zs = [sn['z'] for sn in all_sn]
    sorted_zs, sorted_sn = c.safeSortOneListByAnother(all_zs, [all_zs, all_sn]) 
    #sorted_muDiffs = [sn['muDiff'] - sn['muDiffWMean'] for sn in sorted_sn]
    #sorted_muErrs = [sn['muErr'] for sn in sorted_sn]
    min_z = np.min(sorted_zs)
    max_z = np.max(sorted_zs)
    #sorted_surveys = [sn['survey'] for sn in sorted_sn]
    #sorted_RAs = [sn['RA'] for sn in sorted_sn]
    #sorted_Decs = [sn['Dec'] for sn in sorted_sn]

    cosmo_arch = CosmologicalParameterArchive()
    OmM = cosmo_arch.getOmegaM()[0]
    OmL = cosmo_arch.getOmegaLambda()[0]
    OmR = cosmo_arch.getOmegaR()[0] 
    H0 = cosmo_arch.getH0()[0]
    speedoflight = cosmo_arch.getc()

    #n_z_bins = 10
    #z_bin_starts = np.linspace(min_z, max_z, n_z_bins+1) [0:-1]
    z_bin_starts = [0.0, 0.3]
    n_z_bins = len(z_bin_starts) 

    for k in range(n_z_bins):
        if k < n_z_bins - 1: z_bin_range = [z_bin_starts[k], z_bin_starts[k+1]]
        else: z_bin_range = [z_bin_starts[k], max_z]
        print ('z_bin_range = ' + str(z_bin_range)) 
        binned_sn = [sn for sn in all_sn if sn['z'] > z_bin_range[0] and sn['z'] < z_bin_range[1]]
        if len(binned_sn) <= 1: break 
        binned_zs = [sn['z'] for sn in binned_sn]
        #binned_muDiffs = [sn['muDiff'] - sn['muDiffWMean'] for sn in binned_sn]
        binned_muDiffs = [sn['muDiff'] for sn in binned_sn]
        binned_muErrs = [sn['muErr'] for sn in binned_sn]
        binned_deltaDLs = [10.0 ** (binned_muDiff / 5.0) - 1 for binned_muDiff in binned_muDiffs]
        binned_deltaDLErrs = [abs(binned_muErrs[i] * 10.0 ** (binned_muDiffs[i] / 5.0) * np.log(10) / 5.0) for i in range(len(binned_deltaDLs))]
        binned_RAs = [sn['RA'] for sn in binned_sn]
        binned_Decs = [sn['Dec'] for sn in binned_sn]
        binned_surveys = [sn['survey'] for sn in binned_sn]
        
        comoving_dist = [speedoflight / H0 * integrate.quad(lambda z_int: 1.0 / math.sqrt( OmM * (1.0 + z_int) ** 3.0 + OmL + OmR * (1.0 + z) ** 4.0 ) , 0 , z)[0] for z in binned_zs ]

        delta_DL_products = c.flattenListOfLists([[binned_deltaDLs[i] * binned_deltaDLs[j] for j in range(i+1, len(binned_deltaDLs)) ] for i in range(len(binned_deltaDLs))] )
        delta_DL_product_errs = c.flattenListOfLists([[ math.sqrt((binned_deltaDLs[i] * binned_deltaDLErrs[j]) ** 2.0 + (binned_deltaDLs[j] * binned_deltaDLErrs[i]) ** 2.0)  for j in range(i+1, len(binned_muDiffs)) ] for i in range(len(binned_muDiffs))] )
        delta_DL_product_norms = np.array(delta_DL_products) / np.array(delta_DL_product_errs)
        delta_mu_product_norms = c.flattenListOfLists([[(binned_muDiffs[i] * binned_muDiffs[j]) / (binned_muErrs[i] * binned_muErrs[j]) for j in range(i+1, len(binned_deltaDLs)) ] for i in range(len(binned_deltaDLs))] )
        ang_seps = c.flattenListOfLists([[c.measureAngularSeparationOnSky([binned_RAs[i], binned_Decs[i]], [binned_RAs[j], binned_Decs[j]]) for j in range(i+1, len(binned_RAs))] for i in range(len(binned_RAs))]) 
    

        #comoving_seps = c.flattenListOfLists([[abs(comoving_dist[i] - comoving_dist[j]) for j in range(i+1, len(comoving_dist)) ] for i in range(len(comoving_dist))],  )
        comoving_centers = c.flattenListOfLists([[ 0.5 * (comoving_dist[i] ** 2.0
                                                          + 2.0 * comoving_dist[i] * comoving_dist[j] * np.cos(c.measureAngularSeparationOnSky([binned_RAs[i], binned_Decs[i]], [binned_RAs[j], binned_Decs[j]]) )
                                                          + comoving_dist[j] ** 2.0) ** 0.5
                                                   for j in range(i+1, len(comoving_dist)) ] for i in range(len(comoving_dist))],  )
        comoving_parallels = [0 for center in comoving_centers]
        pair_num = 0
        for i in range(len(comoving_dist)):
            for j in range(i+1, len(comoving_dist)):
                comoving_parallels[pair_num] = 0.5 * abs(comoving_dist[i] ** 2.0 - comoving_dist[j] ** 2.0) / comoving_centers[pair_num] 
                pair_num = pair_num + 1

        comoving_perps = [0 for center in comoving_centers]
        pair_num = 0
        for i in range(len(comoving_dist)):
            for j in range(i+1, len(comoving_dist)):
                comoving_perps[pair_num] = np.abs(np.sin(ang_seps[pair_num])) * comoving_dist[i] * comoving_dist[j] / comoving_centers[pair_num] 
                pair_num = pair_num + 1

        x_lims = [0.0, 1000.0]
        y_lims = [0.0, 300.0]

        #plt.scatter(delta_DL_products, delta_DL_product_errs)
        #plt.show()

        #poly_for_DL_of_ang_sep = np.polyfit(ang_seps, np.array(delta_DL_products) * 10.0 ** 6.0, 2, w = 1.0 / np.array(delta_DL_product_errs) * 10.0 ** 6.0)
        fit_funct = lambda thetas, A, theta0: A / (1.0 + thetas / theta0)
        theta0s = np.arange(0.5, 10.0, 0.5)
        As = np.linspace(-0.5 * 10.0 ** -5.0, 0.5 * 10.0 ** -5.0, 11)
        fit_params = [[0.0 for theta0 in theta0s] for A in As]
        fit_chisqrs = [[0.0 for theta0 in theta0s] for A in As]
        for i in range(len(As)):
            A = As[i]
            print ('Working on A = ' + str(A)) 
            for j in range(len(theta0s)):
                dof = len(delta_mu_product_norms) - 3 
                theta0 = theta0s[j]
                curve = fit_funct(np.array(ang_seps), A, theta0)
                offset = c.weighted_mean(curve - np.array(delta_mu_product_norms), np.array(np.zeros(np.shape(delta_mu_product_errs)) + 1.0 ))
                r_chi_sqr = np.sum(((curve - offset) - np.array(delta_mu_product_norms)) ** 2.0 / (1 ** 2.0)) / dof 
                fit_params[i][j] = [A, theta0, offset]
                fit_chisqrs[i][j] = r_chi_sqr
        plt.imshow(np.array(fit_chisqrs) - np.mean(fit_chisqrs))
        plt.show()
        min_index = np.argmin(c.flattenListOfLists(fit_chisqrs))
        min_rchisqr =c.flattenListOfLists(fit_chisqrs)[min_index]
        min_fit_params = c.flattenListOfLists(fit_params)[min_index]
        print ('[min_index, min_rchisqr, min_fit_params] = ' + str([min_index, min_rchisqr, min_fit_params])) 
        plt.scatter(ang_seps, delta_mu_product_norms)
        plt.errorbar(ang_seps, delta_mu_product_norms, yerr = 1.0, fmt = 'none')
        sorted_ang_seps = sorted(ang_seps)
        #plt.plot(sorted_ang_seps, np.poly1d(poly_for_DL_of_ang_sep)(sorted_ang_seps) * 10.0 ** -6.0, c = 'r')
        plt.show() 
        #n_ang_sep_bins = 40
        #ang_sep_bin_edges = np.linspace(0.0, 180.0, n_ang_sep_bins+1)
        #ang_sep_bins = [(ang_sep_bin_edges[i] + ang_sep_bin_edges[i+1]) / 2.0 for i in range(n_ang_sep_bins)]
        #delta_DL_product_weights = [1.0 / err ** 2.0 for err in delta_DL_product_errs]
         
        #weighted_dLProducts = np.array(delta_DL_products) * np.array(delta_DL_product_weights)
        #wMean_delta_DL_products = [np.sum(weighted_dLProducts  * ((ang_seps > ang_sep_bin_edges[i]) * (ang_seps < ang_sep_bin_edges[i+1])) )
        #                    / np.sum(np.array(delta_DL_product_weights)  * ((ang_seps > ang_sep_bin_edges[i]) * (ang_seps < ang_sep_bin_edges[i+1])) )
        #                    for i in range(n_ang_sep_bins) ]
        #wMean_delta_DL_product_errs = [ np.sqrt(1 / np.sum(np.array(delta_DL_product_weights)  * ((ang_seps > ang_sep_bin_edges[i]) * (ang_seps < ang_sep_bin_edges[i+1])) ))
        #                                   for i in range(n_ang_sep_bins) ]     #/ np.sum((ang_seps > ang_sep_bin_edges[i]) * (ang_seps < ang_sep_bin_edges[i+1])))
                                        
        
        #print ('[ang_sep_bins, wMean_delta_DL_products, wMean_delta_DL_product_errs] = ' + str([ang_sep_bins, wMean_delta_DL_products, wMean_delta_DL_product_errs])) 
        #ang_sep_bins, binned_delta_DL_products = bd.binData(ang_seps, delta_DL_products, y_errs = delta_DL_product_errs)
        
        #plt.errorbar( ang_sep_bins, np.array(wMean_delta_DL_products) * 10.0 ** 6.0, yerr = np.array(wMean_delta_DL_product_errs)* 10.0 ** 6.0, ecolor = 'r', fmt='none')
        #plt.plot( ang_sep_bins, np.array(wMean_delta_DL_products) * 10.0 ** 6.0, c = 'r')
        #plt.show()

        #For temporary reference: z_bin_range = [0.0, 0.2], [ang_sep_bins, binned_delta_DL_products ] = [[4.496745300539829, 13.490192737092542, 22.483640173645256, 31.477087610197966, 40.470535046750676, 49.46398248330339, 58.4574299198561, 67.45087735640881, 76.44432479296152, 85.43777222951424, 94.43121966606695, 103.42466710261967, 112.41811453917236, 121.41156197572508, 130.4050094122778, 139.3984568488305, 148.39190428538322, 157.38535172193593, 166.37879915848865, 175.37224659504136], [[-4.260644868192591e-05, 2.872905605678743e-05, -4.9648795263294904e-05, 8.153467709331572e-05, 2.3297416189716368e-05, -1.5604403094977693e-05, -3.80722209042677e-06, 6.427009011841985e-05, 0.000101907167215507, -3.298702744098934e-05, 0.00010212855619122636, 1.617629813069657e-05, -1.3001362164841915e-05, -5.802560118321389e-05, 8.697676130135586e-05, -7.356110090141447e-06, -2.2866314911848203e-05, 2.2562646881793654e-05, 4.31581375580043e-05, 4.593118205280416e-05], [5.191209018625268e-06, 5.1926132496077226e-06, 4.845107954428657e-06, 4.691414776782246e-06, 5.289313732762643e-06, 5.22683084575302e-06, 5.253898755459909e-06, 5.657720857108593e-06, 5.547187826423691e-06, 5.85592073441274e-06, 5.5098522342895e-06, 5.392698165481146e-06, 6.024438178221094e-06, 5.507375424106925e-06, 5.748565069419812e-06, 7.268802920442577e-06, 7.488401526126141e-06, 8.82550441316428e-06, 9.435344138465517e-06, 1.4544937105615466e-05], [3838, 4592, 5382, 5964, 5990, 5228, 4694, 4281, 4413, 4327, 4608, 4673, 5038, 6018, 4702, 3251, 2919, 2153, 1456, 728]]]

        max_comoving_perps = max(comoving_perps)
        min_comoving_perps = min(comoving_perps)
        max_comoving_parallels = max(comoving_parallels)
        min_comoving_parallels = min(comoving_parallels)
        min_DL_products = min(delta_mu_product_norms)
        max_DL_products = max(delta_mu_product_norms)
        print ('[max_comoving_perps, max_comoving_parallels] = ' + str([max_comoving_perps, max_comoving_parallels])) 
        sc = plt.scatter(comoving_perps, comoving_parallels, c = delta_mu_product_norms, marker = '.', s = 0.5)
        plt.xlabel(r'$r_{\perp}$' + ' [Mpc]')
        plt.ylabel(r'$r_{\parallel}$' + ' [Mpc]')
        plt.title('Uncertainty-normalized luminosity distance correlations for all SNe pairs')
        cbar = plt.colorbar(sc)
        cbar.set_label(r'$(\delta L_i \times \delta L_j)/(\sigma_{\delta L \times \delta L})$', rotation=270)
        #plt.xlim(x_lims)
        #plt.ylim(y_lims) 
        plt.show()



        fitting_funct = lambda rperp, rpar, A, muperp, mupar, sig, el: A * np.exp( -((rperp - muperp) ** 2.0 / (2.0 * sig ** 2.0) + (rpar - mupar) ** 2.0 / (2.0 * (sig * el) ** 2.0)) )
        As = np.linspace(-max(-min_DL_products, max_DL_products), max(-min_DL_products, max_DL_products), 3)
        print ('As = ' + str(As)) 
        muperps = np.linspace(0.0, max_comoving_perps, 3)
        muparallels = np.linspace(0.0, max_comoving_parallels, 3)
        sigs = np.linspace(20.0, 1000.0, 6)
        els = [1.0]

        delta_DL_product_fits = [ [ [ [ [[] for m in range(len(els)) ] for l in range(len(sigs)) ] for k in range(len(muparallels)) ] for j in range(len(muperps)) ] for i in range(len(As)) ]
        delta_DL_rchisqrs = [ [ [ [ [0.0 for m in range(len(els)) ] for l in range(len(sigs)) ] for k in range(len(muparallels)) ] for j in range(len(muperps)) ] for i in range(len(As)) ]
        dof = len(delta_DL_products) - 6 
        for i in range(len(As)):
            A = As[i]
            print ('A = ' + str(A) ) 
            for j in range(len(muperps)):
                muperp = muperps[j]
                print ('muperp = ' + str(muperp) ) 
                for k in range(len(muparallels)):
                    muparallel = muparallels[k]
                    #print ('muparallel= ' + str(muparallel) )
                    for l in range(len(sigs)):
                        sig = sigs[l]
                        for m in range(len(els)):
                            el = els[m]
                            funct_output = fitting_funct(np.array(comoving_perps), np.array(comoving_parallels), A, muperp, muparallel, sig, el) 
                            shift = c.weighted_mean(funct_output - np.array(delta_DL_products), np.array(delta_DL_product_errs))
                            
                            rchisqr = np.sum((funct_output - shift) - np.array(delta_DL_products) ** 2.0 / (np.array(delta_DL_product_errs) ** 2.0)) / dof
                            delta_DL_product_fits[i][j][k][l][m] = [A, muperp, muparallel, el, sig, shift]
                            delta_DL_rchisqrs[i][j][k][l][m] = rchisqr
                            
        print ('delta_DL_product_fits = ' + str(delta_DL_product_fits))
        print ('delta_DL_rchisqrs = ' + str(delta_DL_rchisqrs))
        flattened_rchisqr = c.flattenListOfLists(c.flattenListOfLists(c.flattenListOfLists(c.flattenListOfLists(delta_DL_rchisqrs))))
        flattened_fits = c.flattenListOfLists(c.flattenListOfLists(c.flattenListOfLists(c.flattenListOfLists(delta_DL_product_fits))))
        min_rchisqr_arcmin = np.argmin(flattened_rchisqr)
        min_rchisqr = flattened_rchisqr[min_rchisqr_arcmin]
        min_params = flattened_fits[min_rchisqr_arcmin]
        min_params = [0.00002, 1000.0, 100.0, 500, 1.0, 0.0] 
        print ('[min_rchisqr_arcmin, min_rchisqr, min_params] = ' + str([min_rchisqr_arcmin, min_rchisqr, min_params])) 
        
        max_perp = np.max(comoving_perps) #1000
        max_parallel = np.max(comoving_parallels) #300
        n_comoving_bins = [100, int(100 * max_parallel / max_perp)]
        #print ('n_comoving_bins = ' + str(n_comoving_bins)) 
        comoving_bin_size = [(max_perp - 0.0 ) / n_comoving_bins[0], (max_parallel - 0.0 ) / n_comoving_bins[1]]
        print ('comoving_bin_size = ' + str(comoving_bin_size)) 
        perp_bin_edges = (np.arange(0.0, max_perp, comoving_bin_size[0])).tolist() 
        parallel_bin_edges = (np.arange(0.0, max_parallel, comoving_bin_size[1])).tolist() 
        #print ('[max_perp, max_parallel] = ' + str([max_perp, max_parallel])) 
        #print ('perp_bin_edges = ' + str(perp_bin_edges))
        #print ('parallel_bin_edges = ' + str(parallel_bin_edges)) 
        if parallel_bin_edges[-1] < max_parallel: parallel_bin_edges = parallel_bin_edges + [max_parallel] 
        if perp_bin_edges[-1] < max_perp: perp_bin_edges = perp_bin_edges + [max_perp]
        perp_mesh, parallel_mesh = np.meshgrid([(perp_bin_edges[i] + perp_bin_edges[i+1])/2.0 for i in range(len(perp_bin_edges)-1)],
                                               [(parallel_bin_edges[i] + parallel_bin_edges[i+1])/2.0 for i in range(len(parallel_bin_edges)-1)])
        print ('[perp_mesh, parallel_mesh] = ' + str([perp_mesh, parallel_mesh])) 

        #Bin in such a way that points not in any bin have an index of -1 
        covariance_binning = [[sum([j if (comoving_perps[k] >= perp_bin_edges[j-1] and comoving_perps[k] < perp_bin_edges[j]) else 0 for j in range(1, len(perp_bin_edges)) ]) - 1,
                               sum([i if (comoving_parallels[k] >= parallel_bin_edges[i-1] and comoving_parallels[k] < parallel_bin_edges[i]) else 0 for i in range(1, len(parallel_bin_edges)) ]) - 1]
                              for k in range(len(comoving_perps))]
        #print ('covariance_binning = ' + str(covariance_binning)) 
        deltaDLs_inBins = [[[] for j in range(len(perp_bin_edges)-1) ] for i in range(len(parallel_bin_edges)-1)]
        deltaDLsErrs_inBins = [[[] for j in range(len(perp_bin_edges)-1) ] for i in range(len(parallel_bin_edges)-1)]
        deltaDLsNorms_inBins = [[[] for j in range(len(perp_bin_edges)-1) ] for i in range(len(parallel_bin_edges)-1)]
        print ('I am binning the measured covariances...') 
        for k in range(len(covariance_binning)):
            if k % 10000 == 0: print (str(k / len(covariance_binning) * 100) + '% done.')  
            perp_bin, parallel_bin = covariance_binning[k]
            if perp_bin >= 0 and parallel_bin >= 0: 
            #print ('[perp_bin, parallel_bin] = ' + str([perp_bin, parallel_bin]))
                deltaDL = delta_DL_products[k]
                deltaDLErr = delta_DL_product_errs[k]
                deltaDLNorm = delta_DL_product_norms[k]
                deltaDLs_inBins[parallel_bin][perp_bin] = deltaDLs_inBins[parallel_bin][perp_bin] + [deltaDL]
                deltaDLsErrs_inBins[parallel_bin][perp_bin] = deltaDLsErrs_inBins[parallel_bin][perp_bin] + [deltaDLErr]
                deltaDLsNorms_inBins[parallel_bin][perp_bin] = deltaDLsNorms_inBins[parallel_bin][perp_bin] + [deltaDLNorm]
        print ('I am done binning the measured covariances.')
        smoothedAutocorr = [[c.weighted_mean(deltaDLs_inBins[i][j], deltaDLsErrs_inBins[i][j])  for j in range(len(perp_bin_edges)-1) ] for i in range(len(parallel_bin_edges)-1)]
        smoothedAutocorr_errs = [[math.sqrt(len(deltaDLs_inBins[i][j]) / sum([err ** -2.0 for err in deltaDLsErrs_inBins[i][j]])) /np.sqrt(len(deltaDLs_inBins[i][j]))  if len(deltaDLs_inBins[i][j]) >=1 else 0.0 for j in range(len(perp_bin_edges)-1) ] for i in range(len(parallel_bin_edges)-1)] #standard error calculation
        smoothedAutocorr_norms = [[np.mean(deltaDLsNorms_inBins[i][j])  for j in range(len(perp_bin_edges)-1) ] for i in range(len(parallel_bin_edges)-1)]
        #smoothedAutocorr_norms = [[np.mean(deltaDLsNorms_inBins[i][j]) * np.sqrt(len(deltaDLsNorms_inBins[i][j])) for j in range(len(perp_bin_edges)-1) ] for i in range(len(parallel_bin_edges)-1)]
        smoothedAutocorrDeviations = np.array(smoothedAutocorr) / np.array(smoothedAutocorr_errs)

        fitted_dl_product = fitting_funct(perp_mesh, parallel_mesh, *(min_params[0:-1])) - min_params[-1]
        n_levels = 21
        levels = (np.linspace(0, np.max(fitted_dl_product), n_levels)).tolist()
        f, axarr = plt.subplots(2,1, squeeze = False, figsize = [4, 8])
        #axarr[0,2].imshow(smoothedAutocorr, vmin = levels[0], vmax = levels[-1])
        #axarr[0,1].imshow(smoothedAutocorr_errs)
        print ('[np.nanmean(smoothedAutocorr_norms) - 0.5 * np.nanstd(smoothedAutocorr_norms), np.nanmean(smoothedAutocorr_norms) + 0.5 * np.nanstd(smoothedAutocorr_norms)] = ' + str([np.nanmean(smoothedAutocorr_norms) - 0.5 * np.nanstd(smoothedAutocorr_norms), np.nanmean(smoothedAutocorr_norms) + 0.5 * np.nanstd(smoothedAutocorr_norms)])) 
        im1 = axarr[0,0].imshow(smoothedAutocorr_norms, vmin = np.nanmean(smoothedAutocorr_norms) - 0.5 * np.nanstd(smoothedAutocorr_norms), vmax = np.nanmean(smoothedAutocorr_norms) + 0.5 * np.nanstd(smoothedAutocorr_norms) )
        axarr[0,0].set_xlabel(r'$r_{\perp}$ [Mpc]')
        axarr[0,0].set_ylabel(r'$r_{\parallel}$ [Mpc]')
        axarr[0,0].set_xticks(np.linspace(0, n_comoving_bins[0], 6) )
        axarr[0,0].set_yticks(np.linspace(0, n_comoving_bins[1], 6) )
        axarr[0,0].set_xticklabels([c.round_to_n(tick, 3) for tick in axarr[0,0].get_xticks() * comoving_bin_size[0]])
        axarr[0,0].set_yticklabels( [c.round_to_n(tick, 3) for tick in axarr[0,0].get_yticks() * comoving_bin_size[0]] )
        print ('[axarr[0,0].get_xticks(), axarr[0,0].get_yticks()] = ' + str([axarr[0,0].get_xticks(), axarr[0,0].get_yticks()])) 
        #axarr[0,0].imshow(smoothedAutocorrDeviations)
        n_contours = 11
        #axarr[0].contourf(perp_mesh, parallel_mesh, smoothedAutocorr, levels = [np.min(smoothedAutocorr)] + (np.linspace(np.median([corr for corr in c.flattenListOfLists(smoothedAutocorr) if corr < 0.0]), np.median([corr for corr in c.flattenListOfLists(smoothedAutocorr) if corr > 0.0]), n_contours-2 )).tolist() + [np.max(smoothedAutocorr)])
        #im0 = axarr[1,2].contourf(perp_mesh, parallel_mesh, smoothedAutocorr, levels = levels)
        #axarr[1,1].contourf(perp_mesh, parallel_mesh, smoothedAutocorr_errs, levels = n_contours) # np.linspace(0.0, 2.0 * np.median([err for err in c.flattenListOfLists(smoothedAutocorr_errs)]), n_contours ))
        axarr[1,0].contourf(perp_mesh, parallel_mesh, smoothedAutocorr_norms, levels = np.linspace(np.nanmean(smoothedAutocorr_norms) - 0.5 * np.nanstd(smoothedAutocorr_norms), np.nanmean(smoothedAutocorr_norms) + 0.5 * np.nanstd(smoothedAutocorr_norms), 21))
        axarr[1,0].set_xlabel(r'$r_{\perp}$ [Mpc]')
        axarr[1,0].set_ylabel(r'$r_{\parallel}$ [Mpc]')
        f.suptitle(r'Mean $(\delta L_i \times \delta L_j)/(\sigma_{\delta L \times \delta L})$ in bins')
        #im1 = axarr[1,0].contourf(perp_mesh, parallel_mesh, smoothedAutocorrDeviations) #, levels = np.linspace(np.median([corr for corr in c.flattenListOfLists(smoothedAutocorrDeviations) if corr < 0.0]), np.median([corr for corr in c.flattenListOfLists(smoothedAutocorrDeviations) if corr > 0.0]), n_contours ))
        print ('min_params = ' + str(min_params)) 
        #axarr[2,1].contourf(perp_mesh, parallel_mesh, fitted_dl_product, levels = levels) 
        f.subplots_adjust(right=0.7)
        cbar_ax = f.add_axes([0.75, 0.15, 0.05, 0.7])
        cbar = f.colorbar(im1, cax=cbar_ax)
        cbar.set_label(r'Mean $(\delta L_i \times \delta L_j)/(\sigma_{\delta L \times \delta L})$ in bin', rotation=270, labelpad=20)
        
        plt.show()
        

    #plt.scatter(ang_seps, mu_diff_products)
    #plt.errorbar(ang_seps, mu_diff_products, yerr = mu_diff_product_errs, fmt = 'none')
    #plt.show() 
                 

    #plt.scatter(comoving_seps, ang_seps, c = mu_diff_products, marker = '.')
    #plt.show() 

    
    

    #plt.hist(mu_diff_products, bins= 200)
    #plt.show()
    #plt.hist(comoving_seps, bins= 200)
    #plt.show()
    #plt.hist(ang_seps, bins= 1000)
    #plt.show() 


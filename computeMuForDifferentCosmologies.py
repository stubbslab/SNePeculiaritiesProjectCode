import math
import numpy as np
from CosmologicalParameterArchive import CosmologicalParameterArchive
import scipy.integrate as integrate

def unnormalizedGauss(z, A, mu, sigma):
    return A * math.exp(-(z - mu)**2.0 / (2.0 * sigma**2.0))

def lumDistCanonical(z, H0 = 70.0, OmM = 0.3, OmL = 0.7, Om0 = 1.0, use_canon_params = 1):
    #dp = c * integr.quad(1 / a(t) ) 
    cosmo_archive = CosmologicalParameterArchive()
    if use_canon_params:
        H0 = cosmo_archive.getH0()[0]
        OmM = cosmo_archive.getOmegaM()[0]
        OmL = cosmo_archive.getOmegaLambda()[0]
        Om0 = cosmo_archive.getOmega0()[0]
    c = cosmo_archive.getc()
    #print 'Om0 = ' + str(Om0)
    #print type(Om0) 
    #print 'OmL = ' + str(OmL)
    #print type(OmL) 
    #print 'OmM = ' + str(OmM)
    #print type(OmM) 
    #print 'z = ' + str(z)
    #print type(z) 
    scaling = c / H0 #unite of Mpc 
    dl = scaling * (1.0+z) * integrate.quad(lambda x: 1.0 / math.sqrt( OmM * ((1.0+x) ** (3.0)) + OmL + (1.0 - Om0) * ((1.0+x) ** (2.0)) ), 0.0, z)[0]
    #print 'dl = ' + str(dl) 
    #print 'for z = ' + str(z) + ', dl = ' + str(dl) 
    return dl 

def lumDistOscilDE1(z, A, omega, phi, H0 = 70.0, OmM = 0.3, OmL = 0.7, Om0 = 1.0, use_canon_params = 1):
    cosmo_archive = CosmologicalParameterArchive()
    c = cosmo_archive.getc()
    scaling = c / H0 #unite of MpC
    aML = (OmM / (1.0 - OmM)) ** (1.0 / 3.0)
    a = 1.0 /(1.0 + z) 
    t = 1.0 / H0 * 2.0 / (3.0 * math.sqrt(1.0 - OmM)) * np.log( (a / aML)**(3.0/2.0) + math.sqrt(1.0 + (a / aML)**(3.0/2.0)))
    dl = scaling * (1.0+z) * integrate.quad(lambda x: 1.0 / math.sqrt( OmM * ((1.0+x) ** (3.0)) + OmL / (1.0 - A * np.sin(phi) ) * (1.0 - A * np.sin(omega * t + phi) ) + (1.0 - Om0) * ((1+x) ** (2.0)) ), 0.0, z)[0]
    return dl 

def lumDistHillValleyHillMatter(z, A, rat, mu, sigma, H0 = 70.0, OmM = 0.3, OmL = 0.7, Om0 = 1.0):
    cosmo_archive = CosmologicalParameterArchive()
    c = cosmo_archive.getc()
    scaling = c / H0 #unite of MpC
    matter_perturbation = lambda z_var, A_var, rat_var, mu_var, sigma_var:( unnormalizedGauss(z_var,rat_var * A_var, mu_var + 2.0 * sigma_var, sigma_var) 
                                                                            + unnormalizedGauss(z_var, A_var, mu_var, sigma_var)  
                                                                            + unnormalizedGauss(z_var,rat_var * A_var, mu_var - 2.0 * sigma_var, sigma_var) )
    #print 'matter_perturbation = ' + str(matter_perturbation (z, A, B, mu, sigma))
    dl = scaling * (1.0+z) * integrate.quad(lambda x: 1.0 / math.sqrt( OmM * ((1.0 + x) ** 3.0)
                                                                       * (1.0 + matter_perturbation (x, A, rat, mu, sigma)) /(1.0 + matter_perturbation (0.0, A, rat, mu, sigma))
                                                                       + OmL + (1.0 - Om0) * ((1+x) ** (2.0)) ), 0.0, z)[0]
    #print 'dl = ' + str(dl)

    return dl

def lumDistHillValleyMatter(z, A, mu1, mu2, sigma1, sigma2, H0 = 70.0, OmM = 0.3, OmL = 0.7, Om0 = 1.0):
    cosmo_archive = CosmologicalParameterArchive()
    c = cosmo_archive.getc()
    scaling = c / H0 #unite of MpC
    matter_perturbation = lambda z_var, A_var, mu1_var, mu2_var, sigma1_var, sigma2_var: ( unnormalizedGauss(z_var, A_var, mu1_var, sigma1_var)
                                                                                           + unnormalizedGauss(z_var, -1.0  *  A_var * sigma1_var / sigma2_var, mu2_var, sigma2_var)
                                                                                           if unnormalizedGauss(z_var, mu1_var, sigma1_var)
                                                                                           + unnormalizedGauss(z_var, -1.0  *  A_var * sigma1_var / sigma2_var, mu2_var, sigma2_var) > -1.0
                                                                                           else 1.0 ) 
    #print 'matter_perturbation = ' + str(matter_perturbation (z, A, B, mu, sigma))
    if matter_perturbation (z, A, mu1, mu2, sigma1, sigma2) <=-1.0:
        print ('matter purturbation at ' + str([])) 
    dl = scaling * (1.0+z) * integrate.quad(lambda x: 1.0 / math.sqrt( OmM * ((1.0 + x) ** 3.0)
                                                                       * (1.0 + matter_perturbation (x, A, mu1, mu2, sigma1, sigma2)) /(1.0 + matter_perturbation (0.0, A, mu1, mu2, sigma1, sigma2))
                                                                       + OmL + (1.0 - Om0) * ((1+x) ** (2.0)) ), 0.0, z)[0]
    #print 'dl = ' + str(dl)

    return dl

def lumDistHillValleyDE(z, A, mu1, mu2, sigma1, sigma2, H0 = 70.0, OmM = 0.3, OmL = 0.7, Om0 = 1.0):
    #print 'entering DE hillValley function' 
    cosmo_archive = CosmologicalParameterArchive()
    c = cosmo_archive.getc()
    scaling = c / H0 #unite of MpC
    matter_perturbation = lambda z_var, A_var, mu1_var, mu2_var, sigma1_var, sigma2_var:  ( unnormalizedGauss(z_var, A_var, mu1_var, sigma1_var)
                                                                                           + unnormalizedGauss(z_var, -1.0  *  A_var * sigma1_var / sigma2_var, mu2_var, sigma2_var)
                                                                                           if unnormalizedGauss(z_var, mu1_var, sigma1_var)
                                                                                           + unnormalizedGauss(z_var, -1.0  *  A_var * sigma1_var / sigma2_var, mu2_var, sigma2_var) > -1.0
                                                                                           else 1.0 )  
    #print 'matter_perturbation = ' + str(matter_perturbation (z, A, B, mu, sigma))
    #if matter_perturbation (z, A, mu1, mu2, sigma1, sigma2) <=-1.0:
    #    print 'matter purturbation at ' + str([])
    dl = scaling * (1.0+z) * integrate.quad(lambda x: 1.0 / math.sqrt( OmM * ((1.0 + x) ** 3.0) 
                                                                       + OmL * (1.0 + matter_perturbation (x, A, mu1, mu2, sigma1, sigma2)) /(1.0 + matter_perturbation (0.0,A,mu1,mu2,sigma1,sigma2))
                                                                       + (1.0 - Om0) * ((1+x) ** (2.0)) ), 0.0, z)[0]
    #print 'dl = ' + str(dl)

    return dl

def lumDistGaussDEandMatter(z, A_DE, A_mat, mu_DE, mu_mat, sigma_DE, sigma_mat, H0 = 70.0, OmM = 0.3, OmL = 0.7, Om0 = 1.0):
    #print 'entering DE hillValley function' 
    cosmo_archive = CosmologicalParameterArchive()
    c = cosmo_archive.getc()
    scaling = c / H0 #unite of MpC
    gauss_perturbation = lambda z_var, A_var, mu_var, sigma_var:(  unnormalizedGauss(z_var, A_var, mu_var, sigma_var)
                                                                    if unnormalizedGauss(z_var, A_var, mu_var, sigma_var) > -1.0
                                                                                           else -1.0 ) 
    #print 'matter_perturbation = ' + str(matter_perturbation (z, A, B, mu, sigma))
    #if matter_perturbation (z, A, mu1, mu2, sigma1, sigma2) <=-1.0:
    #    print 'matter purturbation at ' + str([])
    dl = scaling * (1.0+z) * integrate.quad(lambda x: 1.0 / math.sqrt( OmM * ((1.0 + x) ** 3.0) * (1.0 + gauss_perturbation (x, A_mat, mu_mat, sigma_mat)) /(1.0 + gauss_perturbation (0.0, A_mat, mu_mat, sigma_mat))
                                                                       + OmL * (1.0 + gauss_perturbation (x, A_DE, mu_DE, sigma_DE)) /(1.0 + gauss_perturbation (0.0, A_DE, mu_DE, sigma_DE))
                                                                       + (1.0 - Om0) * ((1+x) ** (2.0)) ), 0.0, z)[0]
    #print 'dl = ' + str(dl)

    return dl

def lumDistConstDEandMatter(z, A_DE, A_mat, start_DE, start_mat, width_DE, width_mat, H0 = 70.0, OmM = 0.3, OmL = 0.7, Om0 = 1.0):
    #print 'entering DE hillValley function' 
    cosmo_archive = CosmologicalParameterArchive()
    c = cosmo_archive.getc()
    scaling = c / H0 #unite of MpC
    const_perturbation = lambda z_var, A_var, start_var, width_var:(  A_var if (z_var < start_var and z_var > start_var - width_var) else 0.0) 
                                                                    
    #print 'matter_perturbation = ' + str(matter_perturbation (z, A, B, mu, sigma))
    #if matter_perturbation (z, A, mu1, mu2, sigma1, sigma2) <=-1.0:
    #    print 'matter purturbation at ' + str([])
    dl = scaling * (1.0+z) * integrate.quad(lambda x: 1.0 / math.sqrt( OmM * ((1.0 + x) ** 3.0) * (1.0 + const_perturbation (x, A_mat, start_mat, width_mat)) /(1.0 + const_perturbation (0.0, A_mat, start_mat, width_mat))
                                                                       + OmL * (1.0 + const_perturbation (x, A_DE, start_DE, width_DE)) /(1.0 + const_perturbation (0.0, A_DE, start_DE, width_DE))
                                                                       + (1.0 - Om0) * ((1+x) ** (2.0)) ), 0.0, z)[0]
    #print 'dl = ' + str(dl)

    return dl

def perturbativeFunction (z, params, function, funct_name = '', min_val = -np.inf, max_val = np.inf, ):
    perturb_value = function(z, *params)
    if perturb_value < min_val:
        print ('Warning: ' + funct_name + ' below min allowed value of ' + str(min_val) + ' at params: ' + str([z] + params) + '. Setting to min_val.') 
        perturb_value = min_val
    if perturb_value > max_val:
        print ('Warning: ' + funct_name + ' above max allowed value of ' + str(max_val) + ' at params: ' + str([z] + params) + '. Setting to max_val.') 
        perturb_value = min_val
    return perturb_value

def lumDistFromSpecifiedFunction(z, extra_params, matter_function, de_function, other_function, matter_param_indeces, de_param_indeces, other_param_indeces,
                                                  H0 = 70.0, OmM = 0.3, OmL = 0.7, Om0 = 1.0, use_canon_params = 1 ):
    cosmo_archive = CosmologicalParameterArchive()
    if use_canon_params:
        H0 = cosmo_archive.getH0()[0]
        OmM = cosmo_archive.getOmegaM()[0]
        OmL = cosmo_archive.getOmegaLambda()[0]
        Om0 = cosmo_archive.getOmega0()[0]
    c = cosmo_archive.getc()
    scaling = c / H0 #unite of MpC
    if matter_function is None:
        matter_params = []
        matter_perturbation = lambda z_var, matter_params_var: 0.0
    else: 
        matter_params = [extra_params[i] for i in matter_param_indeces]
        matter_perturbation = lambda z_var, matter_params_var: perturbativeFunction (z_var, matter_params_var, matter_function, funct_name = 'matter perturbation', min_val = -1.0) 
    if de_function is None:
        de_params = []
        de_perturbation = lambda z_var, de_params_var: 0.0
    else:
        de_params = [extra_params[i] for i in de_param_indeces]
        de_perturbation = lambda z_var, de_params_var: perturbativeFunction(z_var, de_params_var, de_function, funct_name = 'de perturbation', min_val = -1.0) 
    if other_function is None:
        other_params = []
        other_perturbation = lambda z_var, other_params_var: 0.0
    else:
        other_params = [extra_params[i] for i in other_param_indeces]
        other_perturbation = lambda z_var, other_params_var: perturbativeFunction(z_var, other_params_var, other_function, funct_name = 'extra perturbation') 
    
    dl = scaling * (1.0+z) * integrate.quad(lambda x: 1.0 / math.sqrt( OmM * ((1.0 + x) ** 3.0)
                                                                       * (1.0 + matter_perturbation (x, matter_params)) / (1.0 + matter_perturbation (0.0, matter_params))
                                                                       + OmL
                                                                       * (1.0 + de_perturbation (x, de_params)) / (1.0 + de_perturbation (0.0, de_params))
                                                                       + other_perturbation (x, other_params)  
                                                                       + (1.0 - Om0) * ((1+x) ** (2.0)) ), 0.0, z)[0]    
    return dl 
    
def computeMuForCosmology(z, extra_params = [], H0 = 70.0, OmM = 0.3, OmL = 0.7, Om0 = 1.0, use_canon_params = 1, pre_loaded_function = 'canon',
                          matter_function = None, matter_param_indeces = [], de_function = None, de_param_indeces = [], other_function = None, other_param_indeces = []):
    cosmo_archive = CosmologicalParameterArchive()
    if use_canon_params:
        H0 = cosmo_archive.getH0()[0]
        OmM = cosmo_archive.getOmegaM()[0]
        OmL = cosmo_archive.getOmegaLambda()[0]
        Om0 = cosmo_archive.getOmega0()[0]

    if matter_function is None and de_function is None and other_function is None:
        if pre_loaded_function == 'oscillitory_DE_1':
            lum_dist_funct = lambda z, extra_params: lumDistOscilDE1(z, *extra_params, H0=H0, OmM = OmM, OmL = OmL, Om0 = Om0)
        
        elif pre_loaded_function == 'up_down_up_matter' or pre_loaded_function == 'down_up_down_matter' or pre_loaded_function == 'hill_valley_hill_matter' or pre_loaded_function == 'valley_hill_valley_matter':
            lum_dist_funct = lambda z, extra_params: lumDistHillValleyHillMatter(z, *extra_params, H0=H0, OmM = OmM, OmL = OmL, Om0 = Om0)

        elif pre_loaded_function == 'up_down_matter' or pre_loaded_function == 'down_up_matter' or pre_loaded_function == 'hill_valley_matter' or pre_loaded_function == 'valley_hill_matter':
            lum_dist_funct = lambda z, extra_params: lumDistHillValleyMatter(z, *extra_params, H0=H0, OmM = OmM, OmL = OmL, Om0 = Om0)

        elif pre_loaded_function == 'up_down_DE' or pre_loaded_function == 'down_up_DE' or pre_loaded_function == 'hill_valley_DE' or pre_loaded_function == 'valley_hill_DE':
            lum_dist_funct = lambda z, extra_params: lumDistHillValleyDE(z, *extra_params, H0=H0, OmM = OmM, OmL = OmL, Om0 = Om0)

        elif pre_loaded_function == 'gauss' or pre_loaded_function == 'hill' or pre_loaded_function == 'valley':
            lum_dist_funct = lambda z, extra_params: lumDistGaussDEandMatter(z, *extra_params, H0=H0, OmM = OmM, OmL = OmL, Om0 = Om0)

        elif pre_loaded_function == 'const':
            lum_dist_funct = lambda z, extra_params: lumDistConstDEandMatter(z, *extra_params, H0=H0, OmM = OmM, OmL = OmL, Om0 = Om0)
        
        else:
            lum_dist_funct = lambda z, extra_params: lumDistCanonical(z, OmM = OmM, OmL = OmL, H0 = H0, Om0 = Om0, use_canon_params = use_canon_params)

    else:
        lum_dist_funct = lambda z, extra_params: lumDistFromSpecifiedFunction (z, extra_params, matter_function, de_function, other_function,
                                                                               matter_param_indeces, de_param_indeces, other_param_indeces,
                                                                               H0 = H0, OmM = OmM, OmL = OmL, Om0 = Om0, use_canon_params = use_canon_params )

    if type(z) is int or type(z) is float:
        lum_dists = lum_dist_funct(z, extra_params)
    else:
        lum_dists = [lum_dist_funct(z_elem, extra_params) for z_elem in z]
    #print 'lum_dists = ' + str(lum_dists) 
    mu = 5.0 * np.log10(lum_dists) + 25.0 # distance modulus 
    return mu 

def computeMuResidualsForCosmologies(z,
                                      extra_params_1 = [], H0_1 = 70.0, OmM_1 = 0.3, OmL_1 = 0.7, Om0_1 = 1.0, use_canon_params_1 = 1, pre_loaded_function_1 = 'canon',
                                      matter_function_1 = None, matter_param_indeces_1 = [], de_function_1 = None, de_param_indeces_1 = [], other_function_1 = None, other_param_indeces_1 = [],
                                      extra_params_2 = [], H0_2 = 70.0, OmM_2 = 0.3, OmL_2 = 0.7, Om0_2 = 1.0, use_canon_params_2 = 1, pre_loaded_function_2 = 'canon',
                                      matter_function_2 = None, matter_param_indeces_2 = [], de_function_2 = None, de_param_indeces_2 = [], other_function_2 = None, other_param_indeces_2 = []):
    return ( computeMuForCosmology(z,
                                   extra_params = extra_params_1, H0 = H0_1, OmM = OmM_1, OmL = OmL_1, Om0 = Om0_1, use_canon_params = use_canon_params_1, pre_loaded_function = pre_loaded_function_1,
                                   matter_function = matter_function_1, matter_param_indeces = matter_param_indeces_1,
                                   de_function = de_function_1, de_param_indeces = de_param_indeces_1,
                                   other_function = other_function_1, other_param_indeces = other_param_indeces_1)
             - computeMuForCosmology(z,
                                   extra_params = extra_params_2, H0 = H0_2, OmM = OmM_2, OmL = OmL_2, Om0 = Om0_2, use_canon_params = use_canon_params_2, pre_loaded_function = pre_loaded_function_2,
                                   matter_function = matter_function_2, matter_param_indeces = matter_param_indeces_2,
                                   de_function = de_function_2, de_param_indeces = de_param_indeces_2,
                                   other_function = other_function_2, other_param_indeces = other_param_indeces_2) )

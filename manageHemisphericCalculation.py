import math
from loadSN import loadSN
import numpy as np
from AstronomicalParameterArchive import AstronomicalParameterArchive
from randomSimulationFunctions import randomSortData 
from randomSimulationFunctions import simulateNormallyDistributedData

def computeHemisphericResidualForObservations(data_set, axes,
                                              surveys_of_interest = ['all'], pull_extinctions = 1,
                                              opening_angle = math.pi / 2.0, statistic_funct = 'wmean', statistic_err_funct = 'wmean'):
    all_sns = loadSN(data_set, ['all'], pull_extinctions = pull_extinctions)

    if surveys_of_interest[0].lower() == 'all':
        surveys_of_interest = np.unique([sn['survey'] for sn in all_sns]).tolist()

    sn_of_interest = [sn for sn in all_sns if sn['survey'] in surveys_of_interest]

    return [computeHemisphericResidual(sn_of_interest, axis, opening_angle = opening_angle, statistic_funct = statistic_funct, statistic_err_funct = statistic_err_funct) for axis in axes]

def computeHemisphericResidualForSimulations(data_set, axes,
                                             surveys_of_interest = ['all'], pull_extinctions = 1,
                                             opening_angle = math.pi / 2.0, statistic_funct = 'wmean', statistic_err_funct = 'wmean',
                                             simulation_method = 'random_sort'):
    all_sns = loadSN(data_set, ['all'], pull_extinctions = pull_extinctions)

    all_sns = loadSN(data_set, ['all'], pull_extinctions = pull_extinctions)

    if surveys_of_interest[0].lower() == 'all':
        surveys_of_interest = np.unique([sn['survey'] for sn in all_sns]).tolist()

    sn_of_interest = [sn for sn in all_sns if sn['survey'] in surveys_of_interest]

    #Simulation_function takes array 1 and array 2 and randomly associates them.
    # We want to end up with a new array of sn.
    # So we can associate to each sn a random residual and then reassign
    if simulation_method == 'random_sort':
        simulation_function = lambda sns, resids, resid_errs: randomSortData(sns, resids, y_errs = resid_errs)
    elif simulation_method == 'normal_dist':
        simulation_function = lambda sns, resids, resid_errs: simulateNormallyDistributedData(sns, lambda sn: sum(resids) / float(len(resids)), resid_errs)
    else:
        simulation_function = lambda sns, resids, resid_errs: randomSortData(sns, resids, y_errs = resid_errs)
    simulated_sns, randomized_resids, randomized_resid_errs = simulation_function(sn_of_interest, [sn['muDiff'] for sn in sn_of_interest], [sn['muErr'] for sn in sn_of_interest])
    for i in range(len(simulated_sns)):
        simulated_sns[i]['muDiff'] = randomized_resids[i]
        simulated_sns[i]['muErr'] = randomized_resid_errs[i]

    return [computeHemisphericResidual(simulated_sns, axis, opening_angle = opening_angle, statistic_funct = statistic_funct, statistic_err_funct = statistic_err_funct) for axis in axes]

def computeHemisphericResidual(all_sns, axis, opening_angle = math.pi / 2.0, statistic_funct = 'wmean', statistic_err_funct = 'wmean'):
    print 'Computing residual for axis ' + str(axis) 
    astro_archive = AstronomicalParameterArchive()
    deg_to_rad = astro_archive.getDegToRad()

    #make sure axis is a unit vector
    axis_mag = math.sqrt(sum([float(axis_elem) ** 2.0 for axis_elem in axis ]))
    axis = [float(axis_elem) / axis_mag for axis_elem in axis]
    
    upper_hem_sn = []
    lower_hem_sn = []

    for sn in all_sns:
        RA = sn['RA']
        RA = RA * deg_to_rad
        Dec = sn['Dec']
        Dec = Dec * deg_to_rad 
        sn_axis = [math.cos(RA) * math.sin(-1.0 * Dec + math.pi / 2.0), math.sin(RA) * math.sin(-1.0 * Dec + math.pi / 2.0), math.cos(-1.0 * Dec + math.pi / 2.0)]
        #angle between axis and sn can be computed from dot product of two axes
        if sum([ sn_axis[i] * axis[i] for i in range(len(axis)) ]) > 1.0:
            print 'dot_product = ' + str(sum([ sn_axis[i] * axis[i] for i in range(len(axis)) ])) + ' so one vector is not unit. '
            print 'sn_axis = ' + str(sn_axis)
            print 'axis = ' + str(axis) 
        sn_off_angle = math.acos(sum([ sn_axis[i] * axis[i] for i in range(len(axis)) ]))
        #print 'sn_off_angle = ' + str(sn_off_angle) 
        if sn_off_angle <= opening_angle: upper_hem_sn = upper_hem_sn + [sn]
        if sn_off_angle >= (math.pi - opening_angle): lower_hem_sn = lower_hem_sn + [sn]

    if statistic_funct is 'wmean':
        statistic_funct = lambda resids, resid_errs: sum([resids[i] * 1.0 / (resid_errs[i] ** 2.0)  for i in range(len(resids))]) / sum([1.0 / (resid_err ** 2.0) for resid_err in resid_errs]) 
    elif statistic_funct is 'mean':
        statistic_funct = lambda resids, resid_errs: sum(resids) / float(len(resids)) 
    if statistic_err_funct is 'none':
        statistic_err_funct = lambda resids, resid_errs: 1.0
    elif statistic_err_funct is 'wmean':
        statistic_err_funct = lambda resids, resid_errs: math.sqrt(1.0 / sum([1.0 / (resid_err ** 2.0) for resid_err in resid_errs]))
    #print "[sn['muDiff'] for sn in upper_hem_sn] = "
    #print [sn['muDiff'] for sn in upper_hem_sn]
    #print "[sn['muErr'] for sn in upper_hem_sn] = "
    #print [sn['muErr'] for sn in upper_hem_sn]
    #print "[sn['muDiff'] for sn in lower_hem_sn] = "
    #print [sn['muDiff'] for sn in lower_hem_sn]
    #print "[sn['muErr'] for sn in lower_hem_sn] = "
    #print [sn['muErr'] for sn in lower_hem_sn]
    if len(upper_hem_sn) == 0:
        print 'No sn observed in your upper hemisphere.  Cannot continue. '
        return [[0,0], [0,0]]
    if len(lower_hem_sn) == 0:
        print 'No sn observed in your lower hemisphere.  Cannot continue. '
        return [[0,0], [0,0]]
    upper_hem_val = statistic_funct([sn['muDiff'] for sn in upper_hem_sn], [sn['muErr'] for sn in upper_hem_sn])
    upper_hem_err = statistic_err_funct([sn['muDiff'] for sn in upper_hem_sn], [sn['muErr'] for sn in upper_hem_sn])
    lower_hem_val = statistic_funct([sn['muDiff'] for sn in lower_hem_sn], [sn['muErr'] for sn in lower_hem_sn])
    lower_hem_err = statistic_err_funct([sn['muDiff'] for sn in lower_hem_sn], [sn['muErr'] for sn in lower_hem_sn])

    return [[upper_hem_val, upper_hem_err], [lower_hem_val, lower_hem_err]]



    

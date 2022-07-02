from SNDataArchive import SNDataArchive
from RawSNDataStorer import RawSNDataStorer
from computeMuForCosmology import computeMuForCosmology
import numpy as np
#from query3dDustMap import query
import csv
from CosmologicalCalculator import CosmologicalCalculator
from astropy.coordinates import SkyCoord
import astropy.units as units
import cantrips as can
#from dustmaps.bayestar import BayestarWebQuery

def loadCepheidSN(cepheid_data_file = 'calibratorset.txt', cepheid_data_dir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/SNIsotropyProject/OriginalSNDataFiles/'):
    cepheid_data = can.readInColumnsToList(cepheid_data_file, file_dir = cepheid_data_dir, n_ignore = 1, delimiter = ' ', verbose = 0,)
    sn_ids = cepheid_data[0]
    hostmasses = cepheid_data[1]
    uMg_color = cepheid_data[4]
    return sn_ids


def loadSN(data_set, surveys_of_interest = ['all'], data_type = 'real', pull_extinctions = 0, save_cols = [], file_name = 'Sn_data_file.csv', relevent_indeces = 'all', zHD = 1,
           OmM = 0.3, OmL = 0.7, OmR = 0.0, H0 = 70.0, Om0 = 1.0, verbose = 0):

    sn_archive = SNDataArchive()
    sn_storer = RawSNDataStorer(data_set, data_type = data_type, relevent_indeces = relevent_indeces, verbose = verbose)
    sn_ras = sn_storer.getDataArray('RA')
    sn_ras = [float(ra) for ra in sn_ras]
    try:
        sn_decs = sn_storer.getDataArray('DECL')
    except KeyError:
        sn_decs = sn_storer.getDataArray('DEC')
    sn_decs = [float(dec) for dec in sn_decs]
    sn_zHDs = sn_storer.getDataArray('zHD')
    sn_zHDs = [float(z) for z in sn_zHDs]
    sn_zCMBs = sn_storer.getDataArray('zCMB')
    sn_zCMBs = [float(z) for z in sn_zCMBs]
    #print ('[sn_zHDs, sn_zCMBs] = ' + str([sn_zHDs, sn_zCMBs]))
    if zHD:
        sn_zs = sn_zHDs[:]
    else:
        sn_zs = sn_zCMBs[:]
    #print 'Read in zs. '

    sn_mBs = sn_storer.getDataArray('mB')
    sn_mBs = [float(mB) for mB in sn_mBs]
    sn_mBErrs = sn_storer.getDataArray('mBERR')
    sn_mBErrs = [float(mBErr) for mBErr in sn_mBErrs]
    #print 'Read in mBs. '
    #sn_mus = sn_storer.getDataArray('MU')
    sn_mus = sn_storer.getDataArray('mB')
    sn_mus = [float(mu) - sn_storer.getBigM()  for mu in sn_mus]
    #sn_muErrs = sn_storer.getDataArray('MUERR')
    sn_muErrs = sn_storer.getDataArray('mBERR')
    sn_muErrs = [float(muErr) for muErr in sn_muErrs]
    #print 'Read in mus. '
    sn_dls = [10.0 ** ((mu - 25.0) / 5.0) for mu in sn_mus]
    sn_dlErrs = [sn_muErrs[i] * 1 / 5.0 * np.log(10) * 10.0 ** ((sn_mus[i] - 25.0) / 5.0) for i in range(len(sn_mus)) ]
    #sn_mus_pred =  sn_storer.getDataArray('MUMODEL')
    #sn_mus_pred = [float(mu) for mu in sn_mus_pred]
    sn_surveys = sn_storer.getSurveys()

    sn_ids = (sn_storer.getDataArray('CID') if sn_storer.checkIfKeywordInDict('CID') else sn_storer.getDataArray('ID'))
    if verbose: print ('Read in all values. ')

    cosmo_calc = CosmologicalCalculator(zs = sn_zs, verbose = verbose )
    sn_ts = cosmo_calc.ts
    sn_taus = cosmo_calc.taus


    survey_color_map = sn_archive.getSurveyColorMap()
    sn_colors = np.array([ survey_color_map[survey] for survey in sn_surveys ])

    unique_surveys = np.unique(np.array(sn_surveys))
    if surveys_of_interest == ['all']:
        surveys_of_interest = unique_surveys

    #function_to_minimize = lambda OmM: np.sum(np.array(sn_mus) - np.array(sn_mus_pred))
    if verbose: print ('[OmM, OmL, OmR, H0, Om0] = ' + str([OmM, OmL, OmR, H0, Om0]) )
    expected_mus = [computeMuForCosmology(z, OmM = OmM, OmL = OmL, OmR = OmR, H0 = H0, Om0 = Om0, use_canon_params = 0) for z in sn_zs]
    mu_diffs = [sn_mus[i] - expected_mus[i] for i in range(len(sn_mBs))]
    sn_extinctions = []

    #sn_given_extinctions = sn_storer.getDataArray('MWEBV')

    all_surveys = np.unique(sn_surveys)
    sn_weighted_means_by_survey = {}
    for survey in all_surveys:
        muDiffs_by_survey = [mu_diffs[i] for i in range(len(sn_mus)) if sn_surveys[i] == survey]
        muErrs_by_survey = [sn_muErrs[i] for i in range(len(sn_mus)) if sn_surveys[i] == survey]
        weights_by_survey = [ 1.0 / muErrs_by_survey[i] ** 2.0 for i in range(len(muErrs_by_survey)) ]
        weighted_mean_for_survey = sum( [ muDiffs_by_survey[i] * weights_by_survey[i] for i in range(len(muDiffs_by_survey)) ] ) / sum(weights_by_survey)
        sn_weighted_means_by_survey[survey] = weighted_mean_for_survey


    #According to SCHLEGEL et al. 1998 (http://iopscience.iop.org/article/10.1086/305772/pdf),
    # the dust map is 'accurate to within 16%'.
    #So to determine the error, I currently just multiply the measured extinction by that constant.
    ext_err_const = 0.16
    if pull_extinctions:
        #sn_extinctions = query(sn_ras, sn_decs, coordsys = 'equ', mode = 'sfd')[u'EBV_SFD']
        #sn_ext_errs = [ext * ext_err_const for ext in sn_extinctions]
        coords = SkyCoord(sn_ras*units.deg, sn_decs*units.deg, distance=np.array([1.0 for sn in sn_decs])*units.Mpc, frame='icrs')
        bayestar = BayestarWebQuery(version='bayestar2017')
        sn_extinctions = bayestar(coords, mode='median')
        sn_ext_errs = [1.0 for i in range(len(sn_ras))]
    else:
        if verbose: print ('Extinctions not pulled; set to placeholder value.  DO NOT USE THIS ARRAY FOR ANY EXTINCTION CALCULATIONS!' )
        sn_extinctions = [0.0 for i in range(len(sn_ras))]
        sn_ext_errs = [1.0 for i in range(len(sn_ras))]



    sns = [ {'RA':sn_ras[i], 'Dec':sn_decs[i], 'z':sn_zs[i], 'zHD':sn_zHDs[i], 'zCMB':sn_zCMBs[i], 't':sn_ts[i], 'tau':sn_taus[i], 'mu':sn_mus[i], 'muErr':sn_muErrs[i],
             #'muDiff':mu_diffs[i], 'muPred':sn_mus_pred[i],
             'dl':sn_dls[i], 'dlErr':sn_dlErrs[i], 'survey':sn_surveys[i],
             'color':sn_colors[i], 'extinction':sn_extinctions[i], 'extinctionErr':sn_ext_errs[i], 'muDiffWMean':sn_weighted_means_by_survey[sn_surveys[i]],
             #'given_extinction':sn_given_extinctions[i],
              'SNID':sn_ids[i] }
            for i in range(len(sn_zs)) if sn_surveys[i] in surveys_of_interest]

    if len(save_cols) > 0:
        with open (file_name, 'wb') as csvfile:
            csv_writer = csv.writer(csvfile, delimiter = ',')
            csv_writer.writerow([column for column in save_cols])
            for sn in sns:
                csv_writer.writerow([sn[column] for column in save_cols])
    return sns

from SNDataArchive import SNDataArchive
from RawSNDataStorer import RawSNDataStorer
from computeMuForCosmology import computeMuForCosmology
import numpy as np
from query3dDustMap import query
import csv 

def loadArtificialSN(data_set, surveys_of_interest = ['all'], pull_extinctions = 1, save_cols = [], file_name = 'Sn_data_file.csv'):
    sn_archive = SNDataArchive()
    sn_storer = RawSNDataStorer(data_set, data_type = 'artificial')
    sn_ras = sn_storer.getDataArray('RA')
    sn_ras = [float(ra) for ra in sn_ras]
    sn_decs = sn_storer.getDataArray('DECL')
    sn_decs = [float(dec) for dec in sn_decs]
    sn_zs = sn_storer.getDataArray('zHD')
    sn_zs = [float(z) for z in sn_zs]
    sn_mBs = sn_storer.getDataArray('mB')
    sn_mBs = [float(mB) for mB in sn_mBs]
    sn_mBErrs = sn_storer.getDataArray('mBERR')
    sn_mBErrs = [float(mBErr) for mBErr in sn_mBErrs]
    sn_mus = sn_storer.getDataArray('MU')
    sn_mus = [float(mu) for mu in sn_mus]
    sn_muErrs = sn_storer.getDataArray('MUERR') 
    sn_muErrs = [float(muErr) for muErr in sn_muErrs]
    sn_dls = [10.0 ** ((mu - 25.0) / 5.0) for mu in sn_mus]
    sn_dlErrs = [sn_muErrs[i] * 1 / 5.0 * np.log(10) * 10.0 ** ((sn_mus[i] - 25.0) / 5.0) for i in range(len(sn_mus)) ]
    sn_mus_pred =  sn_storer.getDataArray('MUMODEL')
    sn_mus_pred = [float(mu) for mu in sn_mus_pred]
    sn_surveys = sn_storer.getSurveys()

    survey_color_map = sn_archive.getSurveyColorMap() 
    sn_colors = np.array([ survey_color_map[survey] for survey in sn_surveys ])

    unique_surveys = np.unique(np.array(sn_surveys))
    if surveys_of_interest == ['all']:
        surveys_of_interest = unique_surveys 

    expected_mus = [computeMuForCosmology(z) for z in sn_zs]
    mu_diffs = [sn_mus[i] - expected_mus[i] for i in range(len(sn_mBs))]
    sn_extinctions = sn_ras

    if pull_extinctions:
        sn_extinctions = query(sn_ras, sn_decs, coordsys = 'equ', mode = 'sfd')[u'EBV_SFD']
    else:
        print 'Extinctions not pulled; set to placeholder value.  DO NOT USE THIS ARRAY FOR ANY EXTINCTION CALCULATIONS!' 
        sn_extinctions = [0.0 for i in range(len(sn_ras))]



    sns = [ {'RA':sn_ras[i], 'Dec':sn_decs[i], 'z':sn_zs[i], 'mu':sn_mus[i], 'muDiff':mu_diffs[i], 
             'muErr':sn_muErrs[i], 'muPred':sn_mus_pred[i], 'dl':sn_dls[i], 'dlErr':sn_dlErrs[i], 'survey':sn_surveys[i], 'color':sn_colors[i], 'extinction':sn_extinctions[i]}
            for i in range(len(sn_zs)) if sn_surveys[i] in surveys_of_interest]

    if len(save_cols) > 0:
        with open (file_name, 'wb') as csvfile:
            csv_writer = csv.writer(csvfile, delimiter = ',')
            csv_writer.writerow([column for column in save_cols])
            for sn in sns:
                csv_writer.writerow([sn[column] for column in save_cols])

        

    return sns



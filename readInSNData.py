import math
import csv
import numpy as np
from SNDataArchive import SNDataArchive


def readInSNData(fileNum):
    data_archive = SNDataArchive()
    sn_file = data_archive.getSNDataFile(fileNum)

    orig_var_indeces = {'ID':0,            'CIDint':1,      'IDSURVEY':2,      'TYPE':3,              'FIELD':4,
                        'CUTFLAG_SNANA':5, 'zCMB':6,        'zCMBERR':7,       'zHD':8,               'zHDERR':9,
                        'VPEC':10,         'VPEC_ERR':11,   'HOST_LOGMASS':12, 'HOST_LOGMASS_ERR':13, 'SNRMAX1':14,
                        'SNRMAX2':15,      'SNRMAX3':16,    'PKMJD':17,        'PKMJDERR':18,         'x1':19,
                        'x1ERR':20,        'c':21,          'cERR':22,         'mB':23,               'mBERR':24,
                        'x0':25,           'x0ERR':26,      'COV_x1_c':27,     'COV_x1_x0':28,        'COV_c_x0':29,
                        'NDOF':30,         'FITCHI2':31,    'FITPROB':32,      'RA':33,               'DECL':34,
                        'TGAPMAX':35,      'MU':36,         'MUMODEL':37,      'MUERR':38,            'MUERR_RAW':39,
                        'MURES':40,        'MUPULL':41,     'ERRCODE':42,      'biasCor_mu':43,       'biasCorErr_mu':44,
                        'biasCor_mB':45,   'biasCor_x1':46, 'biasCor_c':47,    'biasScale_muCOV':48,  'IDSAMPLE':49       }

    file = open(sn_file, 'rb')
    row_num = 1
    start_row = 42
    start_column = 1
    parameter_arrays = [[] for i in range(len(orig_var_indeces))]
    print ('Reading in data from file: ' + sn_file) 
    with open(sn_file) as csvfile:
        myreader = csv.reader(csvfile, delimiter = ' ')
        for row in myreader:
            if row_num >= start_row:
                row = row[start_column:]
                row = [row_elem for row_elem in row if row_elem != '']
                parameter_arrays = [parameter_arrays[i] + [row[i]] for i in range(len(row))]
            row_num = row_num + 1


    x1 = parameter_arrays[orig_var_indeces['x1']]
    print x1
    return parameter_arrays

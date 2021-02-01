from SNDataArchive import SNDataArchive
import math
import csv
import numpy as np

class RawSNDataStorer:

    def getDataArray(self,key):
        if key in self.orig_var_indeces.keys():
            return self.parameter_arrays [self.orig_var_indeces[key]]
        else:
            print ('key: ' + str(key) + ' not found in loaded SN data file.  Providing placeholder (ideally, do not use for anything). ')
            return [self.catch_vals[key] for i in range(self.param_array_length) ]
            

    def getSurveyNumberMap(self):
        return self.survey_map

    def getSurveys(self,survey_numbers = []):
        if survey_numbers ==[]: survey_numbers = self.getDataArray('IDSURVEY') 
        return [self.survey_map[int(survey_number)] for survey_number in survey_numbers]



    def __init__(self, fileNum, data_type = 'real', relevent_indeces = 'all' ):
        data_archive = SNDataArchive() 
        sn_file = data_archive.getSNDataFile(fileNum, data_type = data_type)
        self.data_type = data_type 
        if data_type.lower() in ['artificial_lsst', 'art_lsst']:
            self.orig_var_indeces = {'CID':0,              'IDSURVEY':1,       'TYPE': 2,              'FIELD':3,            'CUTFLAG_SNANA': 4, 
                                     'zCMB':5,             'zCMBERR':6,        'zHD': 7,               'zHDERR':8,           'VPEC': 9,
                                     'VPEC_ERR':10,        'HOST_LOGMASS':11,  'HOST_LOGMASS_ERR': 12, 'SNRMAX1':13,         'SNRMAX2': 14,
                                     'SNRMAX3':15,         'PKMJD':16,         'PKMJDERR':17,          'x1':18,              'x1ERR': 19,
                                     'c':20,               'cERR':21,          'mB':22,                'mBERR':23,           'x0': 24,
                                     'x0ERR':25,           'COV_x1_c':26,      'COV_x1_x0':27,         'COV_c_x0':28,        'NDOF': 29,
                                     'FITCHI2':30,         'FITPROB':31,       'SIM_TYPE_INDEX':32,    'SIM_NONIA_INDEX':33, 'SIM_LIBID': 34,
                                     'SIM_NGEN_LIBID':35,  'SIM_ZCMB':36,      'SIM_DLMAG':37,         'SIM_PKMJD':38,       'SIM_x1':39,
                                     'SIM_c':40,           'SIM_alpha':41,     'SIM_beta':42,          'SIM_x0':43,          'SIM_mB': 44,
                                     'SIM_AV':45,          'MU':46,            'MUMODEL':47,           'MUERR':48,           'MUERR_RAW': 49,
                                     'MURES':50,           'MUPULL':51,        'ERRCODE':52,           'SIM_MU':53,          'SIM_MB': 54,
                                     'biasCor_mu':55,      'biasCorErr_mu':56, 'biasCor_mB':57,        'biasCor_x1':58,      'biasCor_c': 59,
                                     'biasScale_muCOV':60, 'IDSAMPLE': 61}
            self.catch_vals = {'RA':0.0, 'DECL':0.0}
        elif data_type.lower() in ['art_survey','artifical_survey','art_surveys','artififical_survey']: 
            self.orig_var_indeces = {'CID':0,              'IDSURVEY':1,       'TYPE':2,              'FIELD':3,            'CUTFLAG_SNANA':4,
                                     'zCMB':5,             'zCMBERR':6,        'zHD':7,               'zHDERR':8,           'VPEC':9,
                                     'VPEC_ERR':10,        'HOST_LOGMASS':11,  'HOST_LOGMASS_ERR':12, 'SNRMAX1':13,         'SNRMAX2':14,
                                     'SNRMAX3':15,         'PKMJD':16,         'PKMJDERR':17,         'x1':18,              'x1ERR':19,
                                     'c':20,               'cERR':21,          'mB':22,               'mBERR':23,           'x0':24,
                                     'x0ERR':25,           'COV_x1_c':26,      'COV_x1_x0':27,        'COV_c_x0':28,        'NDOF':29,
                                     'FITCHI2':30,         'FITPROB':31,       'SIM_TYPE_INDEX':32,   'SIM_NONIA_INDEX':33, 'SIM_LIBID':34,
                                     'SIM_NGEN_LIBID':35,  'SIM_ZCMB':36,      'SIM_DLMAG':37,        'SIM_PKMJD':38,       'SIM_x1':39,
                                     'SIM_c':40,           'SIM_alpha':41,     'SIM_beta':42,         'SIM_x0':43,          'SIM_mB':44,
                                     'SIM_AV':45,          'MU':46,            'MUMODEL':47,          'MUERR':48,           'MUERR_RAW':49,
                                     'MURES':50,           'MUPULL':51,        'ERRCODE':52,          'SIM_MU':53,          'SIM_MB':54,
                                     'biasCor_mu':55,      'biasCorErr_mu':56, 'biasCor_mB':57,       'biasCor_x1':58,      'biasCor_c':59,
                                     'biasScale_muCOV':60, 'IDSAMPLE':61 }
            self.catch_vals = {'RA':0.0, 'DECL':0.0}
        #Artificial function, but we have simply replaced the measured mus with some mus of our artificial simulation 
        elif data_type.lower() in ['artificial_funct', 'art_funct']:
            self.orig_var_indeces = {'ID':0,            'CIDint':1,      'IDSURVEY':2,      'TYPE':3,              'FIELD':4,
                                     'CUTFLAG_SNANA':5, 'zCMB':6,        'zCMBERR':7,       'zHD':8,               'zHDERR':9,
                                     'VPEC':10,         'VPEC_ERR':11,   'HOST_LOGMASS':12, 'HOST_LOGMASS_ERR':13, 'SNRMAX1':14,
                                     'SNRMAX2':15,      'SNRMAX3':16,    'PKMJD':17,        'PKMJDERR':18,         'x1':19,
                                     'x1ERR':20,        'c':21,          'cERR':22,         'mB':23,               'mBERR':24,
                                     'x0':25,           'x0ERR':26,      'COV_x1_c':27,     'COV_x1_x0':28,        'COV_c_x0':29,
                                     'NDOF':30,         'FITCHI2':31,    'FITPROB':32,      'RA':33,               'DECL':34,
                                     'TGAPMAX':35,      'MU':36,         'MUMODEL':37,      'MUERR':38,            'MUERR_RAW':39,
                                     'MURES':40,        'MUPULL':41,     'ERRCODE':42,      'biasCor_mu':43,       'biasCorErr_mu':44,
                                     'biasCor_mB':45,   'biasCor_x1':46, 'biasCor_c':47,    'biasScale_muCOV':48,  'IDSAMPLE':49       }
            self.catch_vals = {}
        elif data_type.lower() in ['old_real', 'real_old']: #need to account for old SN files from Dan Scolnic that have slightly fewer fields 
            self.orig_var_indeces = {'ID':0,            'CIDint':1,      'IDSURVEY':2,      'TYPE':3,              'FIELD':4,
                                     'CUTFLAG_SNANA':5, 'zCMB':6,        'zCMBERR':7,       'zHD':8,               'zHDERR':9,
                                     'VPEC':10,         'VPEC_ERR':11,   'HOST_LOGMASS':12, 'HOST_LOGMASS_ERR':13, 'SNRMAX1':14,
                                     'SNRMAX2':15,      'SNRMAX3':16,    'PKMJD':17,        'PKMJDERR':18,         'x1':19,
                                     'x1ERR':20,        'c':21,          'cERR':22,         'mB':23,               'mBERR':24,
                                     'x0':25,           'x0ERR':26,      'COV_x1_c':27,     'COV_x1_x0':28,        'COV_c_x0':29,
                                     'NDOF':30,         'FITCHI2':31,    'FITPROB':32,      'RA':33,               'DECL':34,
                                     'TGAPMAX':35,      'MU':36,         'MUMODEL':37,      'MUERR':38,            'MUERR_RAW':39,
                                     'MURES':40,        'MUPULL':41,     'ERRCODE':42,      'biasCor_mu':43,       'biasCorErr_mu':44,
                                     'biasCor_mB':45,   'biasCor_x1':46, 'biasCor_c':47,    'biasScale_muCOV':48,  'IDSAMPLE':49       }
            self.catch_vals = {}
        else: 
            self.orig_var_indeces = {'ID':0,            'CIDint':1,           'IDSURVEY':2,       'TYPE':3,              'FIELD':4,
                                     'CUTFLAG_SNANA':5, 'zCMB':6,             'zCMBERR':7,        'zHD':8,               'zHDERR':9,
                                     'VPEC':10,         'VPEC_ERR':11,        'HOST_LOGMASS':12,  'HOST_LOGMASS_ERR':13, 'SNRMAX1':14,
                                     'SNRMAX2':15,      'SNRMAX3':16,         'PKMJD':17,         'PKMJDERR':18,         'x1':19,
                                     'x1ERR':20,        'c':21,               'cERR':22,          'mB':23,               'mBERR':24,
                                     'x0':25,           'x0ERR':26,           'COV_x1_c':27,      'COV_x1_x0':28,        'COV_c_x0':29,
                                     'NDOF':30,         'FITCHI2':31,         'FITPROB':32,       'RA':33,               'DECL':34,
                                     'TGAPMAX':35,      'TrestMIN':36,        'TrestMAX':37,      'MWEBV':38,            'MU':39,
                                     'MUMODEL':40,      'MUERR':41,           'MUERR_RAW':42,     'MURES':43,            'MUPULL':44,
                                     'ERRCODE':45,      'biasCor_mu':46,      'biasCorErr_mu':47, 'biasCor_mB':48,       'biasCor_x1':49,
                                     'biasCor_c':50,    'biasScale_muCOV':51, 'IDSAMPLE':52     }
            self.catch_vals = {}

        file = open(sn_file, 'rb')
        row_num = 1
        self.to_be_read_length = -1
        if data_type.lower() in ['artificial_lsst', 'art_lsst']:
            self.start_row = 46
            self.to_be_read_length = 400272
            self.important_indeces = np.linspace(0, self.start_row)
        elif data_type.lower() in ['art_survey','artifical_survey','art_surveys','artififical_survey']: 
            self.start_row = 65
        elif data_type.lower() in ['artificial_funct', 'art_funct']:
            self.start_row = 42
        elif data_type.lower() in ['old_real', 'real_old']:
            self.start_row = 42
        else:
            self.start_row = 70
        if relevent_indeces is 'all': relevent_indeces = range( self.start_row )
            
        self.start_column = 1
        if self.to_be_read_length > 0:
            self.parameter_arrays = np.empty(( self.to_be_read_length, len(self.orig_var_indeces) ), dtype = "S10")
        else: 
            self.parameter_arrays = [[] for i in range(len(self.orig_var_indeces))]
        #with open(sn_file) as csvfile:
        #    myreader = csv.reader(csvfile, delimiter = ' ')
        #    for row in myreader:
        #        if row_num >= self.start_row:
        #            row = row[self.start_column:]
        #            row = [row_elem for row_elem in row if row_elem != '']
        #            if row_num % 100 ==1: 
        #                print 'row_num = ' + str(row_num) 
        #        row_num = row_num + 1
        #print 'Done reading in with row_number = ' + str(row_num) 
        with open(sn_file) as csvfile:
            myreader = csv.reader(csvfile, delimiter = ' ')
            for row in myreader:
                if row_num >= self.start_row:
                    row = row[self.start_column:]
                    row = [row_elem for row_elem in row if row_elem != '']
                    if self.to_be_read_length > 0:
                        self.parameter_arrays[row_num - self.start_row,:] = row
                    else: 
                        self.parameter_arrays = [self.parameter_arrays[i] + [row[i]] for i in range(len(self.parameter_arrays))]
                row_num = row_num + 1
        print ('total_read_in = ' + str(row_num - self.start_row))
        if self.to_be_read_length > 0:
            self.parameter_arrays = self.parameter_arrays.transpose().tolist()

        self.param_array_length = len(self.parameter_arrays[0])

        self.survey_map = {1:'SDSS', 3:'ESSENCE', 4:'SNLS', 5:'CSP', 6:'SUBARU',
                      10:'DES', 11:'VIDEO',12:'LSST',14:'XIAN', 15:'PS1MD', 
                      16:'JPAS', 50:'LOWZ', 51:' KAIT', 52:'SNF', 53:'CFA3',
                      54:'CFA4', 55:'SNUO2', 56:'SWIFT',57:'KAITM',58:'KAITW',
                      61:'CFA1', 62:'CFA2', 63:'CFA3S', 64:'CFA3K', 65:'CFA4p1', 
                      66:'CFA4p2', 67:'CFA4p3', 68:'CSP3', 70:'PS1_LOWZ_COMBINED', 71:' PS1_HST_COMBINED',
                      100:'HST', 101:'SNAP', 103:'WFIRST', 104:'JWST', 105:'EUCLID',
                      106:'CANDELS', 107:'CLASH', 110:'GAIA', 111:'HDF', 150:'FOUNDATION',
                      180:'SIMSED', 181:'FLASH', 190:'TEST'      }

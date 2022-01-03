from SNDataArchive import SNDataArchive
import math
import csv
import numpy as np

class RawSNDataStorer:

    def checkIfKeywordInDict(self, key):
        return key in self.orig_var_indeces.keys()

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

    def __init__(self, fileNum, data_type = 'real', relevent_indeces = 'all', verbose = 1 ):
        self.verbose = verbose
        data_archive = SNDataArchive()
        sn_file = data_archive.getSNDataFile(fileNum, data_type = data_type)
        #print ('sn_file = ' + str(sn_file))
        if verbose: self.data_type = data_type
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
        elif data_type.lower() in ['art_pantheon','artifical_pantheon']:
            self.orig_var_indeces = {'CID':0,                 'CIDint':1,            'IDSURVEY':2,        'TYPE':3,           'FIELD':4,
                                     'CUTFLAG_SNANA':5,       'zHEL':6,              'zHELERR':7,         'zCMB':8,           'zCMBERR':9,
                                     'zHD':10,                'zHDERR':11,           'VPEC':12,           'VPECERR':13,       'MWEBV':14,
                                     'HOST_LOGMASS':15,       'HOST_LOGMASS_ERR':16, 'HOST_sSFR':17,      'HOST_sSFR_ERR':18, 'PKMJDINI':19,
                                     'SNRMAX1':20,            'SNRMAX2':21,          'SNRMAX3':22,        'PKMJD':23,         'PKMJDERR':24,
                                     'x1':25,                 'x1ERR':26,            'c':27,              'cERR':28,          'mB':29,
                                     'mBERR':30,              'x0':31,               'x0ERR':32,          'COV_x1_c':33,      'COV_x1_x0':34,
                                     'COV_c_x0':35,           'NDOF':36,             'FITCHI2':37,        'FITPROB':38,       'SIM_TYPE_INDEX':39,
                                     'SIM_TEMPLATE_INDEX':40, 'SIM_LIBID':41,        'SIM_NGEN_LIBID':42, 'SIM_ZCMB':43,      'SIM_ZFLAG':44,
                                     'SIM_VPEC':45,           'SIM_DLMAG':46,        'SIM_PKMJD':47,      'SIM_x1':48,        'SIM_c':49,
                                     'SIM_alpha':50,          'SIM_beta':51,         'SIM_x0':52,         'SIM_mB':53,        'SIM_AV':54,
                                     'SIM_RV':55,             'RA':56,               'DECL':57,            'HOST_RA':58,       'HOST_DEC':59,
                                     'HOST_ANGSEP':60,        'TGAPMAX':61,          'TrestMIN':62,       'TrestMAX':63,      'CUTMASK':64,
                                     'MU':65,                 'MUMODEL':66,          'MUERR':67,          'MUERR_RENORM':68,  'MUERR_RAW':69,
                                     'MUERR_VPEC':70,         'MURES':71,            'MUPULL':72,         'M0DIF':73,         'M0DIFERR':74,
                                     'CHI2':75,               'biasCor_mu':76,       'biasCorErr_mu':77,  'biasCor_mB':78,    'biasCor_x1':79,
                                     'biasCor_c':80,          'biasScale_muCOV':81,  'IDSAMPLE':82,       'IZBIN':83}
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
        elif data_type.lower() in ['realplus', 'real_plus', 'realp', 'real_p', 'pantheon_plus']:
            self.orig_var_indeces = {'CID':0,              'CIDint':1,         'IDSURVEY':2,            'TYPE':3,              'FIELD':4,
                                     'CUTFLAG_SNANA':5,    'ERRFLAG_FIT':6,    'zHEL':7,                'zHELERR':8,           'zCMB':9,
                                     'zCMBERR':10,         'zHD':11,           'zHDERR':12,             'VPEC':13,             'VPECERR':14,
                                     'MWEBV':15,           'HOST_LOGMASS':16,  'HOST_LOGMASS_ERR':17,   'HOST_sSFR':18,        'HOST_sSFR_ERR':19,
                                     'PKMJDINI':20,        'SNRMAX1':21,       'SNRMAX2':22,            'SNRMAX3':23,          'PKMJD':24,
                                     'PKMJDERR':25,        'x1':26,            'x1ERR':27,              'c':28,                'cERR':29,
                                     'mB':30,              'mBERR':31,         'x0':32,                 'x0ERR':33,            'COV_x1_c':34,
                                     'COV_x1_x0':35,       'COV_c_x0':36,      'NDOF':37,               'FITCHI2':38,          'FITPROB':39,
                                     'RA':40,              'DEC':41,           'HOST_RA':42,            'HOST_DEC':43,         'HOST_ANGSEP':44,
                                     'TGAPMAX':45,         'TrestMIN':46,      'TrestMAX':47,           'ELU':48,              'HOSTGAL_SFR':49,
                                     'HOSTGAL_SFR_ERR':50, 'HOSTGAL_sSFR':51,  'HOSTGAL_sSFR_ERR':52,   'CUTMASK':53,          'MU':54,
                                     'MUMODEL':55,         'MUERR':56,         'MUERR_RENORM':57,       'MUERR_RAW':58,        'MUERR_VPEC':59,
                                     'MURES':60,           'MUPULL':61,        'M0DIF':62,              'M0DIFERR':63,         'CHI2':64,
                                     'biasCor_mu':65,      'biasCorErr_mu':66, 'biasCor_muCOVSCALE':67, 'biasCor_muCOVADD':68, 'IDSAMPLE':69,
                                     'IZBIN':70     }
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
        elif data_type.lower() in ['realplus', 'real_plus', 'realp', 'real_p', 'pantheon_plus']:
            self.start_row = 65
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
            print ('sn_file = ' + str(sn_file))
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
        if self.verbose: print ('total_read_in = ' + str(row_num - self.start_row))
        if self.to_be_read_length > 0:
            self.parameter_arrays = self.parameter_arrays.transpose().tolist()

        self.param_array_length = len(self.parameter_arrays[0])

        self.survey_map = {1:'SDSS', 3:'ESSENCE', 4:'SNLS', 5:'CSP', 6:'SUBARU',
                      10:'DES', 11:'VIDEO',12:'LSST',14:'XIAN', 15:'PS1MD',
                      16:'JPAS', 17:'ZTF_MSIP', 18:'ASASSN',
                      49:'JRK07', 50:'LOWZ', 51:'KAIT', 52:'SNF', 53:'CFA3',
                      54:'CFA4', 55:'SNUO2', 56:'SWIFT',57:'KAITM',58:'KAITW', 59:'SWIFTNEW',
                      61:'CFA1', 62:'CFA2', 63:'CFA3S', 64:'CFA3K', 65:'CFA4p1',
                      66:'CFA4p2', 67:'CFA4p3', 68:'CFA5', 69:'CfAIR2',
                      70:'PS1_LOWZ_COMBINED', 71:' PS1_HST_COMBINED', 72:'DES3YR_LOWZ_COMBINED', 73:'JLA_LOWZ_COMBINED',
                      100:'HST', 101:'SNAP', 102:'NGRST', 103:'WFIRST', 104:'JWST', 105:'EUCLID',
                      106:'CANDELS', 107:'CLASH', 110:'GAIA', 111:'HDF', 150:'FOUNDATION',
                      180:'SIMSED', 181:'FLASH', 190:'TEST', 191:'OTHER'   }

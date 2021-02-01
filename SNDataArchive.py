from AstronomicalParameterArchive import AstronomicalParameterArchive
from DirectoryArchive import DirectoryArchive

class SNDataArchive:

    def getSurveyColorMap(self):
        return self.survey_color_map

    def getSNDataFile(self, fileNum, data_type = 'real'):
        if data_type.lower() in ['art_lsst','artificial_lsst']:
            return self.artificial_lsst_sn_data_dir + self.artificial_lsst_sn_data_files[fileNum]
        elif data_type.lower() in ['artificial_funct', 'art_funct']:
            return self.artificial_funct_sn_data_dir + self.artificial_funct_sn_data_files[fileNum]
        elif data_type.lower in ['art_survey','artifical_survey','art_surveys','artififical_survey']:
            return self.artificial_surveys_sn_data_dir + self.artificial_surveys_sn_data_files[fileNum]
        elif data_type.lower() in ['old_real', 'real_old']:
            return self.sn_data_dir + self.old_sn_data_files[fileNum]
        else: 
            return self.sn_data_dir + self.sn_data_files[fileNum] 
   
    
    def __init__(self):
        dir_archive = DirectoryArchive()
        self.artificial_surveys_sn_data_files = {'0c':'SALT2mu_simoptc0_0.fitres',   '1c':'SALT2mu_simoptc0_1.fitres',   '2c':'SALT2mu_simoptc0_2.fitres',
                                                 '3c':'SALT2mu_simoptc0_3.fitres',   '4c':'SALT2mu_simoptc0_4.fitres',   '5c':'SALT2mu_simoptc0_5.fitres',
                                                 '6c':'SALT2mu_simoptc0_6.fitres',   '7c':'SALT2mu_simoptc0_7.fitres',   '8c':'SALT2mu_simoptc0_8.fitres',
                                                 '9c':'SALT2mu_simoptc0_9.fitres',   '10c':'SALT2mu_simoptc0_10.fitres', '11c':'SALT2mu_simoptc0_11.fitres',
                                                 '12c':'SALT2mu_simoptc0_12.fitres', '13c':'SALT2mu_simoptc0_13.fitres', '14c':'SALT2mu_simoptc0_14.fitres',
                                                 '15c':'SALT2mu_simoptc0_15.fitres', '16c':'SALT2mu_simoptc0_16.fitres', '17c':'SALT2mu_simoptc0_17.fitres',
                                                 '18c':'SALT2mu_simoptc0_18.fitres', '19c':'SALT2mu_simoptc0_19.fitres', '20c':'SALT2mu_simoptc0_20.fitres',
                                                 '21c':'SALT2mu_simoptc0_21.fitres', '22c':'SALT2mu_simoptc0_22.fitres', '23c':'SALT2mu_simoptc0_23.fitres',
                                                 '24c':'SALT2mu_simoptc0_24.fitres', '25c':'SALT2mu_simoptc0_25.fitres', '26c':'SALT2mu_simoptc0_26.fitres',
                                                 '27c':'SALT2mu_simoptc0_27.fitres',
                                                 '0g':'SALT2mu_simoptg0_0.fitres',   '1g':'SALT2mu_simoptg0_1.fitres',   '2g':'SALT2mu_simoptg0_2.fitres',
                                                 '3g':'SALT2mu_simoptg0_3.fitres',   '4g':'SALT2mu_simoptg0_4.fitres',   '5g':'SALT2mu_simoptg0_5.fitres',
                                                 '6g':'SALT2mu_simoptg0_6.fitres',   '7g':'SALT2mu_simoptg0_7.fitres',   '8g':'SALT2mu_simoptg0_8.fitres',
                                                 '9g':'SALT2mu_simoptg0_9.fitres',   '10g':'SALT2mu_simoptg0_10.fitres', '11g':'SALT2mu_simoptg0_11.fitres',
                                                 '12g':'SALT2mu_simoptg0_12.fitres', '13g':'SALT2mu_simoptg0_13.fitres', '14g':'SALT2mu_simoptg0_14.fitres',
                                                 '15g':'SALT2mu_simoptg0_15.fitres', '16g':'SALT2mu_simoptg0_16.fitres', '17g':'SALT2mu_simoptg0_17.fitres',
                                                 '18g':'SALT2mu_simoptg0_18.fitres', '19g':'SALT2mu_simoptg0_19.fitres', '20g':'SALT2mu_simoptg0_20.fitres',
                                                 '21g':'SALT2mu_simoptg0_21.fitres', '22g':'SALT2mu_simoptg0_22.fitres', '23g':'SALT2mu_simoptg0_23.fitres',
                                                 '24g':'SALT2mu_simoptg0_24.fitres', '25g':'SALT2mu_simoptg0_25.fitres', '26g':'SALT2mu_simoptg0_26.fitres',
                                                 '27g':'SALT2mu_simoptg0_27.fitres'  }
        
        self.artificial_lsst_sn_data_files = {1:'SALT2mu_fitopt_CT_LSST.fitres'} 
        self.artificial_funct_sn_data_files = {1: 'SALTmu_FITOPT000_MUOPT000_artBumpA2_10_mu2000_sig200_1.FITRES',
                                               2: 'SALTmu_FITOPT000_MUOPT000_artBumpA2_10_mu2000_sig200_2.FITRES',
                                               3: 'SALTmu_FITOPT000_MUOPT000_artBumpA4_10_mu2000_sig400_1.FITRES', 
                                               4: 'SALTmu_FITOPT000_MUOPT000_artBumpA4_10_mu2000_sig400_2.FITRES', 
                                               5: 'SALTmu_FITOPT000_MUOPT000_artZerod.FITRES'}
        self.old_sn_data_files = {1: 'SALT2mu_FITOPT000_MUOPT000.FITRES',
                                  2: 'SALT2mu_FITOPT011_MUOPT000.FITRES'}
        self.sn_data_files = {1: 'Ancillary_G10.FITRES'}
        self.artificial_surveys_sn_data_dir = dir_archive.getArtificialSurveySNDataDirectory()
        self.artificial_lsst_sn_data_dir = dir_archive.getArtificialLSSTSNDataDirectory()
        self.artificial_funct_sn_data_dir = dir_archive.getArtificialFunctSNDataDirectory() 
        self.sn_data_dir = dir_archive.getSNDataDirectory() 

        self.survey_color_map = {'SDSS':'yellowgreen',        'ESSENCE':'azure',        'SNLS':'purple',        'CSP':'blue',
                                 'SUBARU':'blueviolet',       'DES':'coral',            'VIDEO':'crimson',     'LSST':'cyan',
                                 'XIAN':'darkblue',           'PS1MD':'darkgoldenrod',  'JPAS':'darkgreen',    'LOWZ':'darkmagenta',
                                 'KAIT':'darkolivegreen',     'SNF':'darkorange',       'CFA3':'firebrick',    'CFA4':'fuchsia',
                                 'SNUO2':'brown',            'SWIFT':'khaki',          'KAITM':'greenyellow', 'KAITW':'lightgrey',
                                 'CFA1':'lightblue',          'CFA2':'lightpink',       'CFA3S':'green',       'CFA3K':'magenta',
                                 'CFA4p1':'indigo',           'CFA4p2':'orange',        'CFA4p3':'bisque',     'CSP3':'orangered',
                                 'PS1_HST_COMBINED':'silver', 'HST':'lime',             'SNAP':'violet',       'PS1_LOWZ_COMBINED':'red', 
                                 'WFIRST':'aliceblue',      'JWST':'skyblue',         'EUCLID':'tomato',     'CANDELS':'red',
                                 'CLASH':'turquoise',         'GAIA':'teal',            'HDF':'mintcream',     'FOUNDATION':'olive',
                                 'SIMSED':'midnightblue',     'FLASH':'mediumseagreen', 'TEST':'mediumblue'    }  

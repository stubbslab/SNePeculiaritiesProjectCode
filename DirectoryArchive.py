
class DirectoryArchive:

    def getRandomizedBestFitParamsDir(self):
        return self.random_best_fit_params_directory

    def getRandomizedDrawResultsDir(self):
        return self.random_draw_results_directory

    def getSNDataDirectory (self):
        return self.sn_data_directory

    def getPlotDirectory(self):
        return self.plot_directory

    def getArtificialLSSTSNDataDirectory(self):
        return self.artificial_lsst_sn_data_directory

    def getArtificialSurveySNDataDirectory(self):
        return self.artificial_survey_sn_data_directory

    def getArtificialFunctSNDataDirectory(self):
        return self.artificial_funct_sn_data_directory

    def getChiSqrMeasurementsDirectory(self):
        return self.chi_sqr_calculations_directory

    def getRandomizedBestFitParamsDir(self):
        return self.randomized_best_fit_params_dir

    def __init__(self):
        self.chris_SNIsotropy_directory = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/SNIsotropyProject/'
        self.plot_directory = self.chris_SNIsotropy_directory + 'plots/'
        self.sn_data_directory = self.chris_SNIsotropy_directory + 'OriginalSNDataFiles/'
        self.artificial_lsst_sn_data_directory = self.chris_SNIsotropy_directory + 'LSST_data_sets/'
        self.artificial_survey_sn_data_directory = self.chris_SNIsotropy_directory + 'other_art_survey_data_sets/'
        self.artificial_funct_sn_data_directory = self.chris_SNIsotropy_directory + 'data_from_artificial_functs/'
        self.random_draw_results_directory = self.chris_SNIsotropy_directory + 'randomDrawResults/'
        self.chi_sqr_calculations_directory = self.chris_SNIsotropy_directory + 'chiSqrMeasurements/'
        self.randomized_best_fit_params_dir = self.chris_SNIsotropy_directory + 'randomBestFitParams/'
        self.random_draw_results_directory = self.chris_SNIsotropy_directory + 'randomDrawResults/'
        self.artificial_sn_data_directory = self.chris_SNIsotropy_directory + 'randomDrawResults/'

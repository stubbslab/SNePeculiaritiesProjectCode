from loadSN import loadSN


def shuffleAndMeasureBestFitForSNResiduals (data_set, fit_function, surveys_to_fit = ['PS1MD']]):

    all_sn = loadSN(data_set)
    
    

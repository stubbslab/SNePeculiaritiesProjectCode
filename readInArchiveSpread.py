import numpy as np
import pickle 
from DirectoryArchive import DirectoryArchive 

def readInArchiveSpread(file_name = None):
    if file_name is None:
        file_name = 'standard_rand_fit_results_file.npy'
    dir_archive = DirectoryArchive()
    chi_sqr_dir = dir_archive.getChiSqrMeasurementsDirectory()
    chi_sqr_data_file = chi_sqr_dir + file_name
    print 'Loading file '+ chi_sqr_data_file
    fileObject = open(chi_sqr_data_file, 'r')
    return pickle.load(fileObject) 

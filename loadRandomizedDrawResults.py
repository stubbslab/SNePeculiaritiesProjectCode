import numpy as np
from DirectoryArchive import DirectoryArchive 

def loadRandomizedDrawResults(file_name):
    dir_archive = DirectoryArchive()
    load_dir = dir_archive.getRandomizedDrawResultsDir()
    return (np.load(load_dir + file_name).item()) 

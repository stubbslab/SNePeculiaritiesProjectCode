import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from logList import logList
from  RawSNDataStorer import RawSNDataStorer 

def showSNPositions(fileNum, RA_lims = [0.0, 360.0], Dec_lims = [-90.0, 90.0]):
    sn_storer = RawSNDataStorer(fileNum)

    sn_ras = sn_storer.getDataArray('RA')
    sn_decs = sn_storer.getDataArray('DECL') 
    
    plt.scatter(sn_ras, sn_decs)
    plt.xlim( RA_lims )
    plt.ylim( Dec_lims )
    plt.xlabel('RA')
    plt.ylabel('Dec')
    plt.title('SN Positions on Sky')
    plt.grid() 
    plt.show() 

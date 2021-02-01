import numpy as np
import math
from AstronomicalParameterArchive import AstronomicalParameterArchive 

#coordinate transformations pulled from stack exchange at this link: 
    # http://math.stackexchange.com/questions/92301/transforming-from-one-spherical-coordinate-system-to-another

#Simple angular separation equation pulled from 
    #https://en.wikipedia.org/wiki/Great-circle_distance
def calculateAngularSeparation (cent_ras, cent_decs, targ_RA, targ_Dec):
    astro_archive = AstronomicalParameterArchive()
    deg_to_rad = astro_archive.getDegToRad() #math.pi/180.0
    
    #corr_ras = targ_RA - cent_ras
    
    #x1=np.cos(np.array(targ_Dec) * deg_to_rad)*np.cos(corr_ras * deg_to_rad)
    #x2=np.cos(np.array(targ_Dec) * deg_to_rad)*np.sin(corr_ras * deg_to_rad)
    #x3=np.sin(np.array(targ_Dec) * deg_to_rad)

    #x1bar=np.cos(cent_decs * deg_to_rad)*x1+np.sin(cent_decs * deg_to_rad)*x3
    #x2bar=x2
    #x3bar=-np.sin(cent_decs * deg_to_rad)*x1 + np.cos(cent_decs * deg_to_rad)*x3
    #corr_ras = np.angle(x1bar + x2bar*1j, deg = True) 
    #print len(corr_ra)
    #corr_decs = np.arcsin(x3bar) * (1 / deg_to_rad)
    radial_dists = np.arccos(
                             np.sin(cent_decs * deg_to_rad) * np.sin(targ_Dec * deg_to_rad)
                             + np.cos(cent_decs * deg_to_rad) * np.cos(targ_Dec * deg_to_rad)
                              * np.cos((targ_RA - cent_ras) * deg_to_rad)                       )
    
    return radial_dists * 1.0 / deg_to_rad
    

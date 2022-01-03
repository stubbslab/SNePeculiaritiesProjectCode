import numpy as np
import matplotlib.pyplot as plt
import wmom as wmom
#import plotsetup
from matplotlib import gridspec
import matplotlib.ticker as ticker
import linef
#plotsetup.halfpaperfig()
import getsurveyname
#plt.style.use('custom')
import pandas as pd
import yaml
import gzip

survdict = {1:'SDSS',4:'SNLS',5:'CSP',10:'DES',15:'PS1MD',50:'LOWZ/JRK07',51:'KAIT',56:'SWIFT',58:'KAITW',59:'SWIFTNEW',61:'CFA1',62:'CFA2',63:'CFA3S',64:'CFA3K',65:'CFA4p2',66:'CFA4p3',100:'HST',101:'SNAP',106:'CANDELS',150:'FOUND'}
cset = list(pd.read_csv('calibratorset.txt',delim_whitespace=True)['SNID'].array)

def getnuisance(fitres):
 beta = np.nan
 betaerr = np.nan
 alpha = np.nan
 alphaerr = np.nan
 gamma = np.nan
 gammaerr = np.nan


 try:
  lines = open(fitres).readlines()[:200]
  doblines = False
 except:
  lines = gzip.open(fitres).readlines()[:200]
  doblines = True
 #lines = gzip.open(fitres).readlines()[:200]
 for bline in lines:
  if doblines:
   line = bline.decode("utf-8")
  else:
   line = bline
  try:
   if 'beta0          =' in line:
    beta = float(line.split()[3])
    betaerr = float(line.split()[5])
   if 'alpha0         =' in line:
    alpha = float(line.split()[3])
    alphaerr = float(line.split()[5])
   if 'gamma0         =' in line:
    gamma = float(line.split()[3])
    gammaerr = float(line.split()[5])
  except:
   pass
 return alpha,alphaerr,beta,betaerr,gamma,gammaerr


def h0(opt,optzero=None,zmin=.01,zmax=.15,name='',outdf=None,scale=1.):

 covids = pd.read_csv(name+'/full/indiv_covmat/data_wCID.txt',delim_whitespace=True)
 covids['UNIQ'] = covids['CID']+covids['IDSURVEY'].astype(str)

 if optzero is None:
  optzero = opt

 #read in nominal and sys
 df = pd.read_csv(opt,delim_whitespace=True,comment='#')
 dfzero = pd.read_csv(optzero,delim_whitespace=True,comment='#')

 df['UNIQ'] = df['CID']+df['IDSURVEY'].astype(str)
 dfzero['UNIQ'] = dfzero['CID']+dfzero['IDSURVEY'].astype(str)

 dfzero = dfzero[dfzero['UNIQ'].isin(covids['UNIQ'])]

 dfmergouter = pd.merge(df,dfzero,left_on='UNIQ',right_on='UNIQ',suffixes=('','_zero'),how='outer',indicator=True)
 dfmerg = pd.merge(df,dfzero,left_on='UNIQ',right_on='UNIQ',suffixes=('','_zero'),how='inner')

 # Get Nuisance Parameters
 alpha,alphaerr,beta,betaerr,massstep,gammaerr = getnuisance(opt)
 alphazero,alphaerrzero,betazero,betaerrzero,massstepzero,gammaerrzero = getnuisance(optzero)

 #HUBBLE FLOW BOOLEAN
 hf=(dfmerg['zHD'].array>zmin)&(dfmerg['zHD'].array<zmax)&(dfmerg['zHD_zero'].array>zmin)&(dfmerg['zHD_zero'].array<zmax)
 dfhf = dfmerg[hf]

 hfzero=(dfmerg['zHD'].array>zmin)&(dfmerg['zHD'].array<zmax)&(dfmerg['zHD_zero'].array>zmin)&(dfmerg['zHD_zero'].array<zmax) #needs to match up with the systematic version
 hfids = dfmerg['CID'][hfzero].array

 #CALIBRATORS BOOLEAN
 cal=dfmerg['CID'].isin(cset).array
 cal_outer=dfmergouter['CID'].isin(cset).array
 calzero=dfmerg['CID_zero'].isin(cset).array
 calzero_outer=dfmergouter['CID_zero'].isin(cset).array

 #PRINT A WARNING FOR MISSING CALIBRATORS HERE!
 missing = dfmergouter[calzero_outer].query('_merge != "both"')
 both = dfmergouter[calzero_outer].query('_merge == "both"')

 print('NUM CALIBRATORS (with duplicates)',len(both),'should be 74')
 if len(both) != 74:
  print('WARNING'*50)

 if len(missing) > 0:
  print('WARNING '*50)
  print(missing[['CID','CID_zero','IDSURVEY','IDSURVEY_zero','MUPULL','MUPULL_zero','MU_zero','MUMODEL_zero']])

 missing = dfmergouter[cal_outer].query('_merge != "both"')

 if len(missing) > 0:
  print('WARNING '*50)
  print(missing[['CID','CID_zero','IDSURVEY','IDSURVEY_zero','MUPULL','MUPULL_zero','MU_zero','MUMODEL_zero']])


 ####### GRABBING PERTINENT COLUMNS ##################################

 mB_hf,mB_cal = dfmerg['mB'].array[hf],dfmerg['mB'].array[cal]
 mB_hfzero,mB_calzero = dfmerg['mB_zero'].array[hfzero],dfmerg['mB_zero'].array[calzero]

 x1_hf,x1_cal = dfmerg['x1'].array[hf],dfmerg['x1'].array[cal]
 x1_hfzero,x1_calzero = dfmerg['x1_zero'].array[hfzero],dfmerg['x1_zero'].array[calzero]

 c_hf,c_cal = dfmerg['c'].array[hf],dfmerg['c'].array[cal]
 c_hfzero,c_calzero = dfmerg['c_zero'].array[hfzero],dfmerg['c_zero'].array[calzero]

 z_hf,z_cal = dfmerg['zHD'].array[hf],dfmerg['zHD'].array[cal]
 z_hfzero,z_calzero = dfmerg['zHD_zero'].array[hfzero],dfmerg['zHD_zero'].array[calzero]

 mucor_hf,mucor_cal = dfmerg['biasCor_mu'].array[hf],dfmerg['biasCor_mu'].array[cal]
 mucor_hfzero,mucor_calzero = dfmerg['biasCor_mu_zero'].array[hfzero],dfmerg['biasCor_mu_zero'].array[calzero]

 muERR_hf,muERR_cal = dfmerg['MUERR'].array[hf],dfmerg['MUERR'].array[cal]
 muERR_hfzero,muERR_calzero = dfmerg['MUERR_zero'].array[hfzero],dfmerg['MUERR_zero'].array[calzero]

 mass_hf,mass_cal = dfmerg['HOST_LOGMASS'].array[hf],dfmerg['HOST_LOGMASS'].array[cal]
 mass_hfzero,mass_calzero = dfmerg['HOST_LOGMASS_zero'].array[hfzero],dfmerg['HOST_LOGMASS_zero'].array[calzero]
 ##############################################################################


 ##### ADD IN THE MASS STEP TO THE mB VALUES ########################

 highmass  =(mass_hf>10)
 mB_hf[highmass]=mB_hf[highmass]+massstep/2.
 lowmass=(mass_hf<10)
 mB_hf[lowmass]=mB_hf[lowmass]-massstep/2.

 highmass  =(mass_hfzero>10)
 mB_hfzero[highmass]=mB_hfzero[highmass]+massstepzero/2.
 lowmass=(mass_hfzero<10)
 mB_hfzero[lowmass]=mB_hfzero[lowmass]-massstepzero/2.

 highmass  =(mass_cal>10)
 mB_cal[highmass]=mB_cal[highmass]+massstep/2.
 lowmass=(mass_cal<10)
 mB_cal[lowmass]=mB_cal[lowmass]-massstep/2.

 highmass  =(mass_calzero>10)
 mB_calzero[highmass]=mB_calzero[highmass]+massstepzero/2.
 lowmass=(mass_calzero<10)
 mB_calzero[lowmass]=mB_calzero[lowmass]-massstepzero/2.

 #####################################################################

 # Polynomial Cosmology
 q0=-0.55
 j0=1
 #cspeed=3*10**5
 cspeed=2.99792458*10**5

 # m0x values for Equation 4 of https://arxiv.org/pdf/1604.01424.pdf
 m0x_values_hf = mB_hf - mucor_hf + alpha*x1_hf - beta*c_hf
 m0x_values_cal = mB_cal - mucor_cal + alpha*x1_cal - beta*c_cal

 m0x_values_hfzero = mB_hfzero - mucor_hfzero + alphazero*x1_hfzero - betazero*c_hfzero
 m0x_values_calzero = mB_calzero - mucor_calzero + alphazero*x1_calzero - betazero*c_calzero


 #apply the fitopt/muopt scaling ########################################################
 m0x_values_hf = m0x_values_hfzero + scale*(m0x_values_hf-m0x_values_hfzero)
 m0x_values_cal = m0x_values_calzero + scale*(m0x_values_cal-m0x_values_calzero)
 ########################################################################################


 # computing means to report for table
 res_m0x_values_hf = np.average(m0x_values_hf-m0x_values_hfzero,weights=1./muERR_hf**2)
 res_m0x_values_cal = np.average(m0x_values_cal-m0x_values_calzero,weights=1./muERR_cal**2)
 ########################################################################################

 # Intercept ############################################################################
 x = cspeed*z_hf*((1 + 0.5*(1-q0)*z_hf) - (1/6)*(1 - q0 - 3*(q0)**2 + 1) * z_hf**2)
 x = np.log10(x)
 y = 0.2 * m0x_values_hf
 e = muERR_hf*.2

 intercept,intercepterr = wmom.wmom(y-x, 1.0/np.power(e,2), calcerr=True)
 intercept=-1*intercept

 xzero = cspeed*z_hfzero*((1 + 0.5*(1-q0)*z_hfzero) - (1/6)*(1 - q0 - 3*(q0)**2 + 1) * z_hfzero**2)
 xzero = np.log10(xzero)
 yzero = 0.2 * m0x_values_hfzero
 ezero = muERR_hfzero*.2

 interceptzero,intercepterrzero = wmom.wmom(yzero-xzero, 1.0/np.power(ezero,2), calcerr=True)
 interceptzero=-1*interceptzero
 #########################################################################################

 ## SOLVING FOR THE CEPHEID DISTANCES THAT GIVE H0 ZERO IN THE NOMINAL ANALYSIS ########
 muzero = -5*np.log10(70.) + 25 + 5*interceptzero + m0x_values_calzero
 ########################################################################################

 #Get new H0 assuming that the nominal analysis gives H0=70 #############################
 H0vec = 10**((m0x_values_cal-muzero + 5*intercept + 25)/5)
 H0 = np.average(H0vec,weights=1/muERR_cal**2)
 ########################################################################################


 cal_scatter = np.std(m0x_values_cal)
 hf_scatter = np.std((y-x))
 #H0 = np.nan

 if np.abs(H0-70)<10**(-6):
  H0res = 0
 else:
  H0res = H0-70

 nsne = len(muERR_hf)
 ncal = len(muERR_cal)
 outdict = {'Name':[name],'OPT':[opt.split('/')[-1]],'ax':[intercept],'ax-Base':[intercept-interceptzero],'H0-Base':[H0res],'H0':[H0],
            'HFm0x-Base':[res_m0x_values_hf],'CALIBm0x-Base':[res_m0x_values_cal],'CALIBstd':[cal_scatter],'HFstd':[hf_scatter],'NHF':[nsne],'NCAL':[ncal]}

 if outdf is None:
  outdf = pd.DataFrame.from_dict(outdict)
 else:
  thisdf = pd.DataFrame.from_dict(outdict)
  outdf = outdf.append(thisdf)

 return outdf



########## RUN ON DIR OF FITOPTS AND MUOPTS ##########################################

bases = [
 #'ForAdam_v3.01_20210412/full/bs20fits/FITOPT000_MUOPT000.FITRES',
 #'ForAdam_v3.01_20210412/full/fits/FITOPT000_MUOPT000.FITRES',
 #'ForAdam_v2.14_20210402/full/fits/FITOPT000_MUOPT000.FITRES',
 #'ForAdam_v2.13_20210329/full/fits/FITOPT000_MUOPT000.FITRES',
 '/full/fits/FITOPT000_MUOPT000.FITRES',
 #'ForAdam_v2.11_20210303/full/fits/FITOPT000_MUOPT000.FITRES',
 #'ForAdam_v2.10_20210225/full/fits/FITOPT000_MUOPT000.FITRES',
 #'ForAdam_v2.9_20210222/full/fits/FITOPT000_MUOPT000.FITRES',
 #'ForAdam_v2.8_20210216/full/fits/FITOPT000_MUOPT000.FITRES',
 #'ForAdam_v2.7_20210215/full/fits/FITOPT000_MUOPT000.FITRES',
 #'ForAdam_v2.6_20201227/full/fits/FITOPT000_MUOPT000.FITRES',
 #'/scratch/midway2/rkessler/PIPPIN_OUTPUT/PLUS/6_BIASCOR/DATASALT3NEWBANDS/output/OUTPUT_BBCFIT/FITOPT000_MUOPT000.FITRES.gz',
 #'ForAdam_v2.5_20201222/full/fits/FITOPT000_MUOPT000.FITRES',
 #'ForAdam_v2.4_20201203/full/fits/FITOPT000_MUOPT000.FITRES',
 #'ForAdam_v2.3_20201130/full/fits/FITOPT000_MUOPT000.FITRES',
 #'ForAdam_v2.2_20201122/full/fits/FITOPT000_MUOPT000.FITRES',
 #'ForAdam_v2.1_20201120/full/fits/FITOPT000_MUOPT000.FITRES',
 #'ForAdam_v2.0_20201112/full/fits/FITOPT000_MUOPT000.FITRES',
]

names = [
 '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/SNIsotropyProject/OriginalSNDataFiles/ForAdam_v2.12_20210304',
 #'ForAdam_v2.0_20201112',
]

outdf = pd.DataFrame()
for name,base in zip(names,bases):
 print ('[name, base] = ' + str([name, base]))
 outdf = h0(name + base,optzero=None,outdf=outdf,scale=1, name = name)

pd.options.display.float_format = '{:,.4f}'.format
pd.set_option('display.max_rows', 1000)

print(outdf[['Name','ax','ax-Base','H0','H0-Base','HFm0x-Base','CALIBm0x-Base','CALIBstd','HFstd','NHF','NCAL']])

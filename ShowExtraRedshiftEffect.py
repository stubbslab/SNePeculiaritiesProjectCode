import cantrips as c
import numpy as np
import loadSN as lsn
import matplotlib.pyplot as plt
import scipy

if __name__=="__main__":
    OmM = 0.3
    OmL = 0.7
    OmR = 0.0
    H0 = 70.0
    Om0 = 1.0
    pull_extinctions = 0
    surveys_of_interest = ['all']
    surveys_to_excise = ['']
    zHD = 0
    z_include_range = [0.0, np.inf]

    all_sn = lsn.loadSN(1, zHD = zHD, surveys_of_interest = surveys_of_interest, pull_extinctions = pull_extinctions, OmM = OmM, OmL = OmL, OmR = OmR, H0 = H0, Om0 = Om0 )
    all_sn = [sn for sn in all_sn if (sn['RA'] > 0.0 or sn['Dec'] > 0.0)]
    all_sn = [sn for sn in all_sn if not(sn['survey'] in surveys_to_excise)]
    all_sn = [sn for sn in all_sn if (sn['z'] > z_include_range[0] and sn['z'] < z_include_range[1])]
    all_zs = [sn['z'] for sn in all_sn]
    all_ras = [sn['RA'] for sn in all_sn]
    all_decs = [sn['Dec'] for sn in all_sn]
    all_mus = [sn['mu'] for sn in all_sn]
    all_muErrs = [sn['muErr'] for sn in all_sn]
    all_muDiffs = [sn['muDiff'] for sn in all_sn]

    MpcTom = 3.086 * (10 ** 16) * (10 ** 6)
    c = 299792458
    H0 = 70 * (10 ** 3.0)/(MpcTom)
    OmM = 0.3
    OmL = 0.7
    dLInMpc = lambda zf, zg: 1/MpcTom * c / H0 * (1 + zf) * scipy.integrate.quad(lambda z_int: 1/np.sqrt((1 + z_int) ** 3 * OmM + (1 + z_int) ** 0 * OmL + (1 + z_int) ** 4 * OmR), 0, (zf - zg)/(zg + 1))[0]
    mu =  lambda zf, zg: 5 * np.log10(dLInMpc(zf, zg)) + 25

    zg_of_z = lambda z_canon, A_z, mu_z, sig_z: A_z * np.exp(- (mu_z - z_canon) ** 2.0 / (2.0 * sig_z ** 2.0))

    z_space = np.linspace(0.01, 3.01, 301)
    #plt.plot(z_space, [mu(z, 0.0) for z in z_space], c = 'k')
    plt.plot(z_space, [mu(z, zg_of_z(z, 0.0, 0.5, 0.05)) - mu(z, 0.0) for z in z_space], c = 'k')
    plt.plot(z_space, [mu(z, zg_of_z(z, 0.1, 0.5, 0.05)) - mu(z, 0.0) for z in z_space], c = 'r')
    plt.scatter(all_zs, all_muDiffs, marker = '.')
    plt.errorbar(all_zs, all_muDiffs, yerr = all_muErrs, fmt = 'none')
    plt.show( )

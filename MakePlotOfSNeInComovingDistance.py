import cantrips as can
import numpy as np
import matplotlib.pyplot as plt
import loadSN as lsn
from astropy.coordinates import SkyCoord
import makePlotOfPS1MDFieldsClass as mpfc

sn_data_type = 'pantheon_plus'
all_sn = lsn.loadSN(1, ['all'], data_type = sn_data_type)
all_zs = [sn['z'] for sn in all_sn]
all_RAs = [sn['RA'] for sn in all_sn]
all_Decs = [sn['Dec'] for sn in all_sn]
sn_coords = [SkyCoord(all_RAs[i],all_Decs[i], frame='icrs', unit='deg') for i in range(len(all_RAs))]
cmb_cold_gal = [209, -57]
cmb_cold_cel = [48.29990842, -20.43730437]
coord = SkyCoord(cmb_cold_gal[0], cmb_cold_gal[1], frame='galactic', unit='deg')
sky_seps = [sn.separation(coord).rad for sn in sn_coords]
my_seps = [can.measureAngularSeparationOnSky(cmb_cold_cel, [all_RAs[i], all_Decs[i]], return_radian = 1) for i in range(len(all_RAs))]
central_z = 0.155
comoving_bin = 250
solid_angle = 20.0

gal_dens_weighting = 0.0 #can be 0.5
z_range = [-0.1, 3.0]
zHD = 1
cutoff_radius_Mpc = np.inf
resid_profile_funct = 'exp_void'
overdensity_param = 200
randomize_all_sn = 0
sdssdir = '/Users/sashabrownsberger/Documents/Harvard/physics/stubbs/SNIsotropyProject/SDSSGalaxies/'
fulldata_fastRead = can.readInColumnsToList('SDSS_fullCoverage_SDSSGals_pzAll.csv', sdssdir, delimiter = ',', n_ignore = 2, all_np_readable = 1)

field_plotter = mpfc.PanStarsFieldManager(1, full_sdss_gal_data_file = 'SDSS_fullCoverage_SDSSGals_pzAll.csv', preloaded_sdss_gals = fulldata_fastRead, gal_dens_weighting = gal_dens_weighting, z_range = z_range, sn_data_type = sn_data_type, zHD = zHD, cutoff_radius_Mpc = cutoff_radius_Mpc, NFW_overdensity_param = overdensity_param, resid_profile_funct = resid_profile_funct, randomize_all_sn = randomize_all_sn)

sn_within_20Deg = [i for i in range(len(all_sn)) if sky_seps[i] < solid_angle * np.pi / 180.0 ]
sn_within_220Mpc = field_plotter.getSNWithinComovingDistance(central_z, cmb_cold_cel, comoving_bin)

zs = np.linspace(0.0, 1.5, 1501)
comoving_dists = [field_plotter.r_of_z(z) for z in zs]
ang_seps  = np.linspace(0.0, np.pi / 2.0, 201)
comoving_mesh, ang_mesh = np.meshgrid(comoving_dists, ang_seps)
z_mesh, ang_mesh = np.meshgrid(zs, ang_seps)
center_comoving = field_plotter.r_of_z(central_z)
comoving_sep_mesh = field_plotter.r_well_from_geometry(center_comoving, comoving_mesh, ang_mesh)

fig, axarr = plt.subplots(2,1, figsize = (10, 8))
ax0 = axarr[0]
ax1 = axarr[1]
z_bounds = [0.05, 0.25]
ang_sep_bounds = [0.0, 30.0 * np.pi / 180.0 ]
c = ax1.contourf(z_mesh, ang_mesh, comoving_sep_mesh, levels = [0, 50, 100, 150, 200, 250, 300, 350, 400])
cbar = fig.colorbar(c)
ax1.scatter(all_zs, np.array(sky_seps), marker = 'x', c = ['r' if i in sn_within_220Mpc else 'g' if i in sn_within_20Deg else 'k' for i in range(len(all_sn))])
ax1.set_xlim(z_bounds)
ax1.set_ylim(ang_sep_bounds)
cbar.set_label('Comoving sep (Mpc)')
ax1.set_xlabel(r'$z$')
ax1.set_ylabel('Ang. sep (radians)')
ax1.axhline(20 * np.pi / 180.0, c = 'g', linestyle = '--' )

scat_ang_sep = ax0.scatter([all_sn[i]['z'] for i in sn_within_20Deg ], [all_sn[i]['muDiff'] for i in sn_within_20Deg ], c = 'g', alpha = 0.5)
ax0.errorbar([all_sn[i]['z'] for i in sn_within_20Deg ], [all_sn[i]['muDiff'] for i in sn_within_20Deg ], yerr = [all_sn[i]['muErr'] for i in sn_within_20Deg ], c = 'g', alpha = 0.5, fmt = 'none')
scat_comoving = ax0.scatter([all_sn[i]['z'] for i in sn_within_220Mpc ], [all_sn[i]['muDiff'] for i in sn_within_220Mpc ], c = 'r', alpha = 0.5)
ax0.errorbar([all_sn[i]['z'] for i in sn_within_220Mpc ], [all_sn[i]['muDiff'] for i in sn_within_220Mpc ], yerr = [all_sn[i]['muErr'] for i in sn_within_220Mpc ],c = 'r', alpha = 0.5, fmt = 'none')
ax0.axvline(central_z, c = 'k', linestyle = '--')
ax0.set_xlim([0.1, 0.5])
ax0.set_ylim([-0.5, 0.5])
ax0.set_xscale('log')
ax0.set_xlabel(r'$z$')
ax0.set_ylabel(r'$\Delta \mu$ (mag)')
ax0.legend([scat_ang_sep, scat_comoving], [r'SNe within $\Omega$ = ' + str(solid_angle) + ' deg', r'SNe within $r$ = ' + str(comoving_bin) + ' Mpc'])
plt.tight_layout()
plt.show()

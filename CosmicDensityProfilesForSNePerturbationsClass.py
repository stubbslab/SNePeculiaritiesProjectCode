import AstronomicalParameterArchive as apa
import numpy as np

class CosmicDensityProfile:

    def radialMassFunctPoint(self, rs, profile_params, extra_params = None):
        """
        The radial mass profile of an point mass, Returns the point mass in
            10^12 M_sun. No extra parameters are required.
        """
        point_mass = self.param_conversion_functs(profile_params)
        M_encs = [point_mass for r in rs]
        return M_encs

    def radialMassFunctNFW(self, rs, profile_params, extra_params = None):
        """
        The radial mass profile of an exponential void, with total excess mass
            0.  Returns mass in 10^12 M_sun.  Note that, because the varied
            parameter is excess over background, this mass function requires
            knowledge of the density at the central redshift.
        """
        if extra_params == None:
            print ('Exponential void mass profile requires you provide 2 additional params: the critical density at the overdensity redshift and OmM at the overdensity redshift. Returning None.')
            return None
        critical_density, OmM_at_z = extra_params
        unitless_overdensity, r_scale = self.param_conversion_functs(profile_params)
        M_encs = np.pi * 4.0 / 3.0 * critical_density * OmM_at_z * unitless_overdensity * rs ** 3.0 * np.exp(-(rs / r_scale) ** 2.0)
        return M_encs

    def radialMassFunctExpVoid(self, rs, profile_params, extra_params = None):
        """
        The radial mass profile of an exponential void, with total excess mass
            0.  Returns mass in 10^12 M_sun.  Note that, because the varied
            parameter is excess over background, this mass function requires
            knowledge of the density at the central redshift.
        Based on function: delta(r) = - delta_0 ( 1 - (2 r^2) / (3 r_0 ^ 2) )  e ^ (- r^2 / r_0^2)
        """
        if extra_params == None:
            print ('Exponential void mass profile requires you provide 2 additional params: the critical density at the overdensity redshift and OmM at the overdensity redshift. Returning None.')
            return None
        critical_density, OmM_at_z = extra_params
        unitless_overdensity, r_scale = self.param_conversion_functs(profile_params)
        M_encs = np.pi * 4.0 / 3.0 * critical_density * OmM_at_z * unitless_overdensity * rs ** 3.0 * np.exp(-(rs / r_scale) ** 2.0)
        
        return M_encs

    def nfw_unitless_over_int(self, radius_per_rs):
        nfw_int = np.log(radius_per_rs + 1) - radius_per_rs / (radius_per_rs+1)
        return nfw_int

    def computeConcentrationForNFW(self, OmM_at_z, critical_density, halo_mass, scale_radius, overdensity_param):
        """
        Compute the concentration parameter, c, under the assumption
           that the halo density falls discontinuously to 0 for radii
           beyond c * r_s (r_s being the NFW halo scale radius).
        """
        concentration = ((3.0 * halo_mass) / (overdensity_param * OmM_at_z * critical_density * 4.0 * np.pi)) ** (1.0 /3.0) / scale_radius
        return concentration

    def radialMassFunctNFW(self, rs, profile_params, extra_params = None):
        """
        The radial mass profile of an NFW halo.  Returns mass in 10^12 M_sun.
            Note that, because the varied parameter is excess over background,
            this mass function requires knowledge of the density at the
            central redshift.
        """
        if extra_params == None:
            print ('NFW halo mass profile requires you provide 1 additional param: the cosmic overdensity param, usually denoted "Delta".')
            return None
        overdensity_param = extra_params[0]
        halo_mass, r_scale, r_cutoff = self.param_conversion_functs(profile_params)
        M_encs = halo_mass * self.nfw_unitless_over_int(r_wells / r_scale )
        return M_encs

    def radialMassFunctUniform(self, rs, profile_params, extra_params = None):
        """
        The radial mass profile of a tophat halo of Uniform mass.  Returns
            mass in 10^12 M_sun.   Note that, because the varied parameter
            is excess over background, this mass function requires knowledge
            of the density at the central redshift.
        """
        halo_mass, r_cutoff = self.param_conversion_functs(profile_params)
        M_encs = np.where(r_wells < cutoff_radius, halo_mass * (r_wells / cutoff_radius) ** 3.0, halo_mass )
        return M_encs

    def getParamInfoDict(self, density_profile_type):
        param_info_dict_dict = { 'nfw': {'n_params':3, 'units':[r'TM$_{\cdot}$', 'Mpc', 'Mpc'], 'plain_strs':['M', 'rs', 'rc'], 'label_strs':[r'$M$', r'$r_s$', r'$r_c$'] },
                                 'point': {'n_params':1, 'units':[r'TM$_{\cdot}$'], 'plain_strs':['M'], 'label_strs':[r'$M$'] },
                                 'uniform':{'n_params':2, 'units':[r'TM$_{\cdot}$', 'Mpc'], 'plain_strs':['M', 'rc'], 'label_strs':[r'$M$', r'$r_c$'] },
                                 'exp_void':{'n_params':2, 'units':['', 'Mpc'], 'plain_strs':['delta0', 'rs'], 'label_strs':[r'$\delta_0$', r'$r_s$'] }}
        return param_info_dict_dict[density_profile_type]

    def getHaloMassFromVariedParam(self, halo_mass_power, linear_scaling_term = 10 ** 3.0):
        """
        We vary a stand-in parameter for the mass.  Here is where
           we convert back to the physical mass, in 10^12 Mpc.
        """
        halo_mass = np.sign(halo_mass_power) * (10.0 ** abs(halo_mass_power) - 1.0) #Allow the sampling to go from -whatever to + whatever, passing through 0 when halo_mass_power = 0
        halo_mass = halo_mass * linear_scaling_term
        return halo_mass

    def getRadiusFromVariedParam(self, radius_power, linear_scaling_term = 10 ** 0.0):
        """
        We vary a stand-in parameter for the radii parameters.  Here is where
           we convert back to the physical mass, in 10^12 Mpc.
        """
        radius = 10.0 ** radius_power
        radius = radius * linear_scaling_term
        return radius

    def getDensityFromVariedParam(self, density_param, linear_scaling_term = 10 ** 0.0):
        """
        We vary a stand-in parameter for the unitless over/under density.  Should
           we cut off the bottom of the scaling to keep the density above -1?
        """
        density = density_param
        density = density * linear_scaling_term
        return density

    def paramConversionFunctExpVoid(self, params, density_scaling = 10.0 ** 0.0, radial_scaling = 10.0 ** 0.0):
        scaled_params = [self.getDensityFromVariedParam(params[0], linear_scaling_term = density_scaling),
                          self.getRadiusFromVariedParam(params[1], linear_scaling_term = radial_scaling)]
        return scaled_params

    def paramConversionFunctUniform(self, params, density_scaling = 10.0 ** 0.0, radial_scaling = 10.0 ** 0.0):
        scaled_params = [self.getHaloMassFromVariedParam(params[0], linear_scaling_term = density_scaling),
                          self.getRadiusFromVariedParam(params[1], linear_scaling_term = radial_scaling)]
        return scaled_params

    def paramConversionFunctNFW(self, params, density_scaling = 10.0 ** 0.0, radial_scaling = 10.0 ** 0.0):
        scaled_params = [self.getHaloMassFromVariedParam(params[0], linear_scaling_term = density_scaling),
                          self.getRadiusFromVariedParam(params[1], linear_scaling_term = radial_scaling),
                           self.getRadiusFromVariedParam(params[2], linear_scaling_term = radial_scaling)]
        return scaled_params

    def paramConversionFunctPoint(self, params, mass_scaling = 10.0 ** 3.0):
        scaled_params = [self.getHaloMassFromVariedParam(params[0], linear_scaling_term = mass_scaling)]
        return scaled_params

    def getNFWProfile(self):
        #Note: for mass (first) and scale radius (second) params, the numbers are logarithmic
        n_profile_fit_params = 3
        param_bounds = [[-10, 10], [0, 4], [0, 4]]
        param_init_guess = [0, 2, 3]
        radial_mass_funct = self.radialMassFunctNFW
        param_conversion_funct = self.paramConversionFunctNFW
        param_info_dict = self.getParamInfoDict('nfw')
        return [n_profile_fit_params, param_bounds, radial_mass_funct, param_conversion_funct, param_info_dict, param_init_guess]

    def getExponentialVoidProfile(self):
        #Note: for mass (first) and scale radius (second) params, the numbers are logarithmic
        n_profile_fit_params = 2
        param_bounds = [[-200, 200], [0, 4]]
        param_init_guess = [0, 2]
        radial_mass_funct = self.radialMassFunctExpVoid
        param_conversion_funct = self.paramConversionFunctExpVoid
        param_info_dict = self.getParamInfoDict('exp_void')
        return [n_profile_fit_params, param_bounds, radial_mass_funct, param_conversion_funct, param_info_dict, param_init_guess]

    def getUniformProfile(self):
        #Note: for mass (first) and scale radius (second) params, the numbers are logarithmic
        n_profile_fit_params = 2
        param_bounds = [[-10, 10], [0, 4]]
        param_init_guess = [0, 2]
        radial_mass_funct = self.radialMassFunctUniform
        param_conversion_funct = self.paramConversionFunctUniform
        param_info_dict = self.getParamInfoDict('uniform')
        return [n_profile_fit_params, param_bounds, radial_mass_funct, param_conversion_funct, param_info_dict, param_init_guess]

    def getPointMassProfile(self):
        #Note: for mass (first) and scale radius (second) params, the numbers are logarithmic
        n_profile_fit_params = 1
        param_bounds = [[-10, 10]]
        param_init_guess = [0.0]
        radial_mass_funct = self.radialMassFunctPoint
        param_conversion_funct = self.paramConversionFunctPoint
        param_info_dict = self.getParamInfoDict('point')
        return [n_profile_fit_params, param_bounds, radial_mass_funct, param_conversion_funct, param_info_dict, param_init_guess]

    def initializeDensityFunction(self, density_profile_type = None):
        if density_profile_type == None:
            density_profile_type = self.density_profile_type
        if density_profile_type == 'nfw':
            self.n_profile_fit_params, self.param_bounds, self.radial_mass_funct, self.param_conversion_functs, self.param_info_dict, self.param_init_guess = self.initializeNFWProfile()
        elif density_profile_type == 'point':
            self.n_profile_fit_params, self.param_bounds, self.radial_mass_funct, self.param_conversion_functs, self.param_info_dict, self.param_init_guess = self.initializePointMassProfile()
        elif density_profile_type == 'uniform':
            self.n_profile_fit_params, self.param_bounds, self.radial_mass_funct, self.param_conversion_functs, self.param_info_dict, self.param_init_guess = self.initializeUniformMassProfile()
        elif density_profile_type == 'exp_void':
            self.n_profile_fit_params, self.param_bounds, self.radial_mass_funct, self.param_conversion_functs, self.param_info_dict, self.param_init_guess = self.getExponentialVoidProfile()
        return 1

    def getProfileFittingPieces(self):
        return self.radial_mass_funct, self.n_profile_fit_params, self.param_bounds,  self.param_conversion_functs, self.param_info_dict, self.param_init_guess

    def __init__(self, density_profile_type, ):
        self.astro_arch = apa.AstronomicalParameterArchive()
        self.density_profile_type = density_profile_type

        self.initializeDensityFunction(self.density_profile_type)

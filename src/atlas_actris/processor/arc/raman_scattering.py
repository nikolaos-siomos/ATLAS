''' This module calculates the (rotational) Raman scattering induced in N2 and O2 molecules in
the atmosphere. It is an application of the method found in:

A. Behrendt and T. Nakamura, "Calculation of the calibration constant of polarization lidar 
and its dependency on atmospheric temperature," Opt. Express, vol. 10, no. 16, pp. 805-817, 2002.
 
The molecular parameters gamma and epsilon are wavelength dependent, and this 
makes the results of the original paper valid only for 532nm. Some new formulas
have been implemented from:

Tomasi, C., Vitale, V., Petkov, B., Lupi, A. & Cacciari, A. Improved 
algorithm for calculations of Rayleigh-scattering optical depth in standard 
atmospheres. Applied Optics 44, 3320 (2005).

and

Chance, K. V. & Spurr, R. J. D. Ring effect studies: Rayleigh scattering, 
including molecular parameters for rotational Raman scattering, and the 
Fraunhofer spectrum. Applied Optics 36, 5224 (1997).

It is not thoroughly tested, so use with care.
'''
import numpy as np
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt

from .constants import hc, k_b, eps_o
from .utilities import number_density_at_pt
from .make_gas import N2, O2, Ar, CO2, H2O
import sys

def rotational_energy(J, molecular_parameters):
    """ Calculates the rotational energy of a homonuclear diatomic molecule for
    quantum number J. The molecule is specified by passing a dictionary with
    parameters.

    Parameters
    ----------
    J : int
       Rotational quantum number.
    molecular_parameters : dict
       A dictionary containing molecular parameters (specifically, B0 and D0).

    Returns
    -------
    E_rot : float
       Rotational energy of the molecule (J)
    """

    B0 = molecular_parameters['B0']
    D0 = molecular_parameters['D0']

    E_rot = (B0 * J * (J + 1) - D0 * J ** 2 * (J + 1) ** 2) * hc
    return E_rot

def raman_shift_stokes(J, molecular_parameters):
    """ Calculates the rotational Raman shift  (delta en) for the Stokes branch for
    quantum number J.

    Parameters
    ----------
    J : int
       Rotational quantum number
    molecular_parameters : dict
       A dictionary containing molecular parameters (specifically, B0 and D0)

    Returns
    -------
    delta_n: float
       Rotational Raman shift [m-1]
    """

    B0 = molecular_parameters['B0']
    D0 = molecular_parameters['D0']

    delta_n = -B0 * 2 * (2 * J + 3) + D0 * (3 * (2 * J + 3) + (2 * J + 3) ** 3)
    return delta_n

def raman_shift_antistokes(J, molecular_parameters):
    """ Calculates the rotational Raman shift (delta en) for the anti-Stokes branch for
    quantum number J.

    Parameters
    ----------
    J: int
       Rotational quantum number
    molecular_parameters: dict
       A dictionary containing molecular parameters (specifically, B0 and D0)

    Returns
    -------
    delta_n: float
       Rotational Raman shift [m-1]
    """
    B0 = molecular_parameters['B0']
    D0 = molecular_parameters['D0']

    # if J > 1:
    #     delta_n = B0 * 2 * (2 * J - 1) - D0 * (3 * (2 * J - 1) + (2 * J - 1) ** 3)
    # else:
    #     delta_n = np.nan
    
    delta_n = B0 * 2 * (2 * J - 1) - D0 * (3 * (2 * J - 1) + (2 * J - 1) ** 3)
        
    return delta_n

# def qm_xsection_stokes(n_incident, J, temperature, molecular_parameters, istotal = False):
#     """ Calculates the rotational Raman backsattering cross section for the Stokes
#     branch for quantum number J at a temperature T.

#     Parameters
#     ----------
#     n_incident : float
#        Wavenumber of incident light [m-1]
       
#     J : int
#        Rotational quantum number
       
#     temperature : float
#        The ambient temperature [K]
       
#     molecular_parameters : dict
#        A dictionary containing molecular parameters.
       
#     istotal : bool
#        A scalar. If set to True then the total scattering cross section
#        is returned instead of the backscattering cross section

#     Returns
#     -------
#     b_s : float
#        Backscattering [m^{2}sr^{-1}] or total scattering cross section [m^{2}]
       
#     """
    
#     B0 = molecular_parameters['B0']

#     # Check if callable or just a number
#     gamma_square_input = molecular_parameters['gamma_square']

#     gamma_square = gamma_square_input  # Assume a float is provided

#     g_index = np.remainder(int(J), 2)

#     g = molecular_parameters['g'][g_index]

#     I = molecular_parameters['I']
    
#     J = float(J)

#     Q_rot = (2.*I + 1.) ** 2 * k_b * temperature / (hc * B0)
    
#     E_factor = np.exp(-rotational_energy(J, molecular_parameters) / (k_b * temperature))
    
#     X = ((J + 1.) * (J + 2.)) / (2.*J + 3.)
    
#     n_shifted = (n_incident + raman_shift_stokes(J, molecular_parameters))
    
#     # Backscattering cross section
#     b_s = 112. * (np.pi ** 4) * g * (n_shifted ** 4) * X * gamma_square * E_factor /\
#         (15. * Q_rot  * (4. * np.pi * eps_o) ** 2) 
    
#     # Converts to the total scattering cross section
#     if istotal:
#         b_s = b_s * (8. * np.pi / 3.) * (10. / 7.)                                

#     return b_s

# def qm_xsection_antistokes(n_incident, J, temperature, molecular_parameters, istotal = False):
#     """ Calculates the rotational Raman backsattering cross section for the Stokes
#     branch for quantum number J at a temperature T.

#     Parameters
#     ----------
#     n_incident : float
#        Wavenumber of incident light [m-1]
       
#     J : int
#        Rotational quantum number
       
#     temperature : float
#        The ambient temperature [K]
       
#     molecular_parameters : dict
#        A dictionary containing molecular parameters.
       
#     istotal : bool
#        A scalar. If set to True then the total scattering cross section
#        is returned instead of the backscattering cross section

#     Returns
#     -------
#     b_s : float
#        Backscattering [m^{2}sr^{-1}] or total scattering cross section [m^{2}]
       
#     """
    
#     B0 = molecular_parameters['B0']

#     # Check if callable or just a number
#     gamma_square_input = molecular_parameters['gamma_square']

#     gamma_square = gamma_square_input  # Assume a float is provided

#     g_index = np.remainder(int(J), 2)

#     g = molecular_parameters['g'][g_index]

#     I = molecular_parameters['I']
    
#     J = float(J)

#     Q_rot = (2.*I + 1.) ** 2 * k_b * temperature / (hc * B0)
    
#     E_factor = np.exp(-rotational_energy(J, molecular_parameters) / (k_b * temperature))
    
#     X = (J * (J - 1.)) / (2.*J - 1.)
    
#     n_shifted = (n_incident + raman_shift_antistokes(J, molecular_parameters))
    
#     # Backscattering cross section
#     b_s = 112. * (np.pi ** 4) * g * (n_shifted ** 4) * X * gamma_square * E_factor /\
#         (15. * Q_rot * (4. * np.pi * eps_o) ** 2) 

#     # Converts to the total scattering cross section
#     if istotal:
#         b_s = b_s * (8. * np.pi / 3.) * (10. / 7.)                                
        
#     return b_s

def qm_xsection_rr_branch(n_incident, J, max_J, temperature, molecular_parameters, branch, istotal = False):
    """ Calculates the rotational Raman backsattering cross section for the Stokes/AntiStokes/Central
    branches for quantum number J at a temperature T.

    Parameters
    ----------
    n_incident : float
       Wavenumber of incident light [m-1]
       
    J : int
       Rotational quantum number
       
    max_J : float
       Maximum rotational quantum number (number of lines considered) 
       
    temperature : float
       The ambient temperature [K]
       
    molecular_parameters : dict
       A dictionary containing molecular parameters.
       
    branch : string
       Select one of Q (central), S (Stokes), O (anti-stokes) 
       
    istotal : bool
       A scalar. If set to True then the total scattering cross section
       is returned instead of the backscattering cross section

    Returns
    -------
    b_s : float
       Backscattering [m^{2}sr^{-1}] or total scattering cross section [m^{2}]
       
    """

    # Check if callable or just a number
    gamma_square_input = molecular_parameters['gamma_square']

    gamma_square = gamma_square_input  # Assume a float is provided

    g_index = np.remainder(int(J), 2)

    g = molecular_parameters['g'][g_index]
    
    J = float(J)

    # Partion function: the sum of all existing rotational states 
    Q_rot = partition_function_by_summing(max_J=max_J, temperature=temperature, molecular_parameters=molecular_parameters)
    # Q_rot = partition_function_ideal_rigid_rotor(temperature=temperature, molecular_parameters=molecular_parameters)

    # Rotational energy that corresponds to quantum number J
    E_factor = np.exp(-rotational_energy(J, molecular_parameters) / (k_b * temperature))
    
    # Conversion factor from CGS to SI units
    CGS_to_SI = (4. * np.pi * eps_o) ** 2
    
    # Population of rotational states corresponding to quantum number J
    P =  g * (2. * J + 1.) * E_factor
    
    # Placzek-Teller coefficients for each branch: The sum equals unity
    if branch == 'Q':
        n_shifted = n_incident
        X = (J * (J + 1.)) / ((2.*J - 1.) * (2. * J + 3.)) 
        
    elif branch == 'S':
        n_shifted = (n_incident + raman_shift_stokes(J, molecular_parameters))
        X = (3. / 2.) * ((J + 1.) * (J + 2.)) / ((2. * J + 1.) * (2. * J + 3.))   
        
    elif branch == 'O':
        n_shifted = (n_incident + raman_shift_antistokes(J, molecular_parameters))

        if J == 0. or J == 1.:
            X = 0.
        if J > 1:
            X = (3. / 2.) * (J * (J - 1.)) / ((2. * J - 1.) * (2. * J + 1.))  
    else:
        sys.exit(f'-- Error! Rotational Raman branch type does not exist ({branch} provided)')        
    
    # b_s = 112. * (np.pi ** 4) * (n_shifted ** 4) * P * X * gamma_square /\
    #     (45. * Q_rot * CGS_to_SI) 
    b_s = (np.pi ** 2) * (n_shifted ** 4) * P * X * 7. * gamma_square /\
        (45. * Q_rot * (eps_o ** 2)) 

    # Converts to the total scattering cross section
    if istotal:
        b_s = b_s * (8. * np.pi / 3.) * (10. / 7.)                                

    return b_s

# def qm_xsection_polarized(n_incident, J, max_J, temperature, molecular_parameters, istotal = False):
#     """ Calculates the backsattering cross section of the isotropically scattered part 
#     of the Cabannes line (excluding the Q branch contribution). Based on Eq. 9 of Adam 2009

#     M. Adam, “Notes on temperature-dependent lidar equations,”
#     J. Atmos. Ocean. Technol. 26, 1021–1039 (2009).
    
#     Parameters
#     ----------
#     n_incident : float
#        Wavenumber of incident light [m-1]
       
#     J : int or 1D array of integers
#        Rotational quantum number
       
#     max_J : float
#        Maximum rotational quantum number (number of lines considered)     
       
#     temperature : float
#        The ambient temperature [K]
       
#     molecular_parameters : dict
#        A dictionary containing molecular parameters.
       
#     istotal : bool
#        A scalar. If set to True then the total scattering cross section
#        is returned instead of the backscattering cross section

#     Returns
#     -------
#     b_s : float
#        Backscattering [m^{2}sr^{-1}] or total scattering cross section [m^{2}]
       
#     """

#     # Check if callable or just a number
#     alpha_square_input = molecular_parameters['alpha_square']

#     alpha_square = alpha_square_input  # Assume a float is provided

#     g_index = np.remainder(int(J), 2)

#     g = molecular_parameters['g'][g_index]
    
#     J = float(J)

#     # Partion function: the sum of all existing rotational states 
#     Q_rot = partition_function_by_summing(max_J=max_J, temperature=temperature, molecular_parameters=molecular_parameters)
#     # Q_rot = partition_function_ideal_rigid_rotor(temperature=temperature, molecular_parameters=molecular_parameters)

#     # Rotational energy that corresponds to quantum number J
#     E_factor = np.exp(-rotational_energy(J, molecular_parameters) / (k_b * temperature))
    
#     # Conversion factor from CGS to SI units
#     CGS_to_SI = (4. * np.pi * eps_o) ** 2
    
#     # Population of rotational states corresponding to quantum number J
#     P =  g * (2. * J + 1.) * E_factor

#     # Converts to the total scattering cross section    
#     # b_s = 112. * (np.pi ** 4) * (n_incident ** 4) * g * (2. * J + 1.) * E_factor * 45. * alpha_square /\
#     #     (45. * Q_rot * 7. * CGS_to_SI) 
#     b_s = (np.pi ** 2) * (n_incident ** 4) * P * alpha_square /\
#         (Q_rot * (eps_o ** 2 ))

#     # Converts to the total scattering cross section
#     if istotal:
#         b_s = b_s * (8. * np.pi / 3.)                                                          

#     return b_s

# def em_xsection_cabannes(n_incident, alpha_square, epsilon, istotal = False):
#     """ Calculates the backsattering cross section of the isotropically scattered part 
#     of the Cabannes line for a spherical molecule or atom (She et al. 2001)
    
#     Parameters
#     ----------
#     n_incident : float
#        Wavenumber of incident light [m-1]
    
#     alpha_square : float or 1D array of floats
#        The isotropic polarizability square (F * m)^2
       
#     epsilon : float or 1D array of floats
#        The square normalized polarizability

#     istotal : bool
#        A scalar. If set to True then the total scattering cross section
#        is returned instead of the backscattering cross section

#     Returns
#     -------
#     b_s : float
#        Backscattering [m^{2}sr^{-1}] or total scattering cross section [m^{2}]
       
#     """
    
#     CGS_to_SI = (4. * np.pi * eps_o) ** 2 
    
#     b_i = 16. * (np.pi ** 4) * (n_incident ** 4) * alpha_square / CGS_to_SI

#     b_q = 16. * (np.pi ** 4) * (n_incident ** 4) * (7. * epsilon / 180.) / CGS_to_SI
    
#     # Converts to the total scattering cross section
#     if istotal:
#         b_i = b_i * (8. * np.pi / 3.)

#         b_q = b_q * (8. * np.pi / 3.) * (10. / 7.)                                                                                          
    
#     b_s = b_i + b_q
    
#     return b_s

# def em_xsection_rayleigh(n_incident, alpha_square, epsilon, istotal = False):
#     """ Calculates the backsattering cross section of the rayleigh spectrum
#     (polarized + depolarized --> Cabanes + O and S branches
#      for a linear molecule or atom (She et al. 2001)
    
#     Parameters
#     ----------
#     n_incident : float
#        Wavenumber of incident light [m-1]
    
#     alpha_square : float or 1D array of floats
#        The polarized polarizability square (F * m)^2
       
#     epsilon : float or 1D array of floats
#        The square normalized polarizability

#     istotal : bool
#        A scalar. If set to True then the total scattering cross section
#        is returned instead of the backscattering cross section

#     Returns
#     -------
#     b_s : float
#        Backscattering [m^{2}sr^{-1}] or total scattering cross section [m^{2}]
       
#     """
    
#     CGS_to_SI = (4. * np.pi * eps_o) ** 2 
    
#     b_i = 16. * (np.pi ** 4) * (n_incident ** 4) * alpha_square / CGS_to_SI

#     b_q = 16. * (np.pi ** 4) * (n_incident ** 4) * (2. * epsilon / 9.) / CGS_to_SI
    
#     # Converts to the total scattering cross section
#     if istotal:
#         b_i = b_i * (8. * np.pi / 3.)

#         b_q = b_q * (8. * np.pi / 3.) * (10. / 7.)                                                                                          
    
#     b_s = b_i + b_q
    
#     return b_s

def em_xsection_polarized(n_incident, alpha_square, istotal = False):
    """ Calculates the backsattering cross section of the isotropically scattered part 
    of the Cabannes line for a molecule or atom (She et al. 2001). The equation
    holds for any kind of molecule (linear, assymetrical etc.)
    
    Parameters
    ----------
    n_incident : float
       Wavenumber of incident light [m-1]
    
    alpha_square : float or 1D array of floats
       The isotropic polarizability square (F * m)^2
       
    epsilon : float or 1D array of floats
       The square normalized polarizability
       
    istotal : bool
       A scalar. If set to True then the total scattering cross section
       is returned instead of the backscattering cross section

    Returns
    -------
    b_s : float or 1D array of floats (same as J)
       Scattering cross section [m^{2}sr^{-1}]
       
    """
    
    CGS_to_SI = (4. * np.pi * eps_o) ** 2 
    
    b_s = (np.pi ** 2) * (n_incident ** 4) * alpha_square / (eps_o ** 2) 

    # Converts to the total scattering cross section
    if istotal:
        b_s = b_s * (8. * np.pi / 3.)     

    return b_s

def partition_function_ideal_rigid_rotor(molecular_parameters, temperature = 293.15):
    """ Asymptotic (analytic) formula for the partition function 
    (see Apendix of Adam 2009). Applicable if the centrifugal distortion D0 = 0

    M. Adam, “Notes on temperature-dependent lidar equations,”
    J. Atmos. Ocean. Technol. 26, 1021–1039 (2009).
    
    Parameters
    ----------
    I : float
       Nuclear spin
       
    molecular_parameters : dict
       A dictionary containing molecular parameters (specifically, B0, I).

    temperature : float
       Gas temperature in Kelvin

    Returns
    -------
    Q_rot : float
       The partition function for J--> infinity
       
    """    
    
    B0 = molecular_parameters['B0']
    I = molecular_parameters['I']
    
    # Analytical partition function, valid for ideal rigid rotor (D = 0)
    Q_rot = (2. * I + 1.) ** 2 * k_b * temperature / (2. * hc * B0)
       
    return(Q_rot)

def partition_function_by_summing(molecular_parameters, temperature = 293.15, max_J = 40):
    """ Maxwell-Boltzman formula for the partition function 
    (see Apendix of Adam 2009, Eq A3). This method is also applied in Long 2002

    M. Adam, “Notes on temperature-dependent lidar equations,”
    J. Atmos. Ocean. Technol. 26, 1021–1039 (2009).
        
    Parameters
    ----------
    max_J : float
       Maximum rotational quantum number (number of lines considered)
       
    molecular_parameters : dict
       A dictionary containing molecular parameters (specifically, g, B0 and D0).

    temperature : float
       Gas temperature in Kelvin

    Returns
    -------
    Q_rot : float
       The partition function for J--> infinity
       
    """

    Js = np.arange(0,max_J,1)

    g_index = np.remainder(Js.astype(int), 2)

    g_table = np.array(molecular_parameters['g'])

    g = g_table[(g_index)]
    
    E_rot = rotational_energy(Js, molecular_parameters)
    
    Q_rot_J = g * (2. * Js + 1.) * np.exp(- E_rot / (k_b * temperature))

    Q_rot = np.nansum(Q_rot_J)
    
    return(Q_rot)
    
def delta_mol_by_fraction(n_incident, transmission_cabannes, molecular_parameters, transmissions):
    """ Calculates the depolarization ratio of the molecular signal detected by a lidar.

    Parameters
    ----------
    n_incident : float
       Wavenumber of the incident wave (m-1)
    molecular_parameters : list
       A list of dictionaries that describe the molecules.
    relative_transmissions : list
       The relative fraction of the intensity of the rotational Raman wings of the molecules detected by the lidar.

    Returns
    -------
    delta_m : float
       The apparent molecular depolarization ratio.
    """
    conc_fraction = np.array([parameter['relative_concentration'] for parameter in molecular_parameters])

    # Check if callable or just a number
    gamma_squares = []
    for parameter in molecular_parameters:
        gamma_square_input = parameter['gamma_square']

        gamma_square = gamma_square_input  # Assume a float is provided

        gamma_squares.append(gamma_square)
    
    gamma_squares = np.array(gamma_squares)

    epsilons = []
    for parameter in molecular_parameters:
        epsilon_input = parameter['epsilon']

        if callable(epsilon_input):
            epsilon = epsilon_input(n_incident)
        else:
            epsilon = epsilon_input  # Assume a float is provided

        epsilons.append(epsilon)
    
    epsilons = np.array(epsilons)
    
    alpha_squares = []
    for parameter in molecular_parameters:
        alpha_squares_input = parameter['alpha_square']

        alpha_squares.append(alpha_squares_input)
    
    alpha_squares = np.array(alpha_squares)

    relative_transmissions = np.array(transmissions) / transmission_cabannes
    # numerator = np.nansum([con * gamma_square * (3. * x / transmission_cabannes + 1.) for (con, gamma_square, x)
    #                     in zip(concentrations, gamma_squares, relative_transmissions)])
    # denominator = np.nansum(
    #     [con *  (3. * x * gamma_square / transmission_cabannes + gamma_square + 45. * alpha_square) for (con, gamma_square, alpha_square, x, epsilon)
    #      in zip(concentrations, gamma_squares, alpha_squares, relative_transmissions, epsilons)])
    # delta_m = 3.0 / 4. * numerator / denominator

    factor = (1. + 3. * relative_transmissions)
    
    # delta_m_i = (3./4.) * epsilons * factor / (45. + epsilons * factor)

    # delta_m_mean = np.nansum(conc_fraction * delta_m_i) / np.nansum(conc_fraction)
    
    # numerator = (3./4.) * epsilons * factor * conc_fraction

    # denominator = conc_fraction * (45. + epsilons * factor)

    numerator = (3./4.) * gamma_squares * factor * conc_fraction

    denominator = conc_fraction * (45. * alpha_squares + gamma_squares * factor)
    
    delta_m = np.nansum(numerator) / np.nansum(denominator)
    
    
    return delta_m

def delta_mol_by_summing(ds_polarized, ds_depolarized):
    
    delta_m = 3. * ds_depolarized / (7. * ds_polarized + 4. * ds_depolarized)
    
    return(delta_m)


class RotationalRaman:

    def __init__(self, wavelength, temperature, max_J=40, N2_parameters=None, O2_parameters=None, Ar_parameters=None, CO2_parameters=None, H2O_parameters=None, optical_filter = None, istotal = False):
        """
        This class calculates the volume depolarization ratio of the molecular
        backscatter signal detected with a polarization lidar.


        Parameters
        ----------
        wavelength: float
           The lidar emission wavelength (nm)
       temperature: float
           The atmospheric temperature (K)
        max_J : int
           The number of Raman lines to consider in each branch.

        """
        
        if istotal == True and optical_filter != None:
            print("-- Warning!: A filter transmittion function was provided but the 'istotal' parameter was set to True. Please not that the filter won't be taken into account as it not relevant to the total scattering cross section in lidar applications.")
            optical_filter = None
        
        if N2_parameters:
            self.N2_parameters = N2_parameters
        else:
            self.N2_parameters = N2(wavelength)

        if O2_parameters:
            self.O2_parameters = O2_parameters
        else:
            self.O2_parameters = O2(wavelength)

        if CO2_parameters:
            self.CO2_parameters = CO2_parameters
        else:
            self.CO2_parameters = CO2(wavelength)

        if Ar_parameters:
            self.Ar_parameters = Ar_parameters
        else:
            self.Ar_parameters = Ar(wavelength)

        if H2O_parameters:
            self.H2O_parameters = H2O_parameters
        else:
            self.H2O_parameters = H2O(wavelength)            
            
        self.optical_filter = optical_filter
        self.temperature = temperature
        self.max_J = max_J
        self.Js = np.arange(0, max_J)
        # self.J_stokes = np.arange(2, max_J)

        self.wavelength = float(wavelength)
        self.wavenumber = 10 ** 9 / self.wavelength # in m-1

        # Calculate the Raman shift for the N2 lines
        self.dn_stokes_N2 = np.array([raman_shift_stokes(J, self.N2_parameters) for J in self.Js])
        self.dn_astokes_N2 = np.array([raman_shift_antistokes(J, self.N2_parameters) for J in self.Js])

        # Convert to wavelegnth
        self.dl_stokes_N2 = 1 / (1 / self.wavelength + np.array(self.dn_stokes_N2) * 10 ** -9)
        self.dl_astokes_N2 = 1 / (1 / self.wavelength + np.array(self.dn_astokes_N2) * 10 ** -9)

        # Calculate the Raman shift for the Ο2 lines
        self.dn_stokes_O2 = np.array([raman_shift_stokes(J, self.O2_parameters) for J in self.Js])
        self.dn_astokes_O2 = np.array([raman_shift_antistokes(J, self.O2_parameters) for J in self.Js])

        # Convert to wavelegnth
        self.dl_stokes_O2 = 1 / (1 / self.wavelength + np.array(self.dn_stokes_O2) * 10 ** -9)
        self.dl_astokes_O2 = 1 / (1 / self.wavelength + np.array(self.dn_astokes_O2) * 10 ** -9)

        # Calculate the Raman shift for the Ο2 lines
        self.dn_stokes_CO2 = np.array([raman_shift_stokes(J, self.CO2_parameters) for J in self.Js])
        self.dn_astokes_CO2 = np.array([raman_shift_antistokes(J, self.CO2_parameters) for J in self.Js])

        # Convert to wavelegnth
        self.dl_stokes_CO2 = 1 / (1 / self.wavelength + np.array(self.dn_stokes_CO2) * 10 ** -9)
        self.dl_astokes_CO2 = 1 / (1 / self.wavelength + np.array(self.dn_astokes_CO2) * 10 ** -9)
        
        self.ds_Q_N2 = np.array([qm_xsection_rr_branch(
            self.wavenumber, J, max_J, temperature, self.N2_parameters, branch='Q', istotal = istotal) for J in self.Js])
        self.ds_stokes_N2 = np.array([qm_xsection_rr_branch(
            self.wavenumber, J, max_J, temperature, self.N2_parameters, branch='S', istotal = istotal) for J in self.Js])
        self.ds_astokes_N2 = np.array([qm_xsection_rr_branch(
            self.wavenumber, J, max_J, temperature, self.N2_parameters, branch='O', istotal = istotal) for J in self.Js])

        self.ds_Q_O2 = np.array([qm_xsection_rr_branch(
            self.wavenumber, J, max_J, temperature, self.O2_parameters, branch='Q', istotal = istotal) for J in self.Js])
        self.ds_stokes_O2 = np.array([qm_xsection_rr_branch(
            self.wavenumber, J, max_J, temperature, self.O2_parameters, branch='S', istotal = istotal) for J in self.Js])
        self.ds_astokes_O2 = np.array([qm_xsection_rr_branch(
            self.wavenumber, J, max_J, temperature, self.O2_parameters, branch='O', istotal = istotal) for J in self.Js])

        self.ds_Q_CO2 = np.array([qm_xsection_rr_branch(
            self.wavenumber, J, max_J, temperature, self.CO2_parameters, branch='Q', istotal = istotal) for J in self.Js])
        self.ds_stokes_CO2 = np.array([qm_xsection_rr_branch(
            self.wavenumber, J, max_J, temperature, self.CO2_parameters, branch='S', istotal = istotal) for J in self.Js])
        self.ds_astokes_CO2 = np.array([qm_xsection_rr_branch(
            self.wavenumber, J, max_J, temperature, self.CO2_parameters, branch='O', istotal = istotal) for J in self.Js])

        self.ds_polar_N2 = em_xsection_polarized(self.wavenumber, self.N2_parameters['alpha_square'], istotal = istotal)
        self.ds_polar_O2 = em_xsection_polarized(self.wavenumber, self.O2_parameters['alpha_square'], istotal = istotal)
        self.ds_polar_CO2 = em_xsection_polarized(self.wavenumber, self.CO2_parameters['alpha_square'], istotal = istotal)
        self.ds_polar_Ar = em_xsection_polarized(self.wavenumber, self.Ar_parameters['alpha_square'], istotal = istotal)
        self.ds_polar_H2O = em_xsection_polarized(self.wavenumber, self.H2O_parameters['alpha_square'], istotal = istotal)

        # self.ds_el_cab_N2 = em_xsection_cabannes(self.wavenumber, self.N2_parameters['alpha_square'], self.N2_parameters['epsilon'], istotal = istotal)
        # self.ds_el_cab_O2 = em_xsection_cabannes(self.wavenumber, self.O2_parameters['alpha_square'], self.O2_parameters['epsilon'], istotal = istotal)
        # self.ds_el_cab_CO2 = em_xsection_cabannes(self.wavenumber, self.CO2_parameters['alpha_square'], self.CO2_parameters['epsilon'], istotal = istotal)
        # self.ds_el_cab_Ar = em_xsection_cabannes(self.wavenumber, self.Ar_parameters['alpha_square'], self.Ar_parameters['epsilon'], istotal = istotal)
        # self.ds_el_cab_H2O = em_xsection_cabannes(self.wavenumber, self.H2O_parameters['alpha_square'], self.H2O_parameters['epsilon'], istotal = istotal)

        # self.ds_el_ray_N2 = em_xsection_rayleigh(self.wavenumber, self.N2_parameters['alpha_square'], self.N2_parameters['epsilon'], istotal = istotal)
        # self.ds_el_ray_O2 = em_xsection_rayleigh(self.wavenumber, self.O2_parameters['alpha_square'], self.O2_parameters['epsilon'], istotal = istotal)
        # self.ds_el_ray_CO2 = em_xsection_rayleigh(self.wavenumber, self.CO2_parameters['alpha_square'], self.CO2_parameters['epsilon'], istotal = istotal)
        # self.ds_el_ray_Ar = em_xsection_rayleigh(self.wavenumber, self.Ar_parameters['alpha_square'], self.Ar_parameters['epsilon'], istotal = istotal)
        # self.ds_el_ray_H2O = em_xsection_rayleigh(self.wavenumber, self.H2O_parameters['alpha_square'], self.H2O_parameters['epsilon'], istotal = istotal)
      

    # def electromagnetic_rayleigh_cross_section(self, pressure=1013.25):
    #     """ Caclulate the backscattering cross section for the 
    #     Rayleigh spectrum(Cabannes + O and S branches)
        
    #     Parameters
    #     ----------
        
    #     pressure: float
    #        The atmospheric pressure [hPa]
    
    #     Returns
    #     -------
        
    #     cross_section:
    #        The backscattering or total scattering cross section of the 
    #        Rayleigh spectrum (Cabannes + O and S branches) for a linear 
    #        molecule or atom. Units are either [m2sr-1] or [m2].  
    
    #        Cross section type depends on istotal given as input to the
    #        RotationalRaman class
           
    #        The calculation for a single molecule is based on 
    #        She et al. 2001, Eq. 6.
           
    #     cross_sections:
    #        Same as cross_sections but given separately for each gas
        
    #     sigma:
    #        The backscattering or total volumetric 
    #        scattering coefficient (crossection * number density)
            
    #     """
        
    #     ds_N2 = np.nansum(self.ds_el_ray_N2) 
        
    #     ds_O2 = np.nansum(self.ds_el_ray_O2) 
        
    #     ds_CO2 = np.nansum(self.ds_el_ray_CO2)

    #     ds_Ar = np.nansum(self.ds_el_ray_Ar) 
    
    #     ds_H2O = np.nansum(self.ds_el_ray_H2O) 
        
    #     cross_section_gas = np.array([ds_N2, ds_O2, ds_Ar, ds_CO2, ds_H2O])
        
    #     N = number_density_at_pt(pressure, self.temperature, relative_humidity=0., ideal=True)
    
    #     c_N2 = self.N2_parameters['relative_concentration']
    #     c_O2 = self.O2_parameters['relative_concentration']
    #     c_Ar = self.Ar_parameters['relative_concentration']
    #     c_CO2 = self.CO2_parameters['relative_concentration']
    #     c_H2O = self.H2O_parameters['relative_concentration']
        
    #     c = np.array([c_N2, c_O2, c_Ar, c_CO2, c_H2O])

    #     gases = {'N2' : '', 'O2' : '', 'Ar' : '', 'CO2' : '', 'H2O' : ''} 
        
    #     cross_section = np.nansum(c * cross_section_gas) 
    
    #     sigma = N * cross_section
        
    #     cross_sections = dict(zip(gases, cross_section_gas))

    #     return(cross_section, cross_sections, sigma)


    # def electromagnetic_cabannes_cross_section(self, pressure = 1013.25):
        
    #     """ Caclulate the backscattering cross section for the cabannes line. 
    
    #     Parameters
    #     ----------
        
    #     pressure: float
    #        The atmospheric pressure [hPa]
           
    #     Returns
    #     -------
        
    #     cross_section:
    #        The backscattering or total scattering cross section of the 
    #        Cabannes line (isotropic + Q branch) for a linear 
    #        molecule or atom. Units are either [m2sr-1] or [m2].  
           
    #        Cross section type depends on istotal given as input to the
    #        RotationalRaman class

    #     cross_sections:
    #        Same as cross_sections but given separately for each gas
        
    #     sigma:
    #        The backscattering or total volumetric 
    #        scattering coefficient (crossection * number density)
            
    #     """
            
    #     ds_N2 = np.nansum(self.ds_el_cab_N2) 
        
    #     ds_O2 = np.nansum(self.ds_el_cab_O2) 
        
    #     ds_CO2 = np.nansum(self.ds_el_cab_CO2)

    #     ds_Ar = np.nansum(self.ds_el_cab_Ar) 
    
    #     ds_H2O = np.nansum(self.ds_el_cab_H2O) 
        
    #     cross_section_gas = np.array([ds_N2, ds_O2, ds_Ar, ds_CO2, ds_H2O])
        
    #     N = number_density_at_pt(pressure, self.temperature, relative_humidity=0., ideal=True)
    
    #     c_N2 = self.N2_parameters['relative_concentration']
    #     c_O2 = self.O2_parameters['relative_concentration']
    #     c_Ar = self.Ar_parameters['relative_concentration']
    #     c_CO2 = self.CO2_parameters['relative_concentration']
    #     c_H2O = self.H2O_parameters['relative_concentration']
        
    #     c = np.array([c_N2, c_O2, c_Ar, c_CO2, c_H2O])

    #     gases = {'N2' : '', 'O2' : '', 'Ar' : '', 'CO2' : '', 'H2O' : ''} 
        
    #     cross_section = np.nansum(c * cross_section_gas) 
    
    #     sigma = N * cross_section
        
    #     cross_sections = dict(zip(gases, cross_section_gas))
    
    #     return(cross_section, cross_sections, sigma)

    def rayleigh_cross_section(self, pressure = 1013.25):

        """ Caclulate the backscattering cross section for the 
        Rayleigh spectrum (Cabannes + O and S branches) 
        by summing the lines
    
        Parameters
        ----------
        
        pressure: float
           The atmospheric pressure [hPa]
        
        Returns
        -------
        
        cross_section:
           The backscattering or total scattering cross section of the 
           Rayleigh spectrum (Cabannes + O and S branches) for a linear 
           molecule or atom. Units are either [m2sr-1] or [m2].  
        
           Cross section type depends on istotal given as input to the
           RotationalRaman class
           
        cross_sections:
           Same as cross_sections but given separately for each gas
        
        sigma:
           The backscattering or total volumetric 
           scattering coefficient (crossection * number density)
            
        """
        
        optical_filter = self.optical_filter
        
        if optical_filter:
            ds_N2 = optical_filter(self.wavelength)*\
                (np.nansum(self.ds_Q_N2) + np.nansum(self.ds_polar_N2)) +\
                    np.nansum(optical_filter(self.dl_stokes_N2) * self.ds_stokes_N2) +\
                        np.nansum(optical_filter(self.dl_astokes_N2) * self.ds_astokes_N2)
                        
            ds_O2 = optical_filter(self.wavelength)*\
                (np.nansum(self.ds_Q_O2) + np.nansum(self.ds_polar_O2)) +\
                    np.nansum(optical_filter(self.dl_stokes_O2) * self.ds_stokes_O2) +\
                        np.nansum(optical_filter(self.dl_astokes_O2) * self.ds_astokes_O2)
                
            ds_CO2 = optical_filter(self.wavelength)*\
                (np.nansum(self.ds_Q_CO2) + np.nansum(self.ds_polar_CO2)) +\
                    np.nansum(optical_filter(self.dl_stokes_CO2) * self.ds_stokes_CO2) +\
                        np.nansum(optical_filter(self.dl_astokes_CO2) * self.ds_astokes_CO2)

            ds_Ar = optical_filter(self.wavelength) * np.nansum(self.ds_polar_Ar) 

            ds_H2O = optical_filter(self.wavelength) * np.nansum(self.ds_polar_H2O) 

        else:
            ds_N2 = np.nansum(self.ds_polar_N2) +  np.nansum(self.ds_stokes_N2) + np.nansum(self.ds_astokes_N2) + np.nansum(self.ds_Q_N2)
            
            ds_O2 = np.nansum(self.ds_polar_O2) + np.nansum(self.ds_stokes_O2) + np.nansum(self.ds_astokes_O2) + np.nansum(self.ds_Q_O2)

            ds_CO2 = np.nansum(self.ds_polar_CO2) + np.nansum(self.ds_stokes_CO2) + np.nansum(self.ds_astokes_CO2) + np.nansum(self.ds_Q_CO2)

            ds_Ar = np.nansum(self.ds_polar_Ar) 
        
            ds_H2O = np.nansum(self.ds_polar_H2O) 
            
        
        cross_section_gas = np.array([ds_N2, ds_O2, ds_Ar, ds_CO2, ds_H2O])
        
        N = number_density_at_pt(pressure, self.temperature, relative_humidity=0., ideal=True)
    
        c_N2 = self.N2_parameters['relative_concentration']
        c_O2 = self.O2_parameters['relative_concentration']
        c_Ar = self.Ar_parameters['relative_concentration']
        c_CO2 = self.CO2_parameters['relative_concentration']
        c_H2O = self.H2O_parameters['relative_concentration']
        
        c = np.array([c_N2, c_O2, c_Ar, c_CO2, c_H2O])

        gases = {'N2' : '', 'O2' : '', 'Ar' : '', 'CO2' : '', 'H2O' : ''} 
        
        cross_section = np.nansum(c * cross_section_gas)

        sigma = N * cross_section
                    
        cross_sections = dict(zip(gases, cross_section_gas))

        return cross_section, cross_sections, sigma


    def rr_cross_section(self, pressure = 1013.25, istotal = False):

        """ Caclulate the backscattering cross section for the 
        RR Wings (O and S branches) by summing the lines
    
        Parameters
        ----------
        
        pressure: float
           The atmospheric pressure [hPa]
        
        Returns
        -------
        
        cross_section:
           The backscattering or total scattering cross section of the 
           RR Wings (O and S branches) for a linear molecule or atom. 
           Units are either [m2sr-1] or [m2].  
        
           Cross section type depends on istotal given as input to the
           RotationalRaman class

        cross_sections:
           Same as cross_sections but given separately for each gas
        
        sigma:
           The backscattering or total volumetric 
           scattering coefficient (crossection * number density)
            
        """
        
        optical_filter = self.optical_filter

        if optical_filter:
            ds_N2 = (np.nansum(optical_filter(self.dl_stokes_N2) * self.ds_stokes_N2) +\
                        np.nansum(optical_filter(self.dl_astokes_N2) * self.ds_astokes_N2)) 
            
            ds_O2 = (np.nansum(optical_filter(self.dl_stokes_O2) * self.ds_stokes_O2) +\
                        np.nansum(optical_filter(self.dl_astokes_O2) * self.ds_astokes_O2))
                
            ds_CO2 = (np.nansum(optical_filter(self.dl_stokes_CO2) * self.ds_stokes_CO2) +\
                        np.nansum(optical_filter(self.dl_astokes_CO2) * self.ds_astokes_CO2))
        else:
            ds_N2 = np.nansum(self.ds_stokes_N2) + np.nansum(self.ds_astokes_N2)
            
            ds_O2 = np.nansum(self.ds_stokes_O2) + np.nansum(self.ds_astokes_O2)

            ds_CO2 = np.nansum(self.ds_stokes_CO2) + np.nansum(self.ds_astokes_CO2)
            
        ds_Ar = 0.
                
        ds_H2O = 0.
        
        cross_section_rr_gas = np.array([ds_N2, ds_O2, ds_Ar, ds_CO2, ds_H2O])
        
        N = number_density_at_pt(pressure, self.temperature, relative_humidity=0., ideal=True)

        c_N2 = self.N2_parameters['relative_concentration']
        c_O2 = self.O2_parameters['relative_concentration']
        c_Ar = self.Ar_parameters['relative_concentration']
        c_CO2 = self.CO2_parameters['relative_concentration']
        c_H2O = self.H2O_parameters['relative_concentration']
        
        c = np.array([c_N2, c_O2, c_Ar, c_CO2, c_H2O])
  
        gases = {'N2' : '', 'O2' : '', 'Ar' : '', 'CO2' : '', 'H2O' : ''} 

        cross_section_rr = np.nansum(c * cross_section_rr_gas)

        sigma_rr = N * cross_section_rr
    
        cross_sections_rr = dict(zip(gases, cross_section_rr_gas))
                    
        return cross_section_rr, cross_sections_rr, sigma_rr 
    
    def stokes_cross_section(self, pressure = 1013.25, istotal = False):

        """ Caclulate the backscattering cross section for the 
        S branch by summing the lines
    
        Parameters
        ----------
        
        pressure: float
           The atmospheric pressure [hPa]
        
        Returns
        -------
        
        cross_section:
           The backscattering or total scattering cross section of the 
           RR Wings (O and S branches) for a linear molecule or atom. 
           Units are either [m2sr-1] or [m2].  
        
           Cross section type depends on istotal given as input to the
           RotationalRaman class

        cross_sections:
           Same as cross_sections but given separately for each gas
        
        sigma:
           The backscattering or total volumetric 
           scattering coefficient (crossection * number density)
            
        """
        
        optical_filter = self.optical_filter

        if optical_filter:
            ds_N2 = np.nansum(optical_filter(self.dl_stokes_N2) * self.ds_stokes_N2) 
            
            ds_O2 = np.nansum(optical_filter(self.dl_stokes_O2) * self.ds_stokes_O2)
                
            ds_CO2 = np.nansum(optical_filter(self.dl_stokes_CO2) * self.ds_stokes_CO2)
        else:
            ds_N2 = np.nansum(self.ds_stokes_N2)
            
            ds_O2 = np.nansum(self.ds_stokes_O2)

            ds_CO2 = np.nansum(self.ds_stokes_CO2)
            
        ds_Ar = 0.
                
        ds_H2O = 0.
        
        cross_section_s_gas = np.array([ds_N2, ds_O2, ds_Ar, ds_CO2, ds_H2O])
        
        N = number_density_at_pt(pressure, self.temperature, relative_humidity=0., ideal=True)

        c_N2 = self.N2_parameters['relative_concentration']
        c_O2 = self.O2_parameters['relative_concentration']
        c_Ar = self.Ar_parameters['relative_concentration']
        c_CO2 = self.CO2_parameters['relative_concentration']
        c_H2O = self.H2O_parameters['relative_concentration']
        
        c = np.array([c_N2, c_O2, c_Ar, c_CO2, c_H2O])
  
        gases = {'N2' : '', 'O2' : '', 'Ar' : '', 'CO2' : '', 'H2O' : ''} 

        cross_section_s = np.nansum(c * cross_section_s_gas)

        sigma_s = N * cross_section_s
    
        cross_sections_s = dict(zip(gases, cross_section_s_gas))
                    
        return cross_section_s, cross_sections_s, sigma_s
    
    def astokes_cross_section(self, pressure = 1013.25, istotal = False):

        """ Caclulate the backscattering cross section for the 
        O branch by summing the lines
    
        Parameters
        ----------
        
        pressure: float
           The atmospheric pressure [hPa]
        
        Returns
        -------
        
        cross_section:
           The backscattering or total scattering cross section of the 
           O branch for a linear molecule or atom. 
           Units are either [m2sr-1] or [m2].  
        
           Cross section type depends on istotal given as input to the
           RotationalRaman class

        cross_sections:
           Same as cross_sections but given separately for each gas
        
        sigma:
           The backscattering or total volumetric 
           scattering coefficient (crossection * number density)
            
        """
        
        optical_filter = self.optical_filter

        if optical_filter:
            ds_N2 = np.nansum(optical_filter(self.dl_astokes_N2) * self.ds_astokes_N2) 
            
            ds_O2 = np.nansum(optical_filter(self.dl_astokes_O2) * self.ds_astokes_O2)
                
            ds_CO2 = np.nansum(optical_filter(self.dl_astokes_CO2) * self.ds_astokes_CO2)
        else:
            ds_N2 = np.nansum(self.ds_astokes_N2)
            
            ds_O2 = np.nansum(self.ds_astokes_O2)

            ds_CO2 = np.nansum(self.ds_astokes_CO2)
            
        ds_Ar = 0.
                
        ds_H2O = 0.
        
        cross_section_o_gas = np.array([ds_N2, ds_O2, ds_Ar, ds_CO2, ds_H2O])
        
        N = number_density_at_pt(pressure, self.temperature, relative_humidity=0., ideal=True)

        c_N2 = self.N2_parameters['relative_concentration']
        c_O2 = self.O2_parameters['relative_concentration']
        c_Ar = self.Ar_parameters['relative_concentration']
        c_CO2 = self.CO2_parameters['relative_concentration']
        c_H2O = self.H2O_parameters['relative_concentration']
        
        c = np.array([c_N2, c_O2, c_Ar, c_CO2, c_H2O])
  
        gases = {'N2' : '', 'O2' : '', 'Ar' : '', 'CO2' : '', 'H2O' : ''} 

        cross_section_o = np.nansum(c * cross_section_o_gas)

        sigma_o = N * cross_section_o
    
        cross_sections_o = dict(zip(gases, cross_section_o_gas))
                    
        return cross_section_o, cross_sections_o, sigma_o

    def Q_cross_section(self, pressure = 1013.25, istotal = False):

        """ Caclulate the backscattering cross section for the 
        Q branch by summing the lines
    
        Parameters
        ----------
        
        pressure: float
           The atmospheric pressure [hPa]
        
        Returns
        -------
        
        cross_section:
           The backscattering or total scattering cross section of the 
           Q branch for a linear molecule or atom. 
           Units are either [m2sr-1] or [m2].  
        
           Cross section type depends on istotal given as input to the
           RotationalRaman class

        cross_sections:
           Same as cross_sections but given separately for each gas
        
        sigma:
           The backscattering or total volumetric 
           scattering coefficient (crossection * number density)
            
        """
        
        optical_filter = self.optical_filter

        if optical_filter:
            ds_N2 = np.nansum(optical_filter(self.dl_Q_N2) * self.ds_Q_N2) 
            
            ds_O2 = np.nansum(optical_filter(self.dl_Q_O2) * self.ds_Q_O2)
                
            ds_CO2 = np.nansum(optical_filter(self.dl_Q_CO2) * self.ds_Q_CO2)
        else:
            ds_N2 = np.nansum(self.ds_Q_N2)
            
            ds_O2 = np.nansum(self.ds_Q_O2)

            ds_CO2 = np.nansum(self.ds_Q_CO2)
            
        ds_Ar = 0.
                
        ds_H2O = 0.
        
        cross_section_q_gas = np.array([ds_N2, ds_O2, ds_Ar, ds_CO2, ds_H2O])
        
        N = number_density_at_pt(pressure, self.temperature, relative_humidity=0., ideal=True)

        c_N2 = self.N2_parameters['relative_concentration']
        c_O2 = self.O2_parameters['relative_concentration']
        c_Ar = self.Ar_parameters['relative_concentration']
        c_CO2 = self.CO2_parameters['relative_concentration']
        c_H2O = self.H2O_parameters['relative_concentration']
        
        c = np.array([c_N2, c_O2, c_Ar, c_CO2, c_H2O])
  
        gases = {'N2' : '', 'O2' : '', 'Ar' : '', 'CO2' : '', 'H2O' : ''} 

        cross_section_q = np.nansum(c * cross_section_q_gas)

        sigma_q = N * cross_section_q
    
        cross_sections_q = dict(zip(gases, cross_section_q_gas))
                    
        return cross_section_q, cross_sections_q, sigma_q

    def cabannes_cross_section(self, pressure = 1013.25):

        """ Caclulate the backscattering cross section for the 
        Cabannes line (isotropic + Q branch) by summing the lines
    
        Parameters
        ----------
        
        pressure: float
           The atmospheric pressure [hPa]
        
        Returns
        -------
        
        cross_section:
           The backscattering or total scattering cross section of the 
           Cabannes line (isotropic + Q branch) for a linear molecule or atom. 
           Units are either [m2sr-1] or [m2].  
        
           Cross section type depends on istotal given as input to the
           RotationalRaman class

        cross_sections:
           Same as cross_sections but given separately for each gas
        
        sigma:
           The backscattering or total volumetric 
           scattering coefficient (crossection * number density)
            
        """
        
        optical_filter = self.optical_filter
        
        if optical_filter:
            ds_N2 = optical_filter(self.wavelength)*\
                (np.nansum(self.ds_Q_N2) + np.nansum(self.ds_polar_N2)) 
        
            ds_O2 = optical_filter(self.wavelength)*\
                (np.nansum(self.ds_Q_O2) + np.nansum(self.ds_polar_O2)) 
                
            ds_CO2 = optical_filter(self.wavelength)*\
                (np.nansum(self.ds_Q_CO2) + np.nansum(self.ds_polar_CO2)) 

            ds_Ar = optical_filter(self.wavelength) * np.nansum(self.ds_polar_Ar) 

            ds_H2O = optical_filter(self.wavelength) * np.nansum(self.ds_polar_H2O) 

        else:
            ds_N2 = np.nansum(self.ds_Q_N2) + np.nansum(self.ds_polar_N2)
            
            ds_O2 = np.nansum(self.ds_Q_O2) + np.nansum(self.ds_polar_O2)

            ds_CO2 = np.nansum(self.ds_Q_CO2) + np.nansum(self.ds_polar_CO2)

            ds_Ar = np.nansum(self.ds_polar_Ar) 
        
            ds_H2O = np.nansum(self.ds_polar_H2O) 
        
        cross_section_cab_gas = np.array([ds_N2, ds_O2, ds_Ar, ds_CO2, ds_H2O])
        
        N = number_density_at_pt(pressure, self.temperature, relative_humidity=0., ideal=True)
    
        c_N2 = self.N2_parameters['relative_concentration']
        c_O2 = self.O2_parameters['relative_concentration']
        c_Ar = self.Ar_parameters['relative_concentration']
        c_CO2 = self.CO2_parameters['relative_concentration']
        c_H2O = self.H2O_parameters['relative_concentration']
        
        c = np.array([c_N2, c_O2, c_Ar, c_CO2, c_H2O])

        gases = {'N2' : '', 'O2' : '', 'Ar' : '', 'CO2' : '', 'H2O' : ''} 
        
        cross_section_cab = np.nansum(c * cross_section_cab_gas)

        sigma_cab = N * cross_section_cab
                    
        cross_sections_cab = dict(zip(gases, cross_section_cab_gas))

        return cross_section_cab, cross_sections_cab, sigma_cab 

    def polarized_cross_section(self, pressure = 1013.25):

        """ Caclulate the backscattering cross section for the 
        isotropic part (electromagnetic calculations prefered over 
        line summation for isotropic scattering)
    
        Parameters
        ----------
        pressure: float
           The atmospheric pressure [hPa]
        
        Returns
        -------
        
        cross_section:
           The backscattering or total scattering cross section of the 
           isotropic part for a linear molecule or atom. 
           Units are either [m2sr-1] or [m2].  
        
           Cross section type depends on istotal given as input to the
           RotationalRaman class

        cross_sections:
           Same as cross_sections but given separately for each gas
        
        sigma:
           The backscattering or total volumetric 
           scattering coefficient (crossection * number density)
            
        """
        
        optical_filter = self.optical_filter
        
        if optical_filter:
            ds_N2 = optical_filter(self.wavelength) * np.nansum(self.ds_polar_N2) 
        
            ds_O2 = optical_filter(self.wavelength) * np.nansum(self.ds_polar_O2) 
                
            ds_CO2 = optical_filter(self.wavelength) * np.nansum(self.ds_polar_CO2) 

            ds_Ar = optical_filter(self.wavelength) * np.nansum(self.ds_polar_Ar) 

            ds_H2O = optical_filter(self.wavelength) * np.nansum(self.ds_polar_H2O) 

        else:
            ds_N2 = np.nansum(self.ds_polar_N2)
            
            ds_O2 = np.nansum(self.ds_polar_O2)

            ds_CO2 = np.nansum(self.ds_polar_CO2)

            ds_Ar = np.nansum(self.ds_polar_Ar) 
        
            ds_H2O = np.nansum(self.ds_polar_H2O) 
        
        cross_section_polar_gas = np.array([ds_N2, ds_O2, ds_Ar, ds_CO2, ds_H2O])
        
        N = number_density_at_pt(pressure, self.temperature, relative_humidity=0., ideal=True)
    
        c_N2 = self.N2_parameters['relative_concentration']
        c_O2 = self.O2_parameters['relative_concentration']
        c_Ar = self.Ar_parameters['relative_concentration']
        c_CO2 = self.CO2_parameters['relative_concentration']
        c_H2O = self.H2O_parameters['relative_concentration']
        
        c = np.array([c_N2, c_O2, c_Ar, c_CO2, c_H2O])

        gases = {'N2' : '', 'O2' : '', 'Ar' : '', 'CO2' : '', 'H2O' : ''} 
        
        cross_section_polar = np.nansum(c * cross_section_polar_gas)

        sigma_polar = N * cross_section_polar
        
        cross_sections_polar = dict(zip(gases, cross_section_polar_gas))

        return cross_section_polar, cross_sections_polar, sigma_polar 

    def depolarized_cross_section(self, pressure = 1013.25):
 
        """ Caclulate the backscattering cross section for the 
        RR spectrum (O, S, and Q branches) by summing the lines
    
        Parameters
        ----------
        
        pressure: float
           The atmospheric pressure [hPa]
        
        Returns
        -------
        
        cross_section:
           The backscattering or total scattering cross section of the 
           RR spectrum (O, S, and Q branches) for a linear molecule or atom. 
           Units are either [m2sr-1] or [m2].  
        
           Cross section type depends on istotal given as input to the
           RotationalRaman class

        cross_sections:
           Same as cross_sections but given separately for each gas
        
        sigma:
           The backscattering or total volumetric 
           scattering coefficient (crossection * number density)
            
        """
        
        optical_filter = self.optical_filter

        if optical_filter:
            ds_N2 = (np.nansum(optical_filter(self.dl_stokes_N2) * self.ds_stokes_N2) +\
                        np.nansum(optical_filter(self.dl_astokes_N2) * self.ds_astokes_N2) +\
                            optical_filter(self.wavelength) * np.nansum(self.ds_Q_N2)) 
            
            ds_O2 = (np.nansum(optical_filter(self.dl_stokes_O2) * self.ds_stokes_O2) +\
                        np.nansum(optical_filter(self.dl_astokes_O2) * self.ds_astokes_O2) +\
                            optical_filter(self.wavelength) * np.nansum(self.ds_Q_O2)) 
                
            ds_CO2 = (np.nansum(optical_filter(self.dl_stokes_CO2) * self.ds_stokes_CO2) +\
                        np.nansum(optical_filter(self.dl_astokes_CO2) * self.ds_astokes_CO2) +\
                            optical_filter(self.wavelength) * np.nansum(self.ds_Q_CO2)) 
        else:
            ds_N2 = np.nansum(self.ds_stokes_N2) + np.nansum(self.ds_astokes_N2) + np.nansum(self.ds_Q_N2)
            
            ds_O2 = np.nansum(self.ds_stokes_O2) + np.nansum(self.ds_astokes_O2) + np.nansum(self.ds_Q_O2)

            ds_CO2 = np.nansum(self.ds_stokes_CO2) + np.nansum(self.ds_astokes_CO2) + np.nansum(self.ds_Q_CO2)
            
        ds_Ar = 0.
                
        ds_H2O = 0.
                
        cross_section_depol_gas = np.array([ds_N2, ds_O2, ds_Ar, ds_CO2, ds_H2O])
        
        N = number_density_at_pt(pressure, self.temperature, relative_humidity=0., ideal=True)
    
        c_N2 = self.N2_parameters['relative_concentration']
        c_O2 = self.O2_parameters['relative_concentration']
        c_Ar = self.Ar_parameters['relative_concentration']
        c_CO2 = self.CO2_parameters['relative_concentration']
        c_H2O = self.H2O_parameters['relative_concentration']
        
        c = np.array([c_N2, c_O2, c_Ar, c_CO2, c_H2O])
        
        gases = {'N2' : '', 'O2' : '', 'Ar' : '', 'CO2' : '', 'H2O' : ''} 

        cross_section_depol = np.nansum(c * cross_section_depol_gas)

        sigma_depol = N * cross_section_depol
        
        cross_sections_depol = dict(zip(gases, cross_section_depol_gas))

        return cross_section_depol, cross_sections_depol, sigma_depol 

    def rr_contribution(self):
        """

        Parameters
        ----------
        temperature
        optical_filter

        Returns
        -------

        """
        
        optical_filter = self.optical_filter

        if optical_filter:
            x_N2 = (np.nansum(optical_filter(self.dl_stokes_N2) * self.ds_stokes_N2) +
                    np.nansum(optical_filter(self.dl_astokes_N2) * self.ds_astokes_N2)) / (np.nansum(self.ds_stokes_N2) + np.nansum(self.ds_astokes_N2))
    
            x_O2 = (np.nansum(optical_filter(self.dl_stokes_O2) * self.ds_stokes_O2) +
                    np.nansum(optical_filter(self.dl_astokes_O2) * self.ds_astokes_O2)) / (np.nansum(self.ds_stokes_O2) + np.nansum(self.ds_astokes_O2))
            
            x_CO2 = (np.nansum(optical_filter(self.dl_stokes_CO2) * self.ds_stokes_CO2) +
                    np.nansum(optical_filter(self.dl_astokes_CO2) * self.ds_astokes_CO2)) / (np.nansum(self.ds_stokes_CO2) + np.nansum(self.ds_astokes_CO2))

        else:
           x_N2 = 1. 

           x_O2 = 1. 
           
           x_CO2 = 1. 

        x_Ar = 1.
        
        x_H2O = 1.
        
        gases = {'N2' : '', 'O2' : '', 'Ar' : '', 'CO2' : '', 'H2O' : ''} 

    
        X = dict(zip(gases, [x_N2, x_O2, x_Ar, x_CO2, x_H2O]))

        return X
    
    # def delta_mol_filter(self):
        

    #     ind_cabannes = np.argmin(np.abs(optical_filter.wavelength - self.wavelength))

    #     self.transmission_cabannes = optical_filter.transmittances[ind_cabannes]

    #     delta_m = delta_mol_by_fraction(self.wavenumber, self.transmission_cabannes, 
    #                                     [self.N2_parameters, self.O2_parameters, self.Ar_parameters, self.CO2_parameters, self.H2O_parameters],
    #                                     [x_N2, x_O2, x_Ar, x_CO2, x_H2O])

    #     return delta_m

    def delta_mol_rayleigh(self,method):

        optical_filter = self.optical_filter

        if optical_filter:
            transmission_cabannes = optical_filter(self.wavelength)
        else:
            transmission_cabannes = 1.
        
        X_all = self.rr_contribution()
                
        X_N2 = X_all['N2']
        X_O2 = X_all['O2']
        X_Ar = X_all['Ar']
        X_CO2 = X_all['CO2']
        X_H2O = X_all['H2O']

        if method == 'fractional_contribution':
            delta_m = delta_mol_by_fraction(n_incident = self.wavenumber, transmission_cabannes = transmission_cabannes, 
                                            molecular_parameters = [self.N2_parameters, self.O2_parameters, self.Ar_parameters, self.CO2_parameters, self.H2O_parameters],
                                            transmissions = [X_N2, X_O2, X_Ar, X_CO2, X_H2O])
            
        elif method == 'line_summation':
            cross_section_polar, _, _ = self.polarized_cross_section()
            cross_section_depol, _, _ = self.depolarized_cross_section()
            delta_m = delta_mol_by_summing(ds_polarized = cross_section_polar, ds_depolarized = cross_section_depol)

        else:
            sys.exit("-- Error: Depolarization ratio method does not exit! Choose one of 'fractional_conribution' or 'line_summation'")
        
        return delta_m
    
    def delta_mol_cabannes(self,method):

        if method == 'fractional_contribution':
            delta_m = delta_mol_by_fraction(n_incident = self.wavenumber, transmission_cabannes = 1., 
                                            molecular_parameters = [self.N2_parameters, self.O2_parameters, self.Ar_parameters, self.CO2_parameters, self.H2O_parameters],
                                            transmissions = 5*[0.])
            
        elif method == 'line_summation':
            cross_section_polar, _, _ = self.polarized_cross_section()
            cross_section_depol, _, _ = self.depolarized_cross_section()
            cross_section_rr, _, _ = self.rr_cross_section()

            delta_m = delta_mol_by_summing(ds_polarized = cross_section_polar, ds_depolarized = cross_section_depol - cross_section_rr)

        return delta_m  

    # def plot_spectrum(self, temperature, molecular_parameters, pressure=None, optical_filter=None, xlim=None, ylim=None,
    #                   figsize=(10, 5)):
    #     fig = plt.figure(figsize=figsize)
    #     ax1 = fig.add_subplot(111)

    #     self.draw_plot(fig, ax1, temperature, molecular_parameters, pressure=pressure, optical_filter=optical_filter,
    #                    xlim=xlim, ylim=ylim)

    #     return fig, ax1

    # def draw_plot(self, fig, ax1, temperature, molecular_parameters, pressure=None, optical_filter=None, xlim=None,
    #               ylim=None,
    #               color='blue', title=None):
    #     # Get the spectral parameters for the specific molecule.
    #     dl_astokes, dl_stokes, ds_astokes, ds_stokes = self._single_molecule_spectra(temperature, molecular_parameters)

    #     # If pressure is provided, calculate backscatter values, else only cross-sections.
    #     if pressure:
    #         N = number_density_at_pt(pressure, temperature, relative_humidity=0, ideal=True) * \
    #             molecular_parameters['relative_concentration']
    #         values_stokes = N * ds_stokes * 1e-4  # from cm2sr-1 to m2sr-1
    #         values_astokes = N * ds_astokes * 1e-4  # from cm2sr-1 to m2sr-1
    #     else:
    #         values_stokes = ds_stokes
    #         values_astokes = ds_astokes

    #     # Define the bar width
    #     x_step = np.diff(dl_stokes)[0]
    #     bar_width = x_step / 4.

    #     bar1 = ax1.bar(dl_stokes, values_stokes, width=bar_width, color=color, label='Full')
    #     ax1.bar(dl_astokes, values_astokes, width=bar_width, color=color)
    #     ax1.legend([bar1[0]], [f'T = {temperature} K'], loc='upper left')

    #     if optical_filter:
    #         bar2 = ax1.bar(dl_stokes, optical_filter(dl_stokes)
    #                        * values_stokes, width=bar_width, color='green', alpha=1, label='Filtered')
    #         ax1.bar(dl_astokes, optical_filter(dl_astokes)
    #                 * values_astokes, width=bar_width, color='green', alpha=1)

    #     ax1.set_xlabel('Wavelength [nm])')

    #     if pressure:
    #         ax1.set_ylabel(r'$\beta$ [$m^{-1}sr^{-1}$]')
    #     else:
    #         ax1.set_ylabel(r'$\left( \frac{d\sigma}{d\omega}\right)_{\pi}$ [$cm^{2}sr^{-1}$]')

    #     if xlim:
    #         plt.xlim(xlim)

    #     if ylim:
    #         plt.ylim(ylim)

    #     if optical_filter:
    #         xmin, xmax = ax1.get_xlim()
    #         filter_label, line_1, ax2 = optical_filter.draw_plot(ax1, xmin, xmax)
    #         ax2.legend([bar1[0], bar2[0], line_1], ['Full', 'Filtered', filter_label], loc=1)

    #     if title is None:
    #         fig.suptitle('Rotational raman spectrumn of $%s$' % molecular_parameters['name'])
    #     else:
    #         fig.suptitle(title)

    # def _single_molecule_spectra(self, temperature, molecular_parameters):
    #     """ Calculate the scattering parameters for a single molecule type, as described by the parameters.
    #     """
    #     ds_stokes = np.array([qm_xsection_stokes(
    #         self.wavenumber, J, temperature, molecular_parameters) for J in self.Js])
    #     ds_astokes = np.array([qm_xsection_antistokes(
    #         self.wavenumber, J, temperature, molecular_parameters) for J in self.Js])
    #     dn_stokes = np.array(
    #         [raman_shift_stokes(J, molecular_parameters) for J in self.Js])
    #     dn_astokes = np.array(
    #         [raman_shift_antistokes(J, molecular_parameters) for J in self.Js])
    #     dl_stokes = 1 / (1 / self.wavelength + dn_stokes * 10 ** -9)
    #     dl_astokes = 1 / (1 / self.wavelength + dn_astokes * 10 ** -9)
    #     return dl_astokes, dl_stokes, ds_astokes, ds_stokes


class BaseFilter:
    """ Base class containing only plotting methods. Actual filter functions should be implemented in the
    subclasses. """

    def plot(self, xmin, xmax):
        fig = plt.figure()
        ax = plt.subplot(111)
        self.draw_plot(ax, xmin, xmax, twin_axis=False)
        plt.draw()
        plt.show()

    def draw_plot(self, main_axis, xmin, xmax, twin_axis=True, color_str='-g', label='Filter'):
        if twin_axis:
            ax = main_axis.twinx()
        else:
            ax = main_axis
        filter_wavelengths = np.linspace(xmin, xmax, 1000)

        filter_efficiency = self(filter_wavelengths)

        line_1, = ax.plot(filter_wavelengths, filter_efficiency, color_str,
                          label=label)
        ax.set_ylabel('Filter efficiency')
        ax.yaxis.label.set_color('green')
        ax.tick_params(axis='y', colors='green')
        ax.set_ylim(0, 1.1)
        return label, line_1, ax

    def __call__(self, *args, **kwargs):
        raise NotImplementedError("Filter functions should be implemented in the subclasses. ")


class GaussianFilter(BaseFilter):
    def __init__(self, wavelength, fwhm, transmittance=1, off_band_transmittance=0):
        '''
        This simple class represents a gausian filter function. To generate
        a new filter use::

           my_filter = FilterFunction(wavelegnth, fwhm)

        with

        wavelegnth - The central wavelength of the filter in nm
        fwhm       - The fwhm of the filter in nm

        If the the filter is called with a wavelegnth (in nm) as an argument
        it will return the  efficiency at this wavelength, for example::

           my_filter = FilterFunction(532, 5)
           my_filter(532) # Will return 1.0
           my_filter(535) # Will return 0.3685
        '''
        self.wavelength = wavelength
        self.fwhm = fwhm
        self.c = fwhm / 2.354820045031

        self.tranmittance = transmittance
        self.off_band_transmittance = off_band_transmittance

    def __call__(self, wavelength):
        value = self.tranmittance * np.exp(-(wavelength - self.wavelength)
                                            ** 2 / (2 * self.c ** 2)) + self.off_band_transmittance
        return value

class LorentzianFilter(BaseFilter):
    def __init__(self, wavelength, fwhm, transmittance=1, off_band_transmittance=0):
        '''
        This simple class represents a gausian filter function. To generate
        a new filter use::

           my_filter = FilterFunction(wavelegnth, fwhm)

        with

        wavelegnth - The central wavelength of the filter in nm
        fwhm       - The fwhm of the filter in nm

        If the the filter is called with a wavelegnth (in nm) as an argument
        it will return the  efficiency at this wavelength, for example::

           my_filter = FilterFunction(532, 5)
           my_filter(532) # Will return 1.0
           my_filter(535) # Will return 0.3685
        '''
        self.wavelength = wavelength
        self.fwhm = fwhm

        self.gamma = fwhm / 2.

        self.tranmittance = transmittance
        self.off_band_transmittance = off_band_transmittance

    def __call__(self, wavelength):
        
        # value_max = 1 / (2. * np.pi * self.fwhm)  
        value = (self.tranmittance) * self.gamma**2 / ((wavelength - self.wavelength)** 2 + self.gamma**2)
        
        return value

class DoubleLorentzianFilter(BaseFilter):
    def __init__(self, wavelength, fwhm, transmittance=1, off_band_transmittance=0):
        '''
        This simple class represents a gausian filter function. To generate
        a new filter use::

           my_filter = FilterFunction(wavelegnth, fwhm)

        with

        wavelegnth - The central wavelength of the filter in nm
        fwhm       - The fwhm of the filter in nm

        If the the filter is called with a wavelegnth (in nm) as an argument
        it will return the  efficiency at this wavelength, for example::

           my_filter = FilterFunction(532, 5)
           my_filter(532) # Will return 1.0
           my_filter(535) # Will return 0.3685
        '''
        self.wavelength = wavelength
        self.fwhm = fwhm

        self.gamma = fwhm / 2.

        self.tranmittance = transmittance
        self.off_band_transmittance = off_band_transmittance

    def __call__(self, wavelength):
        
        value_max = 1 / (2. * np.pi * self.fwhm)  
        value = (self.tranmittance) * self.gamma**2 / ((wavelength - self.wavelength)** 2 + self.gamma**2)
        
        return value


class SquareFilter(BaseFilter):
    def __init__(self, wavelength, width, transmittance=1, off_band_transmittance=0):
        '''
        This simple class represents a square filter function. To generate
        a new filter use::

           my_filter = FilterFunction(wavelegnth, width, transmittacnce)

        '''
        self.wavelength = wavelength
        self.width = width
        self.min_wavelength = wavelength - width / 2.
        self.max_wavelength = wavelength + width / 2.
        self.transmittance = transmittance
        self.off_band_transmittance = off_band_transmittance

    def __call__(self, wavelength):
        w = np.array(wavelength)
        values = np.ones_like(w)
        transmitting_bins = (w > self.min_wavelength) & (w < self.max_wavelength)
        values[transmitting_bins] = self.transmittance
        values[~transmitting_bins] = self.off_band_transmittance
        return values


class CombinedFilter(BaseFilter):
    def __init__(self, filters):
        """
        A combination of several filters. The results will be the multiplication of all filters.

        Parameters
        ----------
        filters : list
           A list of filters
        """
        self.filters = filters

    def __call__(self, wavelength):
        filter_values = [f(wavelength) for f in self.filters]
        values = np.prod(filter_values, axis=0)
        return values


class FileFilter(BaseFilter):
    def __init__(self, file_path, interpolation='linear', off_band_transmittance = 0):
        """
        A filter with transmission given by a text file.

        Currently assumes a simple two-column, tab-delimited file. The first columne should be the wavelength
        in nm and the second the transmittance at the specified wavelength (0 - 1).

        Parameters
        ----------
        file_path : str
           Path to filter function.
        interpolation : str
           The kind of interpolation between provided values. One of 'linear', 'nearest', 'zero', 'slinear',
           'quadratic', 'cubic'. Corresponds to 'kind' argument of interp1d.

        """
        self.file_path = file_path
        self.interpolation = interpolation
        self.off_band_transmittance = off_band_transmittance
        self.read_data()

    def read_data(self):
        # TODO: Make this more flexible using pandas?
        data = np.loadtxt(self.file_path, delimiter='\t')
        self.wavelengths = data[:, 0]
        self.transmittance = data[:, 1]
        self.transmission_function = interp1d(self.wavelengths, self.transmittance, kind=self.interpolation,bounds_error=False,fill_value=self.off_band_transmittance)

    def __call__(self, wavelength):
        return self.transmission_function(wavelength)


class CustomFilter(BaseFilter):
    def __init__(self, wavelengths, transmittances, interpolation='linear', off_band_transmittance = 0):
        """
        A filter with transmission given by a two arrays.

        This is a thin wrapper around numpy's interp1d.

        Parameters
        ----------
        wavelengths : numpy.array
           Wavelength [nm]
        transmittances : numpy.array
            Filter transmittance at the corresponding wavelength (from 0 to 1).
        interpolation : str
            The kind of interpolation between provided values. One of 'linear', 'nearest', 'zero', 'slinear',
            'quadratic', 'cubic'. Corresponds to 'kind' argument of interp1d.
        off_band_transmittance : float scalar
            Fill transmittance value for the regions out of the provided filter spectrum
        """
        
        self.wavelength = wavelengths
        self.transmittances = transmittances
        self.interpolation = interpolation
        self.off_band_transmittance = off_band_transmittance
        self.transimssion_function = interp1d(wavelengths, transmittances, kind=interpolation, bounds_error=False,fill_value=off_band_transmittance)

    def __call__(self, wavelength):
        return self.transimssion_function(wavelength)

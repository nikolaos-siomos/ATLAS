"""
Basic physics constants used in all submodules.

Values taken from http://physics.nist.gov/cuu/Constants/index.html (2014 values)
"""
h = 6.626070040E-34  # plank constant in J s
c = 299792458.  # speed of light in m s-1

# Molar gas constant
#R = 8.3144598  # in J/mol/K --
R = 8.314510  # Value in Ciddor 1996.
# R = 8.31446261815324 #in J/mol/K -- from Wikipedia (actully NIST)

# Avogadro number in mol-1
NA = 6.02214076E23 # Value in Wikipedia

# Boltzmann constant in J/K
k_b = 1.38064852E-23  # Boltzmann constant in J/K
# k_b = 1.3806504 * 10**-23  # J/K  - Folker value
# k_b = R / NA

# Electric Permissivity of vacuum in F⋅m−1
eps_o = 8.8541878128E-12 

# plank constant * speed of light in cm * J
hc = h * c  # m * J

# plank constant * speed of light / Boltzmann constant in cm * K
hc_k = h * c / k_b  # m * K



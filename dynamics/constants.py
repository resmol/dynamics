"""Module containing physical constants."""

from scipy.constants import physical_constants

# Boltzmann constant in au Hartree/Kelvin
BOLTZMANN_JULES = physical_constants["Boltzmann constant"][0]
HARTREE = physical_constants["atomic unit of energy"][0]
BOLTZMANN_AU = BOLTZMANN_JULES / HARTREE

# from femtoseconds to au
AU_TIME = 1e15 * physical_constants["atomic unit of time"][0]

# from bohr to angstrom
BOHR2ANGS = 1e10 * physical_constants["atomic unit of length"][0]

# Boltzmann constant in au hartree/Kelvin
AU_ENERGY = physical_constants["atomic unit of energy"][0]  # J
BOLTZMANN = physical_constants["Boltzmann constant"][0]  # J/K^-1
KB = BOLTZMANN / AU_ENERGY

# Atomic mass
UMA = physical_constants["unified atomic mass unit"][0]  # Kg
AU_MASS = physical_constants["atomic unit of mass"][0]  # Kg
UMA_AUMASS_RATIO = UMA / AU_MASS

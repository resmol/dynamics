"""Run molecular dynamics."""

import logging
from typing import NamedTuple, Optional

import numpy as np
from scipy.constants import physical_constants

from .molecule import Molecule

# Starting logger
logger = logging.getLogger(__name__)

# from femtoseconds to au
AU_TIME = 1e15 * physical_constants["atomic unit of time"][0]

# Boltzmann constant in au Hartree/Kelvin
BOLTZMANN_JULES = physical_constants["Boltzmann constant"][0]
HARTREE = physical_constants["atomic unit of energy"][0]
BOLTZMANN_AU = BOLTZMANN_JULES / HARTREE


class Configuration(NamedTuple):
    """NamedTuple that contains the configuration values to run a simulation."""

    molecule: Molecule
    connectivity: np.ndarray
    hessian: np.ndarray
    gradient: Optional[np.ndarray]
    time: int
    dt: float
    temperature: float


class Thermostat:
    """Contains the state of the termostat.

    For a detailed description of the thermostats chain,
    see https://aip.scitation.org/doi/10.1063/1.2013227
    """

    def __init__(self, mol: Molecule, dt: float,
                 relax_factor: float = 220, temperature: float = 300):
        """Initialize the state of the thermostat.

        Parameters
        ----------
        dt
            Time delta for the integration (au)
        mol
            Molecule to simulate
        relax_factor
            (τ) is characteristic time scale for the system and it is
            usually taken such that τ > 20dt
        temperature
            of the simulation
        """
        self.x1 = 0
        self.vx1 = 0
        self.x2 = 0
        self.vx2 = 0
        self.scale = 1.

        self.dt = dt
        self.temperature = temperature
        freq = AU_TIME / relax_factor  # 1 / [atomic units of time]
        numat = len(mol.atoms)
        # Q Units are Energy * (times ^ 2)
        self.q1 = 3 * numat * temperature * BOLTZMANN_AU / (freq ** 2)
        self.q2 = BOLTZMANN_AU / temperature * (freq ** 2)


def run_simulation(config: Configuration) -> None:
    """Run a MD simulation using a given `config`."""
    mol = config.molecule
    # Initialize the velocities
    mol.generate_random_velocities()

    # convert time delta to AU
    dt = config.dt * AU_TIME

    # Initialize the thermostat
    thermo = Thermostat(dt, mol, temperature)

    # run the MD
    # for t in config.time:

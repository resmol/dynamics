"""Run molecular dynamics."""

import logging
from typing import NamedTuple, Optional

import numpy as np

from .constants import AU_TIME, BOLTZMANN_AU
from .molecule import Molecule

# Starting logger
logger = logging.getLogger(__name__)


class Configuration(NamedTuple):
    """NamedTuple that contains the configuration values to run a simulation."""

    molecule: Molecule
    connectivity: np.ndarray
    hessian: np.ndarray
    gradient: Optional[np.ndarray]
    time: int
    equilibration: int
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
        self.T = temperature
        freq = AU_TIME / relax_factor  # 1 / [atomic units of time]
        self.natoms = len(mol.atoms)
        # Q Units are Energy * (times ^ 2)
        self.Q1 = 3 * self.natoms * temperature * BOLTZMANN_AU / (freq ** 2)
        self.Q2 = BOLTZMANN_AU / temperature * (freq ** 2)

    def update(self, kinetic_energy: float) -> None:
        """Update the thermostat state."""
        dt2 = self.dt * 0.5

        self.update_vx2()
        self.update_vx1(kinetic_energy)
        self.scale = np.exp(-self.vx1 * dt2)

        new_kinetic = kinetic_energy * self.scale ** 2.

        self.update_vx1(new_kinetic)
        self.update_vx2()

    def update_vx1(self, ek: float) -> None:
        """Update the vx1 term."""
        dt4 = self.dt * 0.25
        dt8 = self.dt * 0.125

        self.vx1 += np.exp(-self.vx2 * dt8)
        g1 = (2 * ek - 3 * self.natoms * self.T * BOLTZMANN_AU) / self.Q1
        self.vx1 += g1 * dt4
        self.vx1 *= np.exp(-self.vx2 * dt8)

    def update_vx2(self) -> None:
        """Update the vx2 term."""
        dt4 = self.dt * 0.25

        g2 = (self.Q1 * (self.vx1 ** 2) -
              self.T * BOLTZMANN_AU) / self.Q2
        self.vx2 += g2 * dt4


def run_simulation(config: Configuration) -> None:
    """Run a MD simulation using a given `config`.

    Parameters
    ----------
    config
        Configuration to run the simulation. e.g. temperature, time, etc.

    """
    logger.info("Starting the simulation")
    mol = config.molecule
    # Initialize the velocities
    mol.generate_random_velocities()
    logger.info(
        "Random initial velocities velocities have been generated using a Maxwell-Boltzmann distribution")

    # Initialize the thermostat
    thermo = Thermostat(mol, config.dt, config.temperature)
    logger.info("The thermostat has been initialize")

    # run the MD
    logger.info("Equilibrating the system")
    run_dynamics(mol, config, thermo, step="equilibration")
    logger.info("Running MD simulation")
    run_dynamics(mol, config, thermo)


def run_dynamics(mol: Molecule, config: Configuration, thermo: Thermostat, step="simulation") -> Molecule:
    """
    Run a molecular dynamics for the given molecule an configuration.

    Parameters
    ----------
    mol
        Molecule representing the state
    config
        Configuration to run the simulation
    thermo
        Thermostat state

    Returns
    -------
    Final state of the molecule

    """
    time = config.time if step == "simulation" else config.equilibration
    dt = config.dt
    for t in range(time):

        # Remove translation and rotations
        mol.remove_translations()
        mol.remove_rotations()

        # First phase of Nose-Hoover
        kinetic_energy = mol.compute_kinetic()
        thermo.update(kinetic_energy)

        # Scale the velocities
        mol.velocities *= thermo.scale

        # Update the phase space
        mol.update_positions(dt)
        mol.update_velocities(dt)

        # Compute new gradient
        mol.update_gradient()

        # Second phase of Nose-Hoover
        kinetic_energy = mol.compute_kinetic()
        mol.update_velocities(dt)
        thermo.update(mol, kinetic_energy)

    logger.info(f"MD {step} has finished!")

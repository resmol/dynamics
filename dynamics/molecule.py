"""Molecule definition."""

from pathlib import Path
from typing import Optional, Sequence, Tuple
from scipy.constants import physical_constants

import mendeleev
import numpy as np

from .xyz_reader import parser_xyz

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


class Atom: 
    """Atom definition."""

    def __init__(self, symbol: str, coordinates: Tuple[float, float, float]) -> None:
        self.symbol = symbol.capitalize()
        self.coordinates = np.array(coordinates)


class Molecule:
    """Molecule definition."""

    def __init__(self, atoms: Optional[Sequence[Atom]] = None, file_geometry: Optional[Path] = None) -> None:
        if file_geometry is not None:
            self.read_from_file(file_geometry)
        else:
            self.set_atoms(atoms)

        natoms = len(self.atoms)
        self.velocities = np.zeros(natoms)
        self.masses = np.array(
            [mendeleev.element(atom.symbol).atomic_weight * UMA_AUMASS_RATIO
             for atom in self.atoms])

    def set_atoms(self, atoms: Sequence[Atom]) -> None:
        """Set atoms in object."""
        self.atoms = atoms
        self.coordinates = np.array([at.coordinates for at in atoms])

    def read_from_file(self, path: Path, unit="angstrom") -> None:
        """Read coordinates from file."""
        xs = parser_xyz.parseFile(path)
        labels = (atom.label.capitalize() for atom in xs)
        coordinates = [np.array(atom.xyz) for atom in xs]
        if unit.lower() == "angstrom":
            rs = [x / BOHR2ANGS for x in coordinates]

        atoms = [Atom(symbol, xyz) for symbol, xyz in zip(labels, rs)]
        self.set_atoms(atoms)

    def generate_random_velocities(self, temperature: float = 300) -> None:
        """Use a Maxwell–Boltzmann distribution at `temperature`.

        From Statistical thermodynamics it is known that the velocities,
        of the atoms in a classical system are distributed according
        to the Maxwell-Boltzmann distribution.This says that if the temperature
        of the system is T, the probability of each component of the velocity
        of the ith atom  having a value between v and v + dv is:

        f(V) dV = sqrt(Mass /2 pi Kb T) * exp(-Massi V^2 / 2 Kb T) * dV

        The values of the velocities of the atoms can be assigned by treating them
        as independent Gaussian random variables drawm from the above distribution
        with mean value of 0 and standart deviation of sqrt(KbT/Mass) 
        """
        sigmas = np.sqrt(temperature * KB / self.masses)
        self.velocities = np.random.normal(0, sigmas)

"""Molecule definition."""

from pathlib import Path
from typing import Optional, Sequence, Tuple

import mendeleev
import numpy as np

from .constants import UMA_AUMASS_RATIO, BOHR2ANGS, KB
from .xyz_reader import parser_xyz


class Atom:
    """Atom definition."""

    def __init__(self, symbol: str, coordinates: Tuple[float, float, float]) -> None:
        self.symbol = symbol.capitalize()
        self.coordinates = np.array(coordinates)


class Molecule:
    """Molecule definition."""

    def __init__(self, atoms: Optional[Sequence[Atom]] = None,
                 file_geometry: Optional[Path] = None,
                 velocities: Optional[np.ndarray] = None,
                 gradient: Optional[np.ndarray] = None) -> None:

        if file_geometry is not None:
            self.read_from_file(file_geometry)
        else:
            self.set_atoms(atoms)

        self.natoms = len(self.atoms)

        self.gradient = gradient if gradient is not None else np.zeros(
            (self.natoms, 3))

        self.velocities = velocities if velocities is not None else np.zeros(
            (self.natoms, 3))

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
        labels = [atom.label.capitalize() for atom in xs]
        coordinates = [np.array(atom.xyz[:3]).astype(np.float) for atom in xs]
        if unit.lower() == "angstrom":
            rs = [x / BOHR2ANGS for x in coordinates]

        atoms = [Atom(symbol, xyz) for symbol, xyz in zip(labels, rs)]
        self.set_atoms(atoms)

    def generate_random_velocities(self, temperature: float = 300) -> None:
        """Use a Maxwell–Boltzmann distribution at `temperature`.

        Notes
        -----
        From Statistical thermodynamics it is known that the velocities,
        of the atoms in a classical system are distributed according
        to the Maxwell-Boltzmann distribution.This says that if the temperature
        of the system is T, the probability of each component of the velocity
        of the ith atom  having a value between v and v + dv is:

        f(V) dV = sqrt(Mass /2 pi Kb T) * exp(-Massi V^2 / 2 Kb T) * dV

        The values of the velocities of the atoms can be assigned by treating them
        as independent Gaussian random variables drawm from the above distribution
        with mean value of 0 and standart deviation of sqrt(2 Kb T/Mass) or
        sqrt( 2 Kb T / (3 * Mass)) for the x, y and z components
        """
        sigmas = np.sqrt(2 * temperature * KB / (3 * self.masses))
        self.velocities = np.array(
            [np.random.normal(0, scale=sigmas[0], size=3) for i in range(sigmas.size)])

    def update_positions(self, dt: float) -> None:
        """Update the positions using time step `dt` and the current gradient."""
        first_term = self.velocities * dt
        second_term = self.compute_gradient_term(dt)

        # the minus in the last term is due to the use of the gradient instead of the force
        self.coordinates += first_term - second_term

    def update_velocities(self, dt: float) -> None:
        """Update the velocities using time step `dt` and the current gradient."""
        # the minus in the last term is due to the use of the gradient instead of the force
        self.velocities -= self.compute_gradient_term(dt)

    def compute_gradient_term(self, dt: float) -> np.ndarray:
        """Compute the gradient contribution to the integral."""
        return (0.5 * dt ** 2) * self.gradient / self.masses.reshape(self.natoms, 1)

    def update_gradient(self) -> None:
        """Compute the gradient for this molecular configuration."""
        raise NotImplementedError

    def compute_kinetic(self) -> float:
        """Compute the kinetic energy of the system."""
        squared_velocities = np.sum(self.velocities * self.velocities, axis=1)
        return np.sum(0.5 * self.masses * squared_velocities)

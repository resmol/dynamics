"""User input validation module.

API
---
.. autofunction:: process_input
.. autodata:: schema_dynamics
"""

import logging
from numbers import Real
from typing import Mapping, TypeVar

import numpy as np
import yaml
from schema import Optional, Or, Schema, SchemaError

from .constants import AU_TIME
from .dynamics import Configuration
from .molecule import Molecule


T = TypeVar("T")

# Starting logger
logger = logging.getLogger(__name__)

#: schema to validate the user's input
schema_dynamics = Schema({
    # electronic energy (hartree)
    "energy": Real,

    # molecular geometry in xyz format
    "molecule": str,

    # path to the matrix containing the connectivity
    "connectivity": str,

    # path to the Hessian in cartesian coordinates
    "hessian": str,

    # path Gradient Vector in Cartesians
    Optional("gradient", default=None): Or(str, None),

    # Integration delta in femtoseconds
    Optional("dt", default=1.0): Real,

    # Total simulation time in femtoseconds
    Optional("time", default=1000): int,

    # Temperature in Kelvin
    Optional("temperature", default=300): Real,

    # equilibration time in femtoseconds
    Optional("equilibration", default=100): int
})


def validate_input(input_file: str) -> Configuration:
    """Check that the input provided by the user is correct.

    Parameters
    ----------
    input_file
        path to the user's input

    Returns
    -------
    Configuration
        Configuration to run the given workflow

    Raises
    ------
    SchemaError
        If the input is not valid
    """
    with open(input_file, 'r') as f:
        dict_input = yaml.load(f, Loader=yaml.FullLoader)

    try:
        data = schema_dynamics.validate(dict_input)
        return load_data(data)
    except SchemaError as e:
        msg = f"There was an error in the input yaml provided:\n{e}"
        logger.error(msg)
        raise


def load_data(inp: Mapping[str, T]) -> Configuration:
    """Load the required data for the simulation."""
    # Square matrix with the connectivy
    geometry = Molecule(file_geometry=inp["molecule"])
    connectivity = np.loadtxt(inp["connectivity"]).astype(np.int)
    natoms = int(np.sqrt(connectivity.size))
    connectivity = connectivity.reshape(natoms, natoms)
    hessian = read_hessian(inp["hessian"], natoms)
    if inp["gradient"] is not None:
        grad = np.loadtxt(inp["gradient"]).reshape(natoms, 3)
    else:
        grad = np.zeros((natoms, 3))

    # convert times to AU
    dt_au, time_au, equilibration_au = [
        inp[x] * AU_TIME for x in ("dt", "time", "equilibration")]

    return Configuration(geometry, connectivity, hessian,
                         grad, int(time_au), int(equilibration_au), dt_au, inp["temperature"])


def read_hessian(file_path: str, natoms: int) -> np.ndarray:
    """Read a square matrix from plain ascii text."""
    return np.loadtxt(file_path).reshape(natoms * 3, natoms * 3)

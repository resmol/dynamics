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
    Optional("temperature", default=300): Real
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
    geometry = Molecule(file_geometry=inp["geometry"])
    connectivity = np.loadtxt(inp["connectivity"]).astype(np.int)
    natoms = int(np.sqrt(connectivity.size))
    hessian = np.loadtxt(inp["hessian"]).reshape(natoms, natoms)
    if inp["gradient"] is not None:
        grad = np.loadtxt(inp["gradient"])

    return Configuration(geometry, connectivity, hessian,
                         grad, inp["time"], inp["dt"], inp["temperature"])

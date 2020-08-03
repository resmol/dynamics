"""Command line interface to run the MD simulation."""

import argparse
import numpy as np
from .dynamics import Configuration, run_simulation
from .molecule import Molecule

parser = argparse.ArgumentParser()
parser.add_argument('-e', '--energy', required=True,
                    help="electronic energy (hartree)")
parser.add_argument('-m', '--mol', required=True, help="molecular geometry")
parser.add_argument('-c', '--con', required=True,
                    help="Connectivity matrix")
parser.add_argument('-h', '--hess', required=True,
                    help="Hessian matrix in Cartesians")
parser.add_argument(
    '-g', '--grad', help="Gradient Vector in Cartesians", default=None)
parser.add_argument(
    '-t', '--time', help="total simulation time in femtoseconds", default=100, type=int)
parser.add_argument(
    '-d', '--dt', help="integration time in femtoseconds", default=1.0, type=float)


def main():
    """Parse the command line arguments and run workflow."""
    args = parser.parse_args()
    run_dynamics(args)


def check_configuration(config: Configuration) -> None:
    """Check that the configuration has a physical meaning."""
    pass


def run_dynamics(args: argparse.Namespace) -> None:
    """Run the simulation."""
    # Square matrix with the connectivy
    geometry = Molecule(file_geometry=args.mol)
    connectivity = np.loadtxt(args.con).astype(np.int)
    natoms = int(np.sqrt(connectivity.size))
    hessian = np.loadtxt(args.hess).reshape(natoms, natoms)
    if args.grad is not None:
        grad = np.loadtxt(args.grad)

    config = Configuration(geometry, connectivity,
                           hessian, grad, args.time, args.dt)

    check_configuration(config)                           
    run_simulation(config)


if __name__ == "__main__":
    main()

"""Command line interface to run the MD simulation."""

import argparse
import logging
import pkg_resources
from .dynamics import run_simulation
from .input_validation import validate_input


# Starting logger
logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(description="run_dynamic -i input.yml")
parser.add_argument('-i', required=True,
                    help="Input file in YAML format")


def main():
    """Parse the command line arguments and run workflow."""
    args = parser.parse_args()
    config = validate_input(args.i)
    configure_logger()
    run_simulation(config)


def configure_logger() -> None:
    """Print initial configuration."""
    file_log = 'dynamic.log'
    logging.basicConfig(filename=file_log, level=logging.DEBUG,
                        format='%(asctime)s---%(levelname)s\n%(message)s\n',
                        datefmt='[%I:%M:%S]')
    handler = logging.StreamHandler()
    handler.terminator = ""

    version = pkg_resources.get_distribution('dynamics').version
    path = pkg_resources.resource_filename('dynamics', '')

    logger.info(f"Using dynamics version: {version} ")
    logger.info(f"dynamics path is: {path}")


if __name__ == "__main__":
    main()

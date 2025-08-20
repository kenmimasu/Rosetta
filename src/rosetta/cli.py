#!/usr/bin/env python
################################################################################
import argparse
import sys

from rosetta import settings
from rosetta.internal.parsers import set_subparsers

################################################################################
parser = argparse.ArgumentParser(description=("Main Rosetta command-line executable."))

parser._optionals.title = "Global options"

parser.add_argument(
    "-s",
    "--silent",
    action="store_true",
    help=("Suppress all warnings and take default answers to all questions"),
)
parser.add_argument(
    "-v",
    "--verbose",
    action="store_true",
    help=("Activate verbose setting for program output"),
)

parser.add_argument(
    "--force", action="store_true", help=("Take default answers to all questions")
)
subparsers = parser.add_subparsers(title="Arguments", metavar="INTERFACE", required=True)
set_subparsers(subparsers)


def main():
    args = parser.parse_args()

    # Global session settings
    settings.silent = args.silent
    settings.verbose = args.verbose
    settings.force = args.force

    # Call selected Rosetta interface
    args.interface(args)

if __name__ == "__main__":
    sys.exit(main())
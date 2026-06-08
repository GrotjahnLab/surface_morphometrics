#!/usr/bin/env python
"""DEPRECATED shim -> `morphometrics distances_orientations`.

This script has moved into the installable `surface_morphometrics` package and is
now invoked as `morphometrics distances_orientations ...`. This thin wrapper forwards to that
command and will be removed in a future release.
"""
import sys

from surface_morphometrics.cli import main

print(
    "DEPRECATED: run `morphometrics distances_orientations ...` instead of `python measure_distances_orientations.py ...`",
    file=sys.stderr,
)
main(["distances_orientations", *sys.argv[1:]])

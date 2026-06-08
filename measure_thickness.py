#!/usr/bin/env python
"""DEPRECATED shim -> `morphometrics measure_thickness`.

This script has moved into the installable `surface_morphometrics` package and is
now invoked as `morphometrics measure_thickness ...`. This thin wrapper forwards to that
command and will be removed in a future release.
"""
import sys

from surface_morphometrics.cli import main

print(
    "DEPRECATED: run `morphometrics measure_thickness ...` instead of `python measure_thickness.py ...`",
    file=sys.stderr,
)
main(["measure_thickness", *sys.argv[1:]])

#!/usr/bin/env python
"""DEPRECATED shim -> `morphometrics sample_density`.

This script has moved into the installable `surface_morphometrics` package and is
now invoked as `morphometrics sample_density ...`. This thin wrapper forwards to that
command and will be removed in a future release.
"""
import sys

from surface_morphometrics.cli import main

print(
    "DEPRECATED: run `morphometrics sample_density ...` instead of `python sample_density.py ...`",
    file=sys.stderr,
)
main(["sample_density", *sys.argv[1:]])

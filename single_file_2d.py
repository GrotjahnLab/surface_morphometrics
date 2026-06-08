#!/usr/bin/env python
"""DEPRECATED shim -> `morphometrics hist2d`.

This script has moved into the installable `surface_morphometrics` package and is
now invoked as `morphometrics hist2d ...`. This thin wrapper forwards to that
command and will be removed in a future release.
"""
import sys

from surface_morphometrics.cli import main

print(
    "DEPRECATED: run `morphometrics hist2d ...` instead of `python single_file_2d.py ...`",
    file=sys.stderr,
)
main(["hist2d", *sys.argv[1:]])

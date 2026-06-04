#!/usr/bin/env python
"""DEPRECATED shim -> `morphometrics make_meshes`.

This script has moved into the installable `surface_morphometrics` package and is
now invoked as `morphometrics make_meshes ...`. This thin wrapper forwards to that
command and will be removed in a future release.
"""
import sys

from surface_morphometrics.cli import main

print(
    "DEPRECATED: run `morphometrics make_meshes ...` instead of `python segmentation_to_meshes.py ...`",
    file=sys.stderr,
)
main(["make_meshes", *sys.argv[1:]])

#!/usr/bin/env python
"""DEPRECATED shim -> `morphometrics refine_mesh`.

This script has moved into the installable `surface_morphometrics` package and is
now invoked as `morphometrics refine_mesh ...`. This thin wrapper forwards to that
command and will be removed in a future release.
"""
import sys

from surface_morphometrics.cli import main

print(
    "DEPRECATED: run `morphometrics refine_mesh ...` instead of `python refine_mesh.py ...`",
    file=sys.stderr,
)
main(["refine_mesh", *sys.argv[1:]])

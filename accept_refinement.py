#!/usr/bin/env python
"""DEPRECATED shim -> `morphometrics accept_refinement`.

This script has moved into the installable `surface_morphometrics` package and is
now invoked as `morphometrics accept_refinement ...`. This thin wrapper forwards to that
command and will be removed in a future release.
"""
import sys

from surface_morphometrics.cli import main

print(
    "DEPRECATED: run `morphometrics accept_refinement ...` instead of `python accept_refinement.py ...`",
    file=sys.stderr,
)
main(["accept_refinement", *sys.argv[1:]])

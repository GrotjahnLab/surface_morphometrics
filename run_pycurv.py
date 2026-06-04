#!/usr/bin/env python
"""DEPRECATED shim -> `morphometrics pycurv`.

This script has moved into the installable `surface_morphometrics` package and is
now invoked as `morphometrics pycurv ...`. This thin wrapper forwards to that
command and will be removed in a future release.
"""
import sys

from surface_morphometrics.cli import main

print(
    "DEPRECATED: run `morphometrics pycurv ...` instead of `python run_pycurv.py ...`",
    file=sys.stderr,
)
main(["pycurv", *sys.argv[1:]])

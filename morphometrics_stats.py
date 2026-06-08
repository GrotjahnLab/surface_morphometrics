#!/usr/bin/env python
"""DEPRECATED shim -> `morphometrics stats`.

This script has moved into the installable `surface_morphometrics` package and is
now invoked as `morphometrics stats ...`. This thin wrapper forwards to that
command and will be removed in a future release.
"""
import sys

from surface_morphometrics.cli import main

print(
    "DEPRECATED: run `morphometrics stats ...` instead of `python morphometrics_stats.py ...`",
    file=sys.stderr,
)
main(["stats", *sys.argv[1:]])

#!/usr/bin/env python
"""DEPRECATED shim -> `morphometrics histogram`.

This script has moved into the installable `surface_morphometrics` package and is
now invoked as `morphometrics histogram ...`. This thin wrapper forwards to that
command and will be removed in a future release.
"""
import sys

from surface_morphometrics.cli import main

print(
    "DEPRECATED: run `morphometrics histogram ...` instead of `python single_file_histogram.py ...`",
    file=sys.stderr,
)
main(["histogram", *sys.argv[1:]])

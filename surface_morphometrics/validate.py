#! /usr/bin/env python
"""Validate that a config.yml and its input folders are set up consistently.

Run this before the pipeline to catch the common setup mistakes up front:

  * ``seg_dir`` exists and holds segmentation MRCs (an error if not),
  * ``tomo_dir`` exists (a *warning* only -- tomograms are needed for density
    sampling, thickness, and refinement, but not for the curvature-only pipeline),
  * each segmentation has a matching tomogram under the prefix convention the rest
    of the pipeline uses (the tomogram name is a prefix of the segmentation name),
    so a name mismatch surfaces here instead of as a silently skipped tomogram,
  * which of the configured ``segmentation_values`` labels actually appear in each
    segmentation (flagging segmentations where none of them appear, or where a label
    value is present that the config does not map).

Usage:
  morphometrics validate config.yml
"""

__author__ = "Benjamin Barad"
__email__ = "benjamin.barad@gmail.com"
__license__ = "GPLv3"

import os
import sys
from glob import glob

import click
import numpy as np

from .config_utils import load_config


def _present_label_values(path):
    """Set of integer label values present in a segmentation MRC (0/background dropped).

    Reads slice-by-slice via mmap so a tomogram-sized volume does not have to be held
    in memory at once. Returns None if the data could not be read.
    """
    import mrcfile

    present = set()
    try:
        with mrcfile.mmap(path, permissive=True, mode="r") as mrc:
            data = mrc.data
            if data is None:
                return None
            for plane in (data if data.ndim == 3 else [data]):
                for value in np.unique(plane):
                    rounded = round(float(value))
                    if rounded != 0 and abs(float(value) - rounded) < 1e-3:
                        present.add(int(rounded))
    except Exception:
        return None
    return present


@click.command()
@click.argument("configfile", type=click.Path(exists=True))
def validate_cli(configfile):
    """Check that CONFIGFILE and its seg/tomo folders are set up consistently.

    Verifies the input directories exist, that every segmentation has a matching
    tomogram, and reports which configured labels are present in each segmentation.
    Exits non-zero if there are errors (warnings do not fail).
    """
    config = load_config(configfile, require=("seg_dir", "segmentation_values"))
    seg_dir = config["seg_dir"]
    tomo_dir = config.get("tomo_dir")
    seg_values = config["segmentation_values"]
    value_to_name = {int(v): name for name, v in seg_values.items()}

    errors = 0
    warnings = 0

    def error(msg):
        nonlocal errors
        errors += 1
        click.echo(f"  ERROR:   {msg}")

    def warn(msg):
        nonlocal warnings
        warnings += 1
        click.echo(f"  WARNING: {msg}")

    def ok(msg):
        click.echo(f"  OK:      {msg}")

    click.echo(f"Validating {configfile}\n")

    # --- input folders ---
    click.echo("Folders:")
    seg_dir_ok = os.path.isdir(seg_dir)
    if seg_dir_ok:
        ok(f"seg_dir exists: {seg_dir}")
    else:
        error(f"seg_dir does not exist: {seg_dir}")

    tomo_note = "(needed only for density sampling, thickness, and refinement)"
    tomo_present = bool(tomo_dir) and os.path.isdir(tomo_dir)
    if not tomo_dir:
        warn(f"tomo_dir is not set {tomo_note}")
    elif tomo_present:
        ok(f"tomo_dir exists: {tomo_dir}")
    else:
        warn(f"tomo_dir does not exist: {tomo_dir} {tomo_note}")

    if not seg_dir_ok:
        return _finish(errors, warnings)

    seg_files = sorted(glob(seg_dir + "*.mrc"))
    if not seg_files:
        error(f"no .mrc segmentations found in seg_dir: {seg_dir}")
        return _finish(errors, warnings)

    tomo_bases = ([os.path.basename(t)[:-4] for t in sorted(glob(tomo_dir + "*.mrc"))]
                  if tomo_present else [])

    click.echo("\nConfigured labels: "
               + ", ".join(f"{name}={value}" for name, value in seg_values.items()))
    click.echo(f"\nSegmentations ({len(seg_files)} found):")
    configured = set(value_to_name)

    for seg in seg_files:
        seg_base = os.path.basename(seg)[:-4]
        click.echo(f"\n  {os.path.basename(seg)}")

        # Matching tomogram: its basename is a prefix of the segmentation basename.
        if tomo_present:
            matches = sorted((t for t in tomo_bases if seg_base.startswith(t)),
                             key=len, reverse=True)
            if matches:
                ok(f"matching tomogram: {matches[0]}.mrc")
            else:
                warn(f"no matching tomogram in {tomo_dir} (expected one whose name is a "
                     f"prefix of '{seg_base}'); density/thickness/refinement will skip "
                     "this segmentation -- check for a name mismatch")

        # Which configured labels actually appear in the volume.
        present = _present_label_values(seg)
        if present is None:
            warn("could not read MRC data to check label values")
            continue
        present_configured = sorted(configured & present)
        if present_configured:
            ok("labels present: "
               + ", ".join(f"{value_to_name[v]}({v})" for v in present_configured))
        else:
            warn("none of the configured label values appear in this segmentation -- "
                 "check segmentation_values")
        absent = sorted(configured - present)
        if absent:
            click.echo("           configured but absent here: "
                       + ", ".join(f"{value_to_name[v]}({v})" for v in absent))
        extra = sorted(present - configured)
        if extra:
            warn(f"label value(s) present but not in segmentation_values: {extra}")

    return _finish(errors, warnings)


def _finish(errors, warnings):
    click.echo("\n" + "-" * 60)
    if errors:
        click.echo(f"FAILED: {errors} error(s), {warnings} warning(s).")
        sys.exit(1)
    elif warnings:
        click.echo(f"OK with {warnings} warning(s).")
    else:
        click.echo("All checks passed.")


if __name__ == "__main__":
    validate_cli()

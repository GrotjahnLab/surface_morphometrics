#! /usr/bin/env python
"""Summarize what has been computed for each segmentation in a dataset.

``status config.yml`` scans the work_dir and reports, per segmentation and per
membrane surface, which steps have run: meshing, curvature, refinement (and the
accepted iteration), and which distance / orientation / thickness measurements
exist. Everything is *derived* from the files on disk (and the per-surface CSV
column headers), so the report always reflects the current state -- nothing is
maintained separately and there is no state to get stale.

What each item is read from:
  * mesh        -- {surface}.surface.vtp exists
  * curvature   -- {surface}.AVV_rh{rh}.gt exists
  * refinement  -- the {surface}.accepted_iter marker (written by accept_refinement),
                   or leftover _refined_iter* files (run but not yet accepted)
  * self-dist / verticality / inter / thickness -- columns in {surface}.AVV_rh{rh}.csv

Usage:
  morphometrics status config.yml
  morphometrics status config.yml --output status.txt
"""

__author__ = "Benjamin Barad"
__email__ = "benjamin.barad@gmail.com"
__license__ = "GPLv3"

import datetime
import os
import re
from glob import glob

import click

from .config_utils import load_config

DONE = "✓"   # check mark
NONE = "—"   # em dash


def _csv_columns(path):
    """Column names from a CSV header (reads only the first line)."""
    try:
        with open(path) as handle:
            return set(handle.readline().rstrip("\n").split(","))
    except OSError:
        return set()


def _refine_status(work_dir, surface, radius_hit):
    """Human-readable refinement state for one surface."""
    marker = f"{work_dir}{surface}.accepted_iter"
    if os.path.exists(marker):
        try:
            return f"accepted@iter{open(marker).read().strip()}"
        except OSError:
            return "accepted"
    iters = glob(f"{work_dir}{surface}_refined_iter*.surface.vtp")
    if iters:
        steps = sorted(int(m.group(1)) for m in
                       (re.search(r"_refined_iter(\d+)\.surface\.vtp$", p) for p in iters) if m)
        span = f"iter{steps[0]}-{steps[-1]}" if len(steps) > 1 else f"iter{steps[0]}"
        return f"run ({span}, not accepted)"
    if os.path.exists(f"{work_dir}{surface}.AVV_rh{radius_hit}.gt.orig.bak"):
        return "accepted (step unknown)"
    return NONE


@click.command()
@click.argument("configfile", type=click.Path(exists=True))
@click.option("-o", "--output", "output_path", default=None, type=click.Path(),
              help="Also write the summary to this text file.")
def status_cli(configfile, output_path):
    """Report what has been computed for each segmentation (derived from work_dir)."""
    config = load_config(configfile, require=("seg_dir", "work_dir", "segmentation_values"))
    seg_dir = config["seg_dir"]
    work_dir = config["work_dir"]
    tomo_dir = config.get("tomo_dir")
    labels = list(config["segmentation_values"])
    radius_hit = config["curvature_measurements"]["radius_hit"]

    lines = []

    def emit(text=""):
        lines.append(text)

    emit(f"Status {NONE} {work_dir}  (radius_hit {radius_hit})        "
         f"{datetime.date.today().isoformat()}")
    emit()

    seg_files = sorted(glob(f"{seg_dir}*.mrc"))
    tomo_bases = ([os.path.basename(t)[:-4] for t in sorted(glob(f"{tomo_dir}*.mrc"))]
                  if tomo_dir and os.path.isdir(tomo_dir) else [])

    if not seg_files:
        emit(f"(no segmentations found in {seg_dir})")

    for seg in seg_files:
        seg_base = os.path.basename(seg)[:-4]
        tnote = ""
        if tomo_bases:
            tnote = (f"   [tomogram {DONE}]" if any(seg_base.startswith(t) for t in tomo_bases)
                     else "   [tomogram missing]")
        emit(f"{seg_base}{tnote}")

        meshed = []
        for label in labels:
            surface = f"{seg_base}_{label}"
            if not os.path.exists(f"{work_dir}{surface}.surface.vtp"):
                continue
            meshed.append(label)
            curv = os.path.exists(f"{work_dir}{surface}.AVV_rh{radius_hit}.gt")
            cols = _csv_columns(f"{work_dir}{surface}.AVV_rh{radius_hit}.csv") if curv else set()
            refine = _refine_status(work_dir, surface, radius_hit)
            partners = [o for o in labels if o != label and f"{o}_dist" in cols]
            emit(f"  {label:<8} mesh {DONE}  curv {DONE if curv else NONE}  "
                 f"refine {refine}  self-dist {DONE if 'self_dist_min' in cols else NONE}  "
                 f"vert {DONE if 'verticality' in cols else NONE}  "
                 f"inter {', '.join(partners) if partners else NONE}  "
                 f"thickness {DONE if 'thickness' in cols else NONE}")

        if not meshed:
            emit("  (no meshes yet)")
        else:
            missing = [label for label in labels if label not in meshed]
            if missing:
                emit(f"  (no mesh for: {', '.join(missing)})")
        emit()

    emit(f"legend:  {DONE} done {NONE} not done   inter: partner surfaces with "
         "distance/orientation measured")

    report = "\n".join(lines)
    click.echo(report)
    if output_path:
        with open(output_path, "w") as handle:
            handle.write(report + "\n")
        click.echo(f"\nWrote {output_path}")


if __name__ == "__main__":
    status_cli()

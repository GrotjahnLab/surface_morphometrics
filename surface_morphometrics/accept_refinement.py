#! /usr/bin/env python
"""Accept a mesh refinement iteration as the new working surface.

`refine_mesh.py` writes a numbered set of intermediate surfaces
(`{basename}_refined_iter{N}.*`) plus a handful of diagnostic plots and a few
summary files.  This script "commits" one of those iterations: it backs up the
original (pre-refinement) surfaces, promotes the chosen iteration to be the new
main surface that the rest of the pipeline consumes, and cleans up the leftover
intermediate files so downstream steps don't pick them up by accident.

Why the cleanup matters: the thickness and distance steps discover surfaces by
globbing `{basename}*{component}.AVV_rh{radius_hit}.gt`.  The per-iteration
graphs (`{basename}_refined_iter{N}.AVV_rh{N}.gt`) match that pattern, so if
they are left behind they would be processed as duplicate surfaces.  This script
removes them, keeping only the accepted surface and the refinement summaries.

Usage:
  accept_refinement.py config.yml STEP            # accept iteration STEP for all refined surfaces
  accept_refinement.py config.yml STEP --component OMM
  accept_refinement.py config.yml STEP --dry-run  # show what would happen without touching files
"""

__author__ = "Benjamin Barad"
__email__ = "benjamin.barad@gmail.com"
__license__ = "GPLv3"

import os
import re
from glob import glob
from pathlib import Path

import click
import yaml


def _iterations_by_surface(work_dir):
    """Map each surface basename to its sorted available refinement iterations.

    Refinement stops early once a surface converges, so different surfaces may
    have different numbers of `_refined_iter*` files.
    """
    pattern = re.compile(r"^(.*)_refined_iter(\d+)\.surface\.vtp$")
    out = {}
    for path in glob(f"{work_dir}*_refined_iter*.surface.vtp"):
        match = pattern.match(os.path.basename(path))
        if match:
            out.setdefault(match.group(1), []).append(int(match.group(2)))
    return {base: sorted(iters) for base, iters in out.items()}


# The canonical surface files the downstream pipeline consumes.  Only these are
# backed up (from the original) and promoted (from the accepted iteration).
# {rh} is filled in with the configured radius_hit.
CANONICAL_EXTENSIONS = [
    ".surface.vtp",
    ".AVV_rh{rh}.gt",
    ".AVV_rh{rh}.vtp",
    ".AVV_rh{rh}.csv",
]

# Regenerable pycurv / refinement intermediates.  These are removed for the
# accepted surface so only the canonical files above remain (rebuilt by pycurv
# if it is re-run).  {rh} is filled in with the configured radius_hit.
INTERMEDIATE_EXTENSIONS = [
    ".NVV_rh{rh}.gt",
    ".scaled_cleaned.gt",
    ".scaled_cleaned.vtp",
    ".lightweight.gt",
    ".lightweight_sampling.csv",
    ".AVV_rh{rh}_sampling.csv",
]

# Summary files written once per surface by refine_mesh.py.  These are kept.
SUMMARY_SUFFIXES = [
    "_refinement_stats.csv",
    "_refinement_convergence.png",
    "_profile_evolution.png",
]


def _move(src, dst, dry_run):
    """Rename src -> dst, announcing the action (or what would happen)."""
    prefix = "  [dry-run] would rename" if dry_run else "  renamed"
    print(f"{prefix}: {os.path.basename(src)} -> {os.path.basename(dst)}")
    if not dry_run:
        os.replace(src, dst)


def _remove(path, dry_run):
    """Delete path, announcing the action (or what would happen)."""
    prefix = "  [dry-run] would delete" if dry_run else "  deleted"
    print(f"{prefix}: {os.path.basename(path)}")
    if not dry_run:
        os.remove(path)


def accept_one(work_dir, basename, step, radius_hit, dry_run):
    """Accept iteration `step` for a single surface basename.

    Returns one of:
      None  -> no such iteration, skipped
      True  -> accepted and the surface has a curvature graph, ready for downstream
      False -> accepted but the iteration was lightweight (no curvature graph);
               pycurv must be re-run on the promoted surface before continuing
    """
    iter_prefix = f"{basename}_refined_iter{step}"
    accepted_surface = f"{work_dir}{iter_prefix}.surface.vtp"
    if not os.path.exists(accepted_surface):
        print(f"  SKIP: no iteration {step} surface found for {basename} "
              f"(expected {os.path.basename(accepted_surface)})")
        return None

    print(f"Accepting iteration {step} for: {basename}")

    # refine_mesh.py only runs full pycurv (producing an AVV curvature graph) on
    # the final iteration; intermediate xcorr iterations get a lightweight graph
    # instead.  If the chosen iteration has no AVV graph, the promoted surface is
    # not ready for downstream analysis until pycurv is re-run.
    has_avv = os.path.exists(f"{work_dir}{iter_prefix}.AVV_rh{radius_hit}.gt")

    # 1. Back up the canonical original surfaces to *.orig.bak so the promotion
    #    below doesn't clobber them and the original remains recoverable.
    for ext_template in CANONICAL_EXTENSIONS:
        original = f"{work_dir}{basename}{ext_template.format(rh=radius_hit)}"
        if os.path.exists(original):
            backup = f"{original}.orig.bak"
            if os.path.exists(backup):
                print(f"  NOTE: backup already exists, leaving original in place: "
                      f"{os.path.basename(backup)}")
                continue
            _move(original, backup, dry_run)

    # 2. Promote the accepted iteration's canonical files to the main surface
    #    names by stripping the "_refined_iter{step}" tag.
    promoted = set()
    for ext_template in CANONICAL_EXTENSIONS:
        ext = ext_template.format(rh=radius_hit)
        src = f"{work_dir}{iter_prefix}{ext}"
        if os.path.exists(src):
            _move(src, f"{work_dir}{basename}{ext}", dry_run)
            promoted.add(src)

    # 3. Clean up everything else: other iterations, per-iteration diagnostic
    #    plots, and the regenerable pycurv/refinement intermediates for this
    #    surface.  Keep the canonical files (step 2), the *.orig.bak backups,
    #    and the per-surface refinement summaries.
    leftovers = set()
    leftovers.update(glob(f"{work_dir}{basename}_refined_iter*"))
    leftovers.update(glob(f"{work_dir}{basename}_samples_iter*.png"))
    leftovers.update(glob(f"{work_dir}{basename}_profile_iter*.png"))
    for ext_template in INTERMEDIATE_EXTENSIONS:
        leftovers.update(glob(f"{work_dir}{basename}{ext_template.format(rh=radius_hit)}"))
    for path in sorted(leftovers):
        # Skip files promoted in step 2 (already moved on a real run; this keeps
        # --dry-run honest), the backups, and the summaries.
        if path in promoted:
            continue
        name = os.path.basename(path)
        if name.endswith(".orig.bak"):
            continue
        if any(name.endswith(s) for s in SUMMARY_SUFFIXES):
            continue
        _remove(path, dry_run)

    # 4. Warn if the promoted surface still needs a curvature graph.
    if not has_avv:
        print(f"  WARNING: iteration {step} was a lightweight (xcorr) iteration with no "
              f"curvature graph.")
        print(f"           The promoted surface has no .AVV_rh{radius_hit}.gt and is NOT "
              f"ready for")
        print(f"           downstream analysis. Re-run pycurv on it first:")
        print(f"             python run_pycurv.py config.yml {basename}.surface.vtp")

    return has_avv


@click.command()
@click.argument("configfile", type=click.Path(exists=True))
@click.argument("step", type=int)
@click.option("--component_name", "--component", "-c", "component_name", default=None,
              help="Optional: only accept surfaces for this component (e.g. OMM). Default: all components.")
@click.option("--tomogram", "-t", default=None,
              help="Only accept surfaces whose basename starts with this tomogram name. Default: all.")
@click.option("--dry-run", is_flag=True, default=False,
              help="Print the actions without renaming or deleting anything.")
def accept_refinement_cli(configfile, step, component_name, tomogram, dry_run):
    """Accept refinement iteration STEP as the new working surface.

    CONFIGFILE: path to config.yml
    STEP: the refinement iteration number to accept (matches *_refined_iter{STEP}.*)

    Backs up the original surfaces to *.orig.bak, promotes the chosen iteration
    to be the main surface (keeping only the canonical .surface.vtp and
    .AVV_rh*.gt/.vtp/.csv), and removes all other refinement files and
    regenerable pycurv intermediates (.NVV_rh*.gt, .scaled_cleaned.*,
    .lightweight.*, sampling CSVs, per-iteration plots). The per-surface
    refinement summaries (_refinement_stats.csv, _refinement_convergence.png,
    _profile_evolution.png) are kept.

    If the chosen iteration is a lightweight (xcorr) iteration with no curvature
    graph, the command warns that pycurv must be re-run on the promoted surface
    before any downstream analysis.
    """
    with open(configfile) as f:
        config = yaml.safe_load(f)

    work_dir = config.get("work_dir", config.get("seg_dir", "./"))
    if not work_dir.endswith("/"):
        work_dir += "/"
    radius_hit = config.get("curvature_measurements", {}).get("radius_hit", 9)

    print(f"Work directory: {work_dir}")
    print(f"Radius hit: {radius_hit}")
    print(f"Accepting iteration: {step}")
    if dry_run:
        print("DRY RUN - no files will be changed")
    print("-" * 60)

    # Discover every refined surface and the iterations it actually has (a surface
    # that converged early stops before `step`).
    iters_by_surface = _iterations_by_surface(work_dir)
    basenames = sorted(iters_by_surface)

    # Apply optional component / tomogram filters on the basename.
    if tomogram:
        basenames = [b for b in basenames if b.startswith(tomogram)]
    if component_name:
        basenames = [b for b in basenames if b.endswith(f"_{component_name}")]

    if not basenames:
        print(f"No refined surfaces found in {work_dir}"
              + (" matching the given filters." if (component_name or tomogram) else "."))
        return

    accepted = 0
    needs_pycurv = []
    for basename in basenames:
        available = iters_by_surface[basename]
        # Use the requested iteration, or the last available one if refinement
        # converged before it (largest available iteration <= step).
        usable = [i for i in available if i <= step] or available
        effective_step = usable[-1]
        if effective_step != step:
            print(f"  WARNING: {basename}: iteration {step} not available - refinement "
                  f"converged early at iteration {effective_step}; using iteration "
                  f"{effective_step}.")
        result = accept_one(work_dir, basename, effective_step, radius_hit, dry_run)
        if result is not None:
            accepted += 1
            if result is False:
                needs_pycurv.append(basename)
        print()

    print("-" * 60)
    verb = "Would accept" if dry_run else "Accepted"
    print(f"{verb} refinement for {accepted} surface(s) (requested iteration {step}; "
          "early-converged surfaces used their last iteration).")
    if accepted:
        print("Originals backed up with a .orig.bak suffix; only canonical surfaces "
              "and refinement summaries were kept.")
    if needs_pycurv:
        print()
        print("ACTION REQUIRED: the following surface(s) came from a lightweight (xcorr)")
        print("iteration and have no curvature graph. Re-run pycurv before any downstream")
        print("step (distances, density sampling, thickness):")
        for basename in needs_pycurv:
            print(f"  morphometrics pycurv {configfile} {basename}.surface.vtp")


if __name__ == "__main__":
    accept_refinement_cli()

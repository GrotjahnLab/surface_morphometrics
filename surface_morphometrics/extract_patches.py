#! /usr/bin/env python
"""Extract subsets of a membrane graph into their own .gt/.vtp/.csv files.

This generalizes GrotjahnLab/patch_analysis/extract_single_patch.py from
"split by patch_number" into an arbitrary graph subsetter with two modes:

1. Split by a label column (one output file per unique id):
     extract_patches.py config.yml --graph TS1_IMM..._patches.gt --by patch_number
   Produces `..._patch_number<id>.{gt,vtp,csv}` for each id > 0.

2. Filter by a numeric vertex property range (one output file):
     extract_patches.py config.yml --graph TS1_IMM.AVV_rh9.gt --property OMM_dist --max 30
   Keeps triangles where the property is within [min, max] (either bound
   optional) and writes a single `..._OMM_dist_filtered.{gt,vtp,csv}`. Useful
   for, e.g., grabbing only the ER within 30 nm of the OMM.

Any per-triangle vertex property can be used (curvature, thickness, the
`{label}_dist` distance properties, patch ids, etc.).
"""

__author__ = "Benjamin Barad"
__email__ = "benjamin.barad@gmail.com"
__license__ = "GPLv3"

import os
from glob import glob

import click
import numpy as np
import yaml


def mask_from_property(values, vmin=None, vmax=None):
    """Boolean mask of triangles whose property is within [vmin, vmax].

    NaN values are always excluded. Either bound may be None (open-ended).
    """
    values = np.asarray(values, dtype=float)
    mask = ~np.isnan(values)
    if vmin is not None:
        mask &= values >= vmin
    if vmax is not None:
        mask &= values <= vmax
    return mask


def _save_subset(source_graph_file, mask, out_base):
    """Filter a freshly loaded graph by a boolean mask and save gt/vtp/csv.

    Reloads the graph so repeated calls (e.g. per patch) start from a clean,
    unfiltered copy — matching extract_single_patch.py's per-id reload.
    """
    from pycurv import TriangleGraph, io
    from graph_tool import load_graph, GraphView
    from .intradistance_verticality import export_csv

    np.bool = bool  # pycurv/graph-tool compatibility shim (deprecated numpy alias)

    tg = TriangleGraph()
    tg.graph = load_graph(source_graph_file)
    vfilt = tg.graph.new_vertex_property("bool")
    vfilt.a = mask.astype(bool)

    filtered = GraphView(tg.graph, vfilt=vfilt)
    filtered.purge_vertices()
    filtered.save(f"{out_base}.gt")

    # Reload the saved subset to emit a clean vtp + csv.
    tg2 = TriangleGraph()
    tg2.graph = load_graph(f"{out_base}.gt")
    surf = tg2.graph_to_triangle_poly()
    io.save_vtp(surf, f"{out_base}.vtp")
    export_csv(tg2, f"{out_base}.csv")


def extract_by_label(graph_file, label_col, output_dir):
    """Split a graph into one subset per unique value (> 0) of label_col."""
    from pycurv import TriangleGraph
    from graph_tool import load_graph

    np.bool = bool

    tg = TriangleGraph()
    tg.graph = load_graph(graph_file)
    if label_col not in tg.graph.vp:
        print(f"  property '{label_col}' not found in {os.path.basename(graph_file)}; skipping")
        return
    labels = tg.graph.vp[label_col].get_array()
    ids = sorted(int(v) for v in np.unique(labels) if v != 0)
    base = os.path.splitext(os.path.basename(graph_file))[0]
    print(f"  splitting on '{label_col}': {len(ids)} subset(s)")
    for rid in ids:
        mask = labels == rid
        out_base = os.path.join(output_dir, f"{base}_{label_col}{rid}")
        _save_subset(graph_file, mask, out_base)
        print(f"    {os.path.basename(out_base)} ({int(mask.sum())} triangles)")


def extract_by_property(graph_file, prop, vmin, vmax, output_dir):
    """Extract the single subset of a graph where prop is within [vmin, vmax]."""
    from pycurv import TriangleGraph
    from graph_tool import load_graph

    np.bool = bool

    tg = TriangleGraph()
    tg.graph = load_graph(graph_file)
    if prop not in tg.graph.vp:
        print(f"  property '{prop}' not found in {os.path.basename(graph_file)}; skipping")
        return
    values = tg.graph.vp[prop].get_array()
    mask = mask_from_property(values, vmin, vmax)
    if not mask.any():
        print(f"  no triangles in range for '{prop}'; skipping")
        return
    base = os.path.splitext(os.path.basename(graph_file))[0]
    out_base = os.path.join(output_dir, f"{base}_{prop}_filtered")
    _save_subset(graph_file, mask, out_base)
    print(f"  {os.path.basename(out_base)} ({int(mask.sum())} of {len(mask)} triangles)")


@click.command()
@click.argument("configfile", type=click.Path(exists=True))
@click.option("--graph", "graph_file", type=click.Path(exists=True), default=None,
              help="Single .gt graph to process (overrides batch mode).")
@click.option("--by", "label_col", default=None,
              help="Split mode: label column to split on (e.g. patch_number, "
                   "patch_random_number, component_number).")
@click.option("--property", "prop", default=None,
              help="Filter mode: vertex property to threshold (e.g. OMM_dist).")
@click.option("--min", "vmin", type=float, default=None, help="Filter mode: minimum value.")
@click.option("--max", "vmax", type=float, default=None, help="Filter mode: maximum value.")
@click.option("--pattern", default=None,
              help="Batch glob within work_dir (default: *_patches.gt for --by, "
                   "*.AVV_rh{radius_hit}.gt for --property).")
@click.option("--output-dir", "output_dir", default=None,
              help="Output directory (defaults to work_dir from config).")
def extract_patches_cli(configfile, graph_file, label_col, prop, vmin, vmax,
                        pattern, output_dir):
    """Extract graph subsets by label (--by) or by property range (--property).

    CONFIGFILE: path to config.yml.
    """
    if (label_col is None) == (prop is None):
        raise click.UsageError("Specify exactly one of --by <label> or --property <name>.")

    with open(configfile) as f:
        config = yaml.safe_load(f)
    work_dir = config.get("work_dir", config.get("seg_dir", "./"))
    if not work_dir.endswith("/"):
        work_dir += "/"
    if output_dir is None:
        output_dir = work_dir
    elif not output_dir.endswith("/"):
        output_dir += "/"
    os.makedirs(output_dir, exist_ok=True)
    radius_hit = config.get("curvature_measurements", {}).get("radius_hit", 9)

    def process(gf):
        print(f"Processing {os.path.basename(gf)}")
        if label_col is not None:
            extract_by_label(gf, label_col, output_dir)
        else:
            extract_by_property(gf, prop, vmin, vmax, output_dir)

    if graph_file is not None:
        process(graph_file)
        return

    if pattern is None:
        pattern = "*_patches.gt" if label_col is not None else f"*.AVV_rh{radius_hit}.gt"
    graphs = sorted(glob(work_dir + pattern))
    # Avoid re-processing our own extractor outputs.
    graphs = [g for g in graphs if "_filtered" not in g]
    if not graphs:
        print(f"No graphs matching {work_dir}{pattern}")
        return
    print(f"Found {len(graphs)} graph(s) to process")
    for gf in graphs:
        process(gf)


if __name__ == "__main__":
    extract_patches_cli()

#! /usr/bin/env python
"""Label connected components of membrane graphs for per-region statistics.

Each connected component of a triangle graph is given a unique integer id in a
`component_number` vertex property.  This parallels the `patch_number` property
written by generate_patches.py, so the same downstream per-region statistics
(average curvature, thickness, etc.) work on connected components or patches.

Components are numbered 1..N by descending size (largest first); 0 is reserved
for triangles in components smaller than `min_component_size` (so it parallels
the "no patch" sentinel used by generate_patches.py).

Two usage modes:
1. Batch over all AVV graphs in work_dir:
     label_connected_components.py config.yml
2. A single graph:
     label_connected_components.py config.yml --graph TS1_IMM.AVV_rh8.gt
"""

__author__ = "Benjamin Barad"
__email__ = "benjamin.barad@gmail.com"
__license__ = "GPLv3"

import os
from glob import glob

import click
import numpy as np
import yaml

from .config_utils import load_config


def label_component_numbers(component_index, min_size=0):
    """Convert raw component indices into 1-based ids ordered by size.

    Parameters
    ----------
    component_index : (T,) int array of raw component labels (e.g. from
        graph_tool.label_components), arbitrary numbering starting at 0.
    min_size : drop components with fewer than this many triangles (their
        triangles get id 0).

    Returns
    -------
    (T,) int array of component ids: 1..N by descending size, 0 = dropped.
    """
    component_index = np.asarray(component_index)
    labels, counts = np.unique(component_index, return_counts=True)
    # Order surviving components by descending size.
    keep = [(lab, cnt) for lab, cnt in zip(labels, counts) if cnt >= min_size]
    keep.sort(key=lambda lc: (-lc[1], lc[0]))
    remap = {lab: i + 1 for i, (lab, _) in enumerate(keep)}
    out = np.zeros_like(component_index, dtype=np.int64)
    for lab, new_id in remap.items():
        out[component_index == lab] = new_id
    return out


def label_components_single(graph_file, min_component_size, output_dir):
    """Label connected components for a single graph and save .gt/.vtp/.csv."""
    from pycurv import TriangleGraph, io
    from graph_tool import load_graph
    from graph_tool.topology import label_components
    from .intradistance_verticality import export_csv

    np.bool = bool  # pycurv/graph-tool compatibility shim (deprecated numpy alias)

    print(f"Processing graph: {graph_file}")
    tg = TriangleGraph()
    tg.graph = load_graph(graph_file)

    raw, hist = label_components(tg.graph, directed=False)
    component_index = raw.get_array()
    numbers = label_component_numbers(component_index, min_size=min_component_size)
    n_components = int(numbers.max())
    print(f"  {len(hist)} raw components -> {n_components} kept "
          f"(min_component_size={min_component_size})")

    vp = tg.graph.new_vertex_property("int")
    vp.a = numbers
    tg.graph.vp["component_number"] = vp

    base = os.path.splitext(os.path.basename(graph_file))[0]
    out_base = os.path.join(output_dir, f"{base}_components")
    tg.graph.save(f"{out_base}.gt")
    surf = tg.graph_to_triangle_poly()
    io.save_vtp(surf, f"{out_base}.vtp")
    export_csv(tg, f"{out_base}.csv")
    print(f"  Saved: {out_base}.gt / .vtp / .csv")


@click.command()
@click.argument("configfile", type=click.Path(exists=True))
@click.option("--graph", "graph_file", type=click.Path(exists=True), default=None,
              help="Single .gt graph to process (overrides batch mode).")
@click.option("--label", default=None,
              help="Restrict batch mode to this membrane label (e.g. IMM). Default: all.")
@click.option("--output-dir", "output_dir", default=None,
              help="Output directory (defaults to work_dir from config).")
@click.option("--min-size", "min_size", type=int, default=None,
              help="Minimum component size in triangles (defaults to "
                   "patch_analysis.min_component_size, else 0).")
def label_components_cli(configfile, graph_file, label, output_dir, min_size):
    """Label connected components of membrane graphs.

    CONFIGFILE: path to config.yml.
    """
    config = load_config(configfile, require=("work_dir",))

    pa_config = config.get("patch_analysis", {})
    work_dir = config.get("work_dir", config.get("seg_dir", "./"))
    if not work_dir.endswith("/"):
        work_dir += "/"
    if output_dir is None:
        output_dir = work_dir
    elif not output_dir.endswith("/"):
        output_dir += "/"
    os.makedirs(output_dir, exist_ok=True)

    radius_hit = config.get("curvature_measurements", {}).get("radius_hit", 9)
    if min_size is None:
        min_size = pa_config.get("min_component_size", 0)

    print("Connected component labeling settings:")
    print(f"  Work directory: {work_dir}")
    print(f"  Output directory: {output_dir}")
    print(f"  Minimum component size: {min_size}")

    if graph_file is not None:
        label_components_single(graph_file, min_size, output_dir)
        return

    # Batch mode: all AVV graphs (optionally restricted to a membrane label).
    if label:
        pattern = f"{work_dir}*_{label}.AVV_rh{radius_hit}.gt"
    else:
        pattern = f"{work_dir}*.AVV_rh{radius_hit}.gt"
    graphs = sorted(glob(pattern))
    # Don't re-process our own outputs.
    graphs = [g for g in graphs if not g.endswith("_components.gt")]
    if not graphs:
        print(f"No graphs matching {pattern}")
        return
    print(f"Found {len(graphs)} graph(s) to process")
    for gpath in graphs:
        label_components_single(gpath, min_size, output_dir)


if __name__ == "__main__":
    label_components_cli()

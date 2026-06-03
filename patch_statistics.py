#! /usr/bin/env python
"""Per-region statistics for membrane patches or connected components.

This reads the per-triangle CSVs written by `generate_patches.py`
(`*_patches.csv`) or `label_connected_components.py` (`*_components.csv`) and
computes area-weighted summary statistics for each region (each patch or
component), for any set of per-triangle properties.

Because both tools tag triangles with an integer region id (`patch_number`,
`patch_random_number`, or `component_number`), the same aggregation works for
all of them — protein patches, their random controls, and whole connected
components — producing one tidy CSV for downstream plotting/statistics.

Usage:
  patch_statistics.py config.yml                          # batch over work_dir CSVs
  patch_statistics.py config.yml --csv TS1_IMM.AVV_rh9_patches.csv
  patch_statistics.py config.yml --properties curvedness_VV,thickness --output out.csv

Region id 0 (triangles in no patch / dropped components) is always excluded.
"""

__author__ = "Benjamin Barad"
__email__ = "benjamin.barad@gmail.com"
__license__ = "GPLv3"

import os
from glob import glob

import click
import numpy as np
import pandas as pd
import yaml

# Candidate region-label columns and the region_type they map to. Whichever are
# present in a CSV get aggregated (unless --label-col forces one).
LABEL_COLUMNS = {
    "patch_number": "patch",
    "patch_random_number": "random",
    "component_number": "component",
}


def _weighted_mean(values, weights):
    """Area-weighted mean, ignoring NaN weights/values."""
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    total = weights.sum()
    if total == 0:
        return np.nan
    return float(np.sum(values * weights) / total)


def aggregate_regions(df, label_col, value_cols, area_col="area",
                      drop_zero_values=True):
    """Aggregate per-triangle data into per-region rows for one label column.

    Parameters
    ----------
    df : DataFrame of per-triangle properties (one row per triangle).
    label_col : column holding the integer region id (0 = excluded).
    value_cols : list of property columns to summarize.
    area_col : column with per-triangle area for area weighting. If missing,
        falls back to unweighted means (each triangle weight 1).
    drop_zero_values : if True, drop triangles where the value is 0 or NaN
        before summarizing that property (matches the original per-patch
        scripts, where 0 means "unmeasured").

    Returns
    -------
    DataFrame: one row per region with n_triangles, total_area, and
    `{prop}_mean` (area-weighted) + `{prop}_median` (area-weighted) per property.
    """
    from morphometrics_stats import weighted_median

    have_area = area_col in df.columns
    rows = []
    region_ids = sorted(int(r) for r in df[label_col].unique() if r != 0)
    for rid in region_ids:
        region = df[df[label_col] == rid]
        area_all = region[area_col].to_numpy(dtype=float) if have_area else \
            np.ones(len(region))
        row = {
            "region_id": rid,
            "n_triangles": len(region),
            "total_area": float(area_all.sum()),
        }
        for col in value_cols:
            if col not in region.columns:
                row[f"{col}_mean"] = np.nan
                row[f"{col}_median"] = np.nan
                continue
            vals = region[col].to_numpy(dtype=float)
            weights = area_all
            mask = ~np.isnan(vals)
            if drop_zero_values:
                mask &= vals != 0
            vals, weights = vals[mask], weights[mask]
            if len(vals) == 0:
                row[f"{col}_mean"] = np.nan
                row[f"{col}_median"] = np.nan
                continue
            row[f"{col}_mean"] = _weighted_mean(vals, weights)
            row[f"{col}_median"] = float(weighted_median(vals, weights))
        rows.append(row)
    return pd.DataFrame(rows)


def aggregate_csv(csv_file, value_cols, label_col=None, drop_zero_values=True):
    """Aggregate one CSV across all applicable label columns.

    Returns a tidy DataFrame tagged with `source` (the CSV base name) and
    `region_type` (patch / random / component).
    """
    df = pd.read_csv(csv_file)
    base = os.path.basename(csv_file)
    for suffix in ("_patches.csv", "_components.csv", ".csv"):
        if base.endswith(suffix):
            source = base[: -len(suffix)]
            break
    else:
        source = base

    if label_col is not None:
        label_cols = {label_col: LABEL_COLUMNS.get(label_col, label_col)}
    else:
        label_cols = {c: t for c, t in LABEL_COLUMNS.items() if c in df.columns}

    if not label_cols:
        print(f"  no region-label column in {base}; skipping")
        return pd.DataFrame()

    out = []
    for col, region_type in label_cols.items():
        part = aggregate_regions(df, col, value_cols,
                                 drop_zero_values=drop_zero_values)
        if part.empty:
            continue
        part.insert(0, "source", source)
        part.insert(1, "region_type", region_type)
        out.append(part)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


@click.command()
@click.argument("configfile", type=click.Path(exists=True))
@click.option("--csv", "csv_file", type=click.Path(exists=True), default=None,
              help="Single CSV to process (overrides batch mode).")
@click.option("--properties", default=None,
              help="Comma-separated per-triangle properties to summarize "
                   "(defaults to patch_analysis.statistics_properties).")
@click.option("--label-col", "label_col", default=None,
              help="Force a single region-label column (e.g. component_number). "
                   "Default: auto-detect patch_number/patch_random_number/component_number.")
@click.option("--pattern", default=None,
              help="Glob (within work_dir) of CSVs to aggregate in batch mode "
                   "(default: *_patches.csv).")
@click.option("--output", "output_csv", default=None,
              help="Output CSV path (default: work_dir/patch_statistics.csv).")
@click.option("--keep-zeros", is_flag=True, default=False,
              help="Keep triangles whose value is 0 (default drops them, like the "
                   "original per-patch scripts).")
def patch_statistics_cli(configfile, csv_file, properties, label_col, pattern,
                         output_csv, keep_zeros):
    """Compute per-region (patch/component) area-weighted statistics.

    CONFIGFILE: path to config.yml.
    """
    with open(configfile) as f:
        config = yaml.safe_load(f)

    pa_config = config.get("patch_analysis", {})
    work_dir = config.get("work_dir", config.get("seg_dir", "./"))
    if not work_dir.endswith("/"):
        work_dir += "/"

    if properties is not None:
        value_cols = [p.strip() for p in properties.split(",") if p.strip()]
    else:
        value_cols = pa_config.get("statistics_properties",
                                   ["curvedness_VV", "thickness"])

    drop_zero_values = not keep_zeros
    print("Patch statistics settings:")
    print(f"  Properties: {value_cols}")
    print(f"  Drop zero values: {drop_zero_values}")
    print(f"  Region label column: {label_col or 'auto-detect'}")

    if csv_file is not None:
        files = [csv_file]
    else:
        pattern = pattern or "*_patches.csv"
        files = sorted(glob(work_dir + pattern))
        if not files:
            print(f"No CSVs matching {work_dir}{pattern}")
            return
    print(f"  Processing {len(files)} CSV file(s)")

    results = []
    for f in files:
        print(f"Aggregating {os.path.basename(f)}")
        part = aggregate_csv(f, value_cols, label_col=label_col,
                             drop_zero_values=drop_zero_values)
        if not part.empty:
            results.append(part)

    if not results:
        print("No regions aggregated.")
        return
    combined = pd.concat(results, ignore_index=True)

    if output_csv is None:
        output_csv = os.path.join(work_dir, "patch_statistics.csv")
    combined.to_csv(output_csv, index=False)
    print(f"\nWrote {len(combined)} region rows to {output_csv}")
    # Brief summary per region type.
    for region_type, grp in combined.groupby("region_type"):
        print(f"  {region_type}: {len(grp)} regions across {grp['source'].nunique()} source(s)")


if __name__ == "__main__":
    patch_statistics_cli()

#! /usr/bin/env python
"""Maintainer tool: build the small distributable example dataset.

Crops a full tomogram + segmentation pair down to a small sub-volume that still
contains a useful stretch of two membranes (default IMM+OMM), downcasts to keep
the download small (tomogram -> float16, segmentation -> int8), writes a tailored
ready-to-run config.yml, and bundles everything into a .tar.gz to upload to
Zenodo. The segmentation compresses to almost nothing; the tomogram does not, so
its size is set by the crop window and dtype (see --box / --tomo-dtype).

This is a one-off developer utility (not part of the installed package). Run it
against the full data, then upload the resulting tarball and paste its URL +
sha256 into surface_morphometrics/fetch_example.py.

Example:
  python tools/make_example_dataset.py \
      --tomogram   ~/Downloads/test_refinement/tomograms/YTC042_2_lam10_ts_004.mrc \
      --segmentation ~/Downloads/test_refinement/segmentations/YTC042_2_lam10_ts_004_labels.mrc \
      --out ~/Downloads/surface_morphometrics_example.tar.gz
"""

import argparse
import hashlib
import os
import tarfile
import tempfile

import numpy as np
import mrcfile
import yaml

BUNDLE_DIRNAME = "surface_morphometrics_example"
TEMPLATE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "surface_morphometrics", "config_template.yml")


def parse_int_triple(text):
    parts = [int(x) for x in text.split(",")]
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("expected three comma-separated ints, e.g. 160,384,384")
    return np.array(parts)


def parse_labels(text):
    """'IMM=1,OMM=2' -> {'IMM': 1, 'OMM': 2}"""
    out = {}
    for item in text.split(","):
        name, _, value = item.partition("=")
        out[name.strip()] = int(value)
    return out


def crop_bounds(shape, center, box):
    half = box // 2
    mn = np.maximum(center - half, 0)
    mx = np.minimum(center + box - half, np.array(shape))  # exclusive
    return mn, mx


def downcast_tomogram(data, dtype):
    if dtype == "float16":
        return data.astype(np.float16)
    if dtype == "float32":
        return data.astype(np.float32)
    # int8 / int16: rescale the full data range into the integer range.
    lo, hi = float(data.min()), float(data.max())
    span = (hi - lo) or 1.0
    if dtype == "int8":
        scaled = (data - lo) / span * 255.0 - 128.0
        return np.round(scaled).astype(np.int8)
    if dtype == "int16":
        scaled = (data - lo) / span * 65535.0 - 32768.0
        return np.round(scaled).astype(np.int16)
    raise ValueError(f"unsupported tomo dtype {dtype}")


def write_mrc(path, data, voxel_size_angstrom):
    with mrcfile.new(path, overwrite=True) as m:
        m.set_data(data)
        m.voxel_size = voxel_size_angstrom


def build_config(labels, tomo_name, seg_name):
    """Tailor the packaged template to the example layout and label mapping."""
    with open(TEMPLATE) as f:
        cfg = yaml.safe_load(f)
    cfg["seg_dir"] = "./segmentations/"
    cfg["tomo_dir"] = "./tomograms/"
    cfg["work_dir"] = "./morphometrics/"
    cfg["exp_name"] = "example"
    cfg["segmentation_values"] = dict(labels)
    names = list(labels)
    cfg["distance_and_orientation_measurements"]["intra"] = names
    # Inter: first label vs the rest (if at least two labels).
    if len(names) >= 2:
        cfg["distance_and_orientation_measurements"]["inter"] = {names[0]: names[1:]}
    cfg["thickness_measurements"]["components"] = names
    header = (
        "# Surface Morphometrics example configuration (generated).\n"
        f"# Tomogram: tomograms/{tomo_name}\n"
        f"# Segmentation: segmentations/{seg_name}\n"
        "# Run the pipeline from this directory, e.g.:\n"
        "#   morphometrics make_meshes config.yml\n"
        "#   morphometrics pycurv config.yml <name>.surface.vtp\n"
    )
    return header + yaml.safe_dump(cfg, sort_keys=False)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tomogram", required=True)
    ap.add_argument("--segmentation", required=True)
    ap.add_argument("--out", required=True, help="Output .tar.gz path.")
    ap.add_argument("--box", type=parse_int_triple, default=parse_int_triple("160,384,384"),
                    help="Crop size z,y,x (default 160,384,384).")
    ap.add_argument("--center", type=parse_int_triple, default=None,
                    help="Crop center z,y,x (default: centroid of --center-label).")
    ap.add_argument("--center-label", type=int, default=2,
                    help="Label whose centroid is the crop center if --center is omitted (default 2 = OMM).")
    ap.add_argument("--keep-labels", type=parse_labels, default=parse_labels("IMM=1,OMM=2"),
                    help="name=value pairs to keep; all others zeroed (default IMM=1,OMM=2).")
    ap.add_argument("--tomo-dtype", choices=["float16", "int16", "int8", "float32"],
                    default="float16")
    args = ap.parse_args()

    print(f"Reading {args.segmentation}")
    with mrcfile.open(args.segmentation, permissive=True) as m:
        seg = m.data
        seg_vox = m.voxel_size.x
    print(f"Reading {args.tomogram}")
    with mrcfile.open(args.tomogram, permissive=True) as m:
        tomo = m.data
        tomo_vox = m.voxel_size.x
    if seg.shape != tomo.shape:
        raise SystemExit(f"shape mismatch: seg {seg.shape} vs tomo {tomo.shape}")

    keep_values = set(args.keep_labels.values())
    if args.center is not None:
        center = args.center
    else:
        center = np.argwhere(seg == args.center_label).mean(0).astype(int)
    mn, mx = crop_bounds(seg.shape, center, args.box)
    print(f"Crop center {center.tolist()}  ->  [{mn.tolist()} : {mx.tolist()}]  "
          f"dims {(mx - mn).tolist()}")

    sl = (slice(mn[0], mx[0]), slice(mn[1], mx[1]), slice(mn[2], mx[2]))
    tomo_crop = downcast_tomogram(tomo[sl], args.tomo_dtype)
    seg_crop = seg[sl].copy()
    seg_crop[~np.isin(seg_crop, list(keep_values))] = 0
    seg_crop = seg_crop.astype(np.int8)
    for name, value in args.keep_labels.items():
        print(f"  {name} ({value}): {(seg_crop == value).sum()} voxels")

    tomo_name = os.path.basename(args.tomogram)
    seg_name = os.path.basename(args.segmentation)

    with tempfile.TemporaryDirectory() as tmp:
        root = os.path.join(tmp, BUNDLE_DIRNAME)
        os.makedirs(os.path.join(root, "tomograms"))
        os.makedirs(os.path.join(root, "segmentations"))
        write_mrc(os.path.join(root, "tomograms", tomo_name), tomo_crop, tomo_vox)
        write_mrc(os.path.join(root, "segmentations", seg_name), seg_crop, seg_vox)
        with open(os.path.join(root, "config.yml"), "w") as f:
            f.write(build_config(args.keep_labels, tomo_name, seg_name))
        with open(os.path.join(root, "README.txt"), "w") as f:
            f.write(
                "Surface Morphometrics example dataset (cropped sub-volume).\n\n"
                f"Tomogram voxel size: {tomo_vox:.2f} A.  Labels kept: "
                f"{args.keep_labels}.\n\n"
                "Run the pipeline from this directory:\n"
                "  morphometrics make_meshes config.yml\n"
                "  morphometrics pycurv config.yml <name>.surface.vtp\n"
                "  morphometrics distances_orientations config.yml\n"
                "  morphometrics sample_density config.yml\n"
                "  morphometrics measure_thickness config.yml\n"
                "  # optional: morphometrics refine_mesh config.yml ; "
                "morphometrics accept_refinement config.yml <step>\n"
            )
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        with tarfile.open(args.out, "w:gz") as tar:
            tar.add(root, arcname=BUNDLE_DIRNAME)

    size_mb = os.path.getsize(args.out) / 1e6
    digest = sha256(args.out)
    print("\n" + "=" * 64)
    print(f"Wrote {args.out}")
    print(f"  size:   {size_mb:.1f} MB")
    print(f"  sha256: {digest}")
    print("Upload this to Zenodo, then set _ZENODO_URL and _ZENODO_SHA256 in")
    print("surface_morphometrics/fetch_example.py to the published file URL + this hash.")
    print("=" * 64)


if __name__ == "__main__":
    main()

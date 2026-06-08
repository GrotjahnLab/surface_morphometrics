# Changelog

All notable changes to the Surface Morphometrics toolkit are documented here.
This project loosely follows [Keep a Changelog](https://keepachangelog.com/) and
[Semantic Versioning](https://semver.org/).

## [2.0.0b1] — beta

This release turns the collection of scripts into an installable Python package
driven by a single `morphometrics` command, and adds protein-patch /
connected-component analysis and surface export. It is a **breaking change** to
how the toolkit is invoked.

### ⚠️ Breaking changes
- The pipeline is now run through one grouped command, e.g.
  `morphometrics make_meshes config.yml` instead of `python segmentation_to_meshes.py config.yml`.
  The old `python <script>.py …` calls still work for now via **deprecation shims**
  that warn and forward; they will be removed in a future release. See the
  [migration table](README.md#upgrading-from-older-versions).
- Config keys renamed: `data_dir` → `seg_dir`, and `max_triangles` →
  `simplify_max_triangles` (only used when `simplify: true`). `make_meshes` warns
  if it detects the old names.
- Thickness measurement now writes the per-triangle `thickness`/`offset`/
  `average_width` properties **back into the surface `.gt`/`.vtp`/`.csv` in place**
  rather than creating separate `*_refined.gt`/`.vtp` files.

### Added
- Installable `surface_morphometrics` package with a `morphometrics` console entry
  point (lazy subcommands, grouped/ordered `--help`, next-step hints, shell
  completion).
- `morphometrics new_config` — write a starter `config.yml` into the working dir.
- `morphometrics fetch_example` — download a small cropped tomogram+segmentation
  test set from Zenodo (or point at the full data on EMPIAR-12534).
- Density-guided mesh refinement workflow: `refine_mesh` + `accept_refinement`
  (commit a chosen iteration in place, with `.orig.bak` backups).
- Protein-patch and connected-component analysis: `generate_patches`
  (STAR-driven, with random controls and combined-STAR basename matching),
  `label_components`, `extract_patches`, and `patch_statistics`.
- `morphometrics export_obj` — export a quantified surface to a colormapped
  OBJ + MTL for Blender / MeshLab (with `--nan-color` for unmeasured faces).
- `pyproject.toml` packaging; `pip install -e .` is wired into the conda
  environment files; basic `tests/` for the pure-logic helpers.

### Changed
- Research/paper/one-off scripts moved to `old_scripts/`.
- pycurv no longer prints the "run in parallel / this may take an hour" warnings
  (it is now about as fast as mesh generation); the slow-step note moved to
  `refine_mesh`, which is the slowest, optional step.
- README reorganized (Installation / Quick start / Pipeline / Analysis &
  visualization / Reference / Upgrading) with a table of contents.

[2.0.0b1]: https://github.com/GrotjahnLab/surface_morphometrics/releases

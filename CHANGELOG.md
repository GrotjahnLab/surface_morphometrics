# Changelog

All notable changes to the Surface Morphometrics toolkit are documented here.
This project loosely follows [Keep a Changelog](https://keepachangelog.com/) and
[Semantic Versioning](https://semver.org/).

## [2.0.0b2] — beta

A quality and robustness release on top of 2.0.0b1, focused on the density-profile
bilayer fitting (mesh refinement and thickness measurement), a forgiving centralized
config system, and packaging fixes.

### Added
- Per-triangle `bilayer_resolution` score (0–1) written to the thickness CSV — a
  reliability flag distinguishing strictly-resolved leaflets from prior-recovered
  (slightly thin) measurements.
- `morphometrics new_config --simple` — a minimal, fully-runnable starter config
  (the full annotated template is still available via the default `--verbose`).
- Centralized config defaults: omitted parameter sections fall back to documented
  values (a single source of truth, kept in sync with the template by a test), so
  partial/stripped configs run.
- Per-command required-key validation with one clean error message instead of a
  `KeyError` traceback; directory paths no longer need trailing slashes.
- Distance/orientation measurements default to all-vs-all when `intra`/`inter` are
  omitted (announced; set them explicitly empty to opt out).
- Plugins can contribute `morphometrics` subcommands via a `morphometrics.commands`
  entry point.
- GitHub Actions test workflow and unit tests for the fitting and config helpers.

### Changed
- **Dual-Gaussian bilayer fitting overhaul** (mesh refinement and thickness): leaflet
  detection by height above the shared solvent base (fixes near-universal fallback to
  a single Gaussian on clean bilayers); a shared-width, centered fit over a symmetric
  window (removes a ~0.1 nm centering bias and the collapse failure mode); and a
  two-tier gate that recovers locally-merged bilayers from the global-average prior
  while still rejecting genuine single/skewed peaks. The whole-surface average fit and
  its `_fit.svg` plot use the same model.
- Thickness measurement now applies a real quality gate (resolution + R² + physical
  range) and reports `NaN` for unmeasurable triangles instead of a spurious value.
- `accept_refinement` falls back to a surface's last available iteration (with a
  warning) when it converged before the requested step.
- `mesh_refinement` defaults now match the documented template (the code's omitted-key
  defaults had drifted): `average_radius` 25, `average_radius_min` 12,
  `max_total_offset` 8, `laplacian_iterations` 5, `laplacian_lambda` 0.5,
  `convergence_threshold` 0.05, `xcorr_iterations` [1, 2, 3].
- `export_obj`: `--angstroms` renamed to `--scale_to_angstroms`, with config-aware
  unit conversion.

### Fixed
- pycurv fork deadlock on Linux (cap `OMP_NUM_THREADS=1` before importing pycurv).
- `export_obj` reads legacy / empty VTP surfaces instead of crashing.
- NumPy 2.x compatibility (replace the removed `np.trapz`).
- Deprecated `matplotlib.cm.get_cmap` usage.

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

[2.0.0b2]: https://github.com/GrotjahnLab/surface_morphometrics/releases
[2.0.0b1]: https://github.com/GrotjahnLab/surface_morphometrics/releases

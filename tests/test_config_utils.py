"""Tests for config loading: dir normalization, required keys, and defaults."""
import os
import tempfile
import textwrap

import click
import pytest
import yaml

from surface_morphometrics import config_utils as cu

TEMPLATE = os.path.join(os.path.dirname(cu.__file__), "config_template.yml")


def _write(text):
    path = tempfile.mktemp(suffix=".yml")
    with open(path, "w") as handle:
        handle.write(textwrap.dedent(text))
    return path


def test_defaults_match_template():
    # DEFAULTS is the single source of truth and must stay in sync with the template.
    with open(TEMPLATE) as handle:
        tmpl = yaml.safe_load(handle)
    # intra/inter hold example membership in the template but default to empty.
    skip = {("distance_and_orientation_measurements", "intra"),
            ("distance_and_orientation_measurements", "inter")}
    for section, default in cu.DEFAULTS.items():
        if isinstance(default, dict):
            assert section in tmpl, f"{section} missing from template"
            for key, value in default.items():
                if (section, key) in skip:
                    continue
                assert tmpl[section].get(key) == value, f"{section}.{key}"
        else:
            assert tmpl.get(section) == default, section


def test_trailing_slashes_added():
    p = _write("seg_dir: ./s\ntomo_dir: ./t\nwork_dir: ./w\n")
    c = cu.load_config(p)
    assert c["seg_dir"] == "./s/" and c["tomo_dir"] == "./t/" and c["work_dir"] == "./w/"
    os.remove(p)


def test_minimal_config_is_filled_with_defaults():
    p = _write("seg_dir: ./s\nwork_dir: ./w\nsegmentation_values:\n  OMM: 1\n")
    c = cu.load_config(p, require=("seg_dir", "work_dir", "segmentation_values"))
    assert c["cores"] == 6
    assert c["curvature_measurements"]["radius_hit"] == 9
    assert c["mesh_refinement"]["average_radius"] == 25          # template value, not code's old 12
    assert c["mesh_refinement"]["convergence_threshold"] == 0.05
    assert c["surface_generation"]["isotropic_remesh"] is True
    # intra/inter are intentionally NOT defaulted here: omission lets the distance
    # command default them to all-vs-all (vs. an explicit empty meaning "skip").
    assert "intra" not in c["distance_and_orientation_measurements"]
    assert c["distance_and_orientation_measurements"]["mindist"] == 3   # params still filled
    os.remove(p)


def test_user_values_override_defaults():
    p = _write("work_dir: ./w\ncurvature_measurements:\n  radius_hit: 15\n")
    c = cu.load_config(p)
    assert c["curvature_measurements"]["radius_hit"] == 15        # user wins
    assert c["curvature_measurements"]["min_component"] == 30     # rest filled from defaults
    os.remove(p)


def test_density_sampling_legacy_fallback():
    # Old configs put sampling settings under thickness_measurements; keep honoring them.
    p = _write("work_dir: ./w\nthickness_measurements:\n  sample_spacing: 0.5\n")
    c = cu.load_config(p)
    assert c["density_sampling"]["sample_spacing"] == 0.5         # carried from legacy location
    assert c["density_sampling"]["scan_range"] == 10             # default fills the rest
    os.remove(p)


def test_module_defaults_not_mutated():
    p = _write("work_dir: ./w\n")
    c = cu.load_config(p)
    c["mesh_refinement"]["xcorr_iterations"].append(99)
    assert cu.DEFAULTS["mesh_refinement"]["xcorr_iterations"] == [1, 2, 3]
    os.remove(p)


def test_distance_targets_default_to_all_vs_all():
    labels = ["OMM", "IMM", "ER"]
    intra, inter, intra_def, inter_def = cu.resolve_distance_targets({}, labels)
    assert intra == labels and intra_def is True
    # upper triangle, each unordered pair once
    assert inter == {"OMM": ["IMM", "ER"], "IMM": ["ER"]} and inter_def is True


def test_distance_targets_explicit_empty_means_none():
    intra, inter, intra_def, inter_def = cu.resolve_distance_targets(
        {"intra": [], "inter": {}}, ["OMM", "IMM"])
    assert intra == [] and inter == {} and intra_def is False and inter_def is False


def test_distance_targets_partial_value_is_honored():
    intra, inter, intra_def, inter_def = cu.resolve_distance_targets(
        {"intra": ["OMM"]}, ["OMM", "IMM", "ER"])
    assert intra == ["OMM"] and intra_def is False          # subset kept as given
    assert inter == {"OMM": ["IMM", "ER"], "IMM": ["ER"]}   # inter omitted -> all pairs
    assert inter_def is True


def test_distance_targets_single_label_has_no_pairs():
    intra, inter, _, _ = cu.resolve_distance_targets({}, ["OMM"])
    assert intra == ["OMM"] and inter == {}


def test_missing_required_key_raises_clean_error():
    p = _write("seg_dir: ./s\n")
    with pytest.raises(click.ClickException) as exc:
        cu.load_config(p, require=("seg_dir", "work_dir", "tomo_dir"))
    msg = exc.value.format_message()
    assert "work_dir" in msg and "tomo_dir" in msg and "seg_dir" not in msg.split("\n", 1)[1]
    os.remove(p)


def test_empty_value_counts_as_missing():
    p = _write("work_dir: ''\n")
    with pytest.raises(click.ClickException):
        cu.load_config(p, require=("work_dir",))
    os.remove(p)

"""Tests for the `morphometrics validate` command."""
import numpy as np
import pytest

mrcfile = pytest.importorskip("mrcfile")
from click.testing import CliRunner

from surface_morphometrics.validate import validate_cli, _present_label_values


def _write_seg(path, values):
    """Write a tiny segmentation MRC containing exactly the given label values."""
    arr = np.zeros((4, 4, 4), dtype=np.int8)
    for i, v in enumerate(values):
        arr.flat[i] = v
    with mrcfile.new(str(path), overwrite=True) as mrc:
        mrc.set_data(arr)


def _config(path, seg_dir, segvals, tomo_dir=None):
    lines = [f"seg_dir: {seg_dir}/"]
    if tomo_dir is not None:
        lines.append(f"tomo_dir: {tomo_dir}/")
    lines.append("segmentation_values:")
    lines += [f"  {name}: {val}" for name, val in segvals.items()]
    path.write_text("\n".join(lines) + "\n")


def test_present_label_values(tmp_path):
    p = tmp_path / "s.mrc"
    _write_seg(p, [1, 2, 2, 3])
    assert _present_label_values(str(p)) == {1, 2, 3}


def test_validate_matches_tomogram_and_reports_labels(tmp_path):
    seg, tomo = tmp_path / "segs", tmp_path / "tomos"
    seg.mkdir(); tomo.mkdir()
    _write_seg(seg / "TS1_labels.mrc", [1, 2])
    _write_seg(tomo / "TS1.mrc", [0])          # tomogram content is irrelevant here
    cfg = tmp_path / "c.yml"
    _config(cfg, seg, {"IMM": 1, "OMM": 2}, tomo_dir=tomo)
    r = CliRunner().invoke(validate_cli, [str(cfg)])
    assert r.exit_code == 0
    assert "matching tomogram: TS1.mrc" in r.output      # prefix match
    assert "IMM(1), OMM(2)" in r.output
    assert "All checks passed." in r.output


def test_validate_warns_on_name_mismatch(tmp_path):
    seg, tomo = tmp_path / "segs", tmp_path / "tomos"
    seg.mkdir(); tomo.mkdir()
    _write_seg(seg / "AAA_labels.mrc", [1])
    _write_seg(tomo / "BBB.mrc", [0])          # name does not prefix-match the seg
    cfg = tmp_path / "c.yml"
    _config(cfg, seg, {"IMM": 1}, tomo_dir=tomo)
    r = CliRunner().invoke(validate_cli, [str(cfg)])
    assert r.exit_code == 0                     # warnings do not fail
    assert "no matching tomogram" in r.output


def test_validate_warns_on_orphan_tomogram(tmp_path):
    seg, tomo = tmp_path / "segs", tmp_path / "tomos"
    seg.mkdir(); tomo.mkdir()
    _write_seg(seg / "TS1_labels.mrc", [1])
    _write_seg(tomo / "TS1.mrc", [0])         # matched by the segmentation
    _write_seg(tomo / "EXTRA.mrc", [0])       # orphan: no segmentation starts with it
    cfg = tmp_path / "c.yml"
    _config(cfg, seg, {"IMM": 1}, tomo_dir=tomo)
    r = CliRunner().invoke(validate_cli, [str(cfg)])
    assert r.exit_code == 0
    assert "EXTRA.mrc matches no segmentation" in r.output
    assert "matching tomogram: TS1.mrc" in r.output     # the matched one still reported


def test_validate_warns_on_unmapped_label_value(tmp_path):
    seg = tmp_path / "segs"
    seg.mkdir()
    _write_seg(seg / "TS1_labels.mrc", [1, 2])   # value 2 present but not configured
    cfg = tmp_path / "c.yml"
    _config(cfg, seg, {"IMM": 1})
    r = CliRunner().invoke(validate_cli, [str(cfg)])
    assert r.exit_code == 0
    assert "present but not in segmentation_values" in r.output


def test_validate_errors_on_missing_seg_dir(tmp_path):
    cfg = tmp_path / "c.yml"
    cfg.write_text("seg_dir: /no/such/dir/\nsegmentation_values:\n  IMM: 1\n")
    r = CliRunner().invoke(validate_cli, [str(cfg)])
    assert r.exit_code == 1
    assert "seg_dir does not exist" in r.output

"""Tests for the `morphometrics status` command (derived from work_dir artifacts)."""
from click.testing import CliRunner

from surface_morphometrics.status import status_cli, _refine_status


def _touch(path):
    open(path, "w").close()


def _csv(path, columns):
    with open(path, "w") as handle:
        handle.write(",".join(columns) + "\n")


def _dataset(tmp_path):
    """A work_dir with IMM fully processed (+accepted iter 4) and OMM meshed only."""
    seg, work = tmp_path / "segs", tmp_path / "work"
    seg.mkdir(); work.mkdir()
    _touch(seg / "TS1_labels.mrc")
    base = "TS1_labels"
    # IMM: full pipeline
    _touch(work / f"{base}_IMM.surface.vtp")
    _touch(work / f"{base}_IMM.AVV_rh9.gt")
    _csv(work / f"{base}_IMM.AVV_rh9.csv",
         ["index", "kappa_1", "verticality", "self_dist_min",
          "OMM_dist", "OMM_orientation", "thickness", "bilayer_resolution"])
    (work / f"{base}_IMM.accepted_iter").write_text("4\n")
    # OMM: meshed only (no curvature yet)
    _touch(work / f"{base}_OMM.surface.vtp")
    cfg = tmp_path / "c.yml"
    cfg.write_text(
        f"seg_dir: {seg}/\nwork_dir: {work}/\n"
        "segmentation_values:\n  IMM: 1\n  OMM: 2\n  ER: 3\n"
        "curvature_measurements:\n  radius_hit: 9\n")
    return cfg


def test_status_reports_per_surface_state(tmp_path):
    r = CliRunner().invoke(status_cli, [str(_dataset(tmp_path))])
    assert r.exit_code == 0
    out = r.output
    # IMM: everything done, accepted iter shown, inter partner + thickness derived
    assert "accepted@iter4" in out
    assert "inter OMM" in out
    assert "thickness ✓" in out
    # OMM: meshed but no curvature -> no measurements
    assert "OMM" in out and "curv —" in out
    # configured-but-unmeshed label noted
    assert "no mesh for: ER" in out


def test_status_writes_output_file(tmp_path):
    out = tmp_path / "status.txt"
    r = CliRunner().invoke(status_cli, [str(_dataset(tmp_path)), "-o", str(out)])
    assert r.exit_code == 0 and out.exists()
    assert "accepted@iter4" in out.read_text()


def test_refine_status_run_but_not_accepted(tmp_path):
    work = tmp_path / "w"; work.mkdir()
    for n in (1, 2, 3):
        _touch(work / f"S_IMM_refined_iter{n}.surface.vtp")
    state = _refine_status(str(work) + "/", "S_IMM", 9)
    assert "not accepted" in state and "iter1-3" in state


def test_accept_refinement_writes_marker(tmp_path):
    from surface_morphometrics.accept_refinement import accept_one
    work = str(tmp_path) + "/"
    base = "TS1_labels_IMM"
    _touch(tmp_path / f"{base}_refined_iter2.surface.vtp")   # the iteration to accept
    accept_one(work, base, 2, radius_hit=9, dry_run=False)
    assert (tmp_path / f"{base}.accepted_iter").read_text().strip() == "2"


def test_accept_refinement_marker_dry_run(tmp_path):
    from surface_morphometrics.accept_refinement import accept_one
    work = str(tmp_path) + "/"
    base = "TS1_labels_IMM"
    _touch(tmp_path / f"{base}_refined_iter2.surface.vtp")
    accept_one(work, base, 2, radius_hit=9, dry_run=True)
    assert not (tmp_path / f"{base}.accepted_iter").exists()   # dry-run writes nothing

"""Tests for dependency-light CLI helpers: new_config and the fetch checksum."""
import hashlib
import os
import tempfile

import yaml
from click.testing import CliRunner

from surface_morphometrics import cli
from surface_morphometrics import fetch_example


def test_new_config_writes_valid_template():
    runner = CliRunner()
    with runner.isolated_filesystem():
        result = runner.invoke(cli.new_config, [])
        assert result.exit_code == 0
        assert os.path.exists("config.yml")
        cfg = yaml.safe_load(open("config.yml"))
        # a few keys the pipeline relies on
        for key in ("seg_dir", "tomo_dir", "work_dir", "segmentation_values",
                    "curvature_measurements"):
            assert key in cfg


def test_new_config_refuses_overwrite_without_force():
    runner = CliRunner()
    with runner.isolated_filesystem():
        assert runner.invoke(cli.new_config, []).exit_code == 0
        # second call without --force should fail; --force succeeds
        assert runner.invoke(cli.new_config, []).exit_code != 0
        assert runner.invoke(cli.new_config, ["--force"]).exit_code == 0


def test_lazy_cli_lists_all_commands_without_heavy_imports():
    runner = CliRunner()
    out = runner.invoke(cli.main, ["--help"]).output
    for name in ("make_meshes", "pycurv", "measure_thickness", "generate_patches",
                 "export_obj", "new_config", "fetch_example"):
        assert name in out


def test_fetch_example_sha256_matches_hashlib():
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "blob.bin")
        payload = b"surface morphometrics" * 1000
        with open(path, "wb") as f:
            f.write(payload)
        assert fetch_example._sha256(path) == hashlib.sha256(payload).hexdigest()

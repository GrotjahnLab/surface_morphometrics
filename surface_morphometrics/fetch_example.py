#! /usr/bin/env python
"""`morphometrics fetch_example` - download a dataset for testing the pipeline.

Two sources:
  zenodo (default): a small (~tens of MB) cropped tomogram + segmentation that
    exercises the whole pipeline, including density sampling / thickness /
    refinement. Extracts a ready-to-run `surface_morphometrics_example/` folder
    (with its own config.yml).
  empiar: the full-size source data on EMPIAR-12534 (hundreds of MB to GB per
    file) for users who want to run on the real tomograms.
"""

__author__ = "Benjamin Barad"
__email__ = "benjamin.barad@gmail.com"
__license__ = "GPLv3"

import hashlib
import os
import sys
import tarfile
import tempfile
import urllib.request

import click

# Filled in after the example bundle is uploaded to Zenodo (see
# tools/make_example_dataset.py, which prints the sha256 to paste here).
_ZENODO_URL = "https://zenodo.org/records/20547754/files/surface_morphometrics_example.tar.gz?download=1"
_ZENODO_SHA256 = "de35899667a913af249480454dccf6cf107ce5128b7ac205e3f09da3e5ab1a9b"
_BUNDLE_DIRNAME = "surface_morphometrics_example"

# EMPIAR full data (the source the example is cropped from).
_EMPIAR_ID = "12534"
_EMPIAR_BASE = f"https://ftp.ebi.ac.uk/empiar/world_availability/{_EMPIAR_ID}/data/"


def _download(url, dest_path):
    """Stream a URL to dest_path with a simple progress indicator."""
    print(f"Downloading {url}")

    def _hook(block_num, block_size, total_size):
        if total_size > 0:
            pct = min(100, block_num * block_size * 100 // total_size)
            print(f"\r  {pct:3d}%  ({total_size/1e6:.1f} MB)", end="", file=sys.stderr)

    urllib.request.urlretrieve(url, dest_path, reporthook=_hook)
    print("", file=sys.stderr)


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _fetch_zenodo(dest, force):
    if not _ZENODO_URL or not _ZENODO_SHA256:
        raise click.ClickException(
            "The example dataset has not been published to Zenodo yet "
            "(_ZENODO_URL/_ZENODO_SHA256 are unset in fetch_example.py). "
            "Use `--source empiar` for the full data in the meantime."
        )
    target = os.path.join(dest, _BUNDLE_DIRNAME)
    if os.path.exists(target) and not force:
        raise click.ClickException(
            f"{target} already exists; use --force to overwrite."
        )
    os.makedirs(dest, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        archive = os.path.join(tmp, "example.tar.gz")
        _download(_ZENODO_URL, archive)
        digest = _sha256(archive)
        if digest != _ZENODO_SHA256:
            raise click.ClickException(
                f"Checksum mismatch (expected {_ZENODO_SHA256}, got {digest}). "
                "The download may be corrupt or the URL out of date."
            )
        print("Checksum OK; extracting...")
        with tarfile.open(archive) as tar:
            tar.extractall(dest)
    click.echo(f"\nExample dataset ready in {target}")
    click.echo(f"  cd {target} && morphometrics make_meshes config.yml")


def _fetch_empiar(dest, empiar_file, force):
    if not empiar_file:
        click.echo(f"EMPIAR-{_EMPIAR_ID} full data is served from:")
        click.echo(f"  {_EMPIAR_BASE}")
        click.echo("Browse the entry to find the file(s) you need, then re-run with")
        click.echo("  morphometrics fetch_example --source empiar --empiar-file <relative/path.mrc>")
        click.echo("(Each tomogram is hundreds of MB to GB. For a quick test use the default Zenodo source.)")
        return
    os.makedirs(dest, exist_ok=True)
    out = os.path.join(dest, os.path.basename(empiar_file))
    if os.path.exists(out) and not force:
        raise click.ClickException(f"{out} already exists; use --force to overwrite.")
    _download(_EMPIAR_BASE + empiar_file.lstrip("/"), out)
    click.echo(f"\nDownloaded {out}")


@click.command(name="fetch_example")
@click.option("--source", type=click.Choice(["zenodo", "empiar"]), default="zenodo",
              show_default=True, help="Where to fetch data from.")
@click.option("--dest", default=".", show_default=True, type=click.Path(),
              help="Directory to download into.")
@click.option("--empiar-file", default=None,
              help="With --source empiar: relative path of a file within the EMPIAR entry to download.")
@click.option("--force", is_flag=True, default=False, help="Overwrite existing files.")
def fetch_example_cli(source, dest, empiar_file, force):
    """Download example data for testing the pipeline (small Zenodo set or full EMPIAR data)."""
    if source == "zenodo":
        _fetch_zenodo(dest, force)
    else:
        _fetch_empiar(dest, empiar_file, force)


if __name__ == "__main__":
    fetch_example_cli()

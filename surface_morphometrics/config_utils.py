#! /usr/bin/env python
"""Shared helpers for loading the pipeline config.yml."""

__author__ = "Benjamin Barad"
__email__ = "benjamin.barad@gmail.com"
__license__ = "GPLv3"

import click
import yaml

# Directory settings that the pipeline joins with filenames by plain string
# concatenation (e.g. work_dir + "*.gt"), so they must end in a path separator.
DIR_KEYS = ("seg_dir", "tomo_dir", "work_dir")

# Human-readable purpose of each setting, used in the missing-required-key error.
_KEY_HELP = {
    "seg_dir": "directory of segmentation MRC files",
    "work_dir": "directory for pipeline output files",
    "tomo_dir": "directory of tomogram MRC files (needed for density sampling, "
                "thickness, and mesh refinement)",
}


def normalize_dirs(config):
    """Append a trailing "/" to any set directory path that lacks one.

    Mutates and returns ``config``. Only ``seg_dir``/``tomo_dir``/``work_dir`` that
    are non-empty strings are touched; empty or unset values are left alone so each
    command can apply its own defaulting.
    """
    if isinstance(config, dict):
        for key in DIR_KEYS:
            value = config.get(key)
            if isinstance(value, str) and value and not value.endswith("/"):
                config[key] = value + "/"
    return config


def require_keys(config, required, configfile=""):
    """Raise a clean CLI error if any ``required`` top-level key is unset.

    A key counts as unset if it is missing or maps to an empty/None value. All
    missing keys are reported together so the user fixes the config in one pass.
    ``required`` is per-command: each command demands only the settings it uses
    (e.g. tomogram steps require ``tomo_dir``; analysis/export steps need only
    ``work_dir``).
    """
    missing = [key for key in required if not (isinstance(config, dict) and config.get(key))]
    if missing:
        target = configfile or "config.yml"
        details = "\n".join(
            f"  - {key}" + (f"  ({_KEY_HELP[key]})" if key in _KEY_HELP else "")
            for key in missing
        )
        raise click.ClickException(
            f"{target} is missing required setting(s):\n{details}\n"
            "Add them to your config.yml (start one with `morphometrics new_config`)."
        )


def load_config(configfile, require=()):
    """Load a pipeline config YAML, normalizing directory paths and validating.

    A config written without trailing slashes (e.g. ``work_dir: ./morphometrics``)
    behaves identically to one with them, since the pipeline builds paths by string
    concatenation. ``require`` is a tuple of top-level keys this command needs; if any
    is missing the user gets one clean error (see :func:`require_keys`) instead of a
    ``KeyError`` traceback partway through the run.
    """
    with open(configfile) as handle:
        config = yaml.safe_load(handle) or {}
    normalize_dirs(config)
    require_keys(config, require, configfile)
    return config

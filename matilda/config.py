"""
config.py
=========
Site configuration for the Matilda daemon: where the web plots go, where conda
and the pynika / pyirena environments live, which Tiled server to ask, and how
far back to look for scans.

The values live in one JSON file, by default ``~/.matilda/config.json``
(override the location with the ``MATILDA_CONFIG`` environment variable).

* No file: the built-in defaults (``DEFAULTS``, the APS 12-ID-E values) are
  used, so the beamline behaves exactly as before.
* ``ensure_config_file()`` writes the defaults to that file when it does not
  exist yet, so a user has a complete, valid file to edit.  The daemon calls
  it at start-up; ``python -m matilda.config`` does the same and prints the
  values in force.  Merely importing this module never writes anything.
* A key missing from the file takes its default; an unknown key, a value of
  the wrong type or a file that is not valid JSON is logged and the default
  is used — a typo must not stop live reduction.  A broken file is never
  overwritten.

Environment variables keep precedence where they existed before:
``MATILDA_IMAGE_PATH`` overrides ``image_path`` (the literal ``none``
disables the web plots).

Core module: standard library only, safe on the headless daemon path.
"""

from __future__ import annotations

import json
import logging
import os
import socket
import sys

CONFIG_ENV_VAR = "MATILDA_CONFIG"
DEFAULT_CONFIG_PATH = os.path.join(os.path.expanduser("~"), ".matilda", "config.json")

# Built-in defaults: the values used at APS 12-ID-E.  Keep them identical to
# what the daemon used before the config file existed.
DEFAULTS: dict = {
    "image_path": "/home/joule/WEBUSAXS/www_live/",
    "conda_executable": "/APSshare/miniconda/x86_64/bin/conda",
    "pynika_conda_env_path": "/home/beams/USAXS/.conda/envs/pynika",
    "pyirena_conda_env_path": "/home/beams/USAXS/.conda/envs/pyirena",
    "tiled_server": "auto",
    "tiled_port": 8000,
    "tiled_catalog": "usaxs_MongoDB",
    "tiled_timeout_s": 10,
    "days_to_look_back": 1,
    "days_to_look_back_blanks": 5,
    "days_to_look_back_tunes": 1,
    "images_in_graphs": 10,
    "tunes_to_show": 5,
}

# Written into the file next to the values (JSON has no comments).
HELP: dict = {
    "image_path": "Folder for the live JPEG web plots; null disables them. "
                  "The MATILDA_IMAGE_PATH environment variable overrides this.",
    "conda_executable": "Full path to the conda executable used to run pynika and pyirena.",
    "pynika_conda_env_path": "Full path to the pynika conda environment.",
    "pyirena_conda_env_path": "Full path to the pyirena conda environment.",
    "tiled_server": "Tiled server host name. 'auto' = localhost on usaxscontrol, "
                    "usaxscontrol.xray.aps.anl.gov elsewhere.",
    "tiled_port": "Tiled server port.",
    "tiled_catalog": "Tiled catalog holding the USAXS runs.",
    "tiled_timeout_s": "Timeout of one Tiled request, seconds.",
    "days_to_look_back": "How many days back the daemon searches for new scans.",
    "days_to_look_back_blanks": "How many days back the daemon searches for blanks.",
    "days_to_look_back_tunes": "How many days back the daemon searches for tune scans.",
    "images_in_graphs": "Number of most recent scans shown in each web plot.",
    "tunes_to_show": "Number of tune curves per tune plot (keep <= 10).",
}

# Keys whose value may also be null.
_NULLABLE = {"image_path"}

_ABOUT = ("Matilda site configuration. Edit the values; '_about' and '_help' are "
          "ignored. Delete the file to get the defaults back. See docs/configuration.md.")

_BEAMLINE_HOSTS = ("usaxscontrol", "usaxscontrol.xray.aps.anl.gov")
_BEAMLINE_TILED_HOST = "usaxscontrol.xray.aps.anl.gov"


def config_path() -> str:
    """Path of the site configuration file (MATILDA_CONFIG or ~/.matilda/config.json)."""
    return os.path.expanduser(os.environ.get(CONFIG_ENV_VAR) or DEFAULT_CONFIG_PATH)


def _valid(key, value) -> bool:
    if value is None:
        return key in _NULLABLE
    default = DEFAULTS[key]
    if isinstance(default, bool):
        return isinstance(value, bool)
    if isinstance(default, int):
        return isinstance(value, int) and not isinstance(value, bool)
    if isinstance(default, float):
        return isinstance(value, (int, float)) and not isinstance(value, bool)
    return isinstance(value, type(default))


def load_config(path: str | None = None) -> dict:
    """Return the site configuration: DEFAULTS updated with the file's values.

    Never raises and never writes: problems are logged and the affected keys
    (or, for an unreadable file, all keys) keep their defaults.
    """
    path = path or config_path()
    cfg = dict(DEFAULTS)
    if not os.path.isfile(path):
        return cfg
    try:
        with open(path, encoding="utf-8") as f:
            data = json.load(f)
    except (OSError, ValueError) as exc:
        logging.error(f"Matilda config {path} could not be read ({exc}); using defaults.")
        return cfg
    if not isinstance(data, dict):
        logging.error(f"Matilda config {path} is not a JSON object; using defaults.")
        return cfg
    for key, value in data.items():
        if key.startswith("_"):
            continue
        if key not in DEFAULTS:
            logging.warning(f"Matilda config {path}: unknown key {key!r} ignored.")
            continue
        if not _valid(key, value):
            logging.warning(f"Matilda config {path}: {key}={value!r} has the wrong type; "
                            f"using the default {DEFAULTS[key]!r}.")
            continue
        cfg[key] = value
    return cfg


def ensure_config_file(path: str | None = None) -> str | None:
    """Write the defaults to the config file if it does not exist yet.

    Returns the path written, or None when the file already existed or could
    not be written (read-only home, ...; logged, not raised).
    """
    path = path or config_path()
    if os.path.exists(path):
        return None
    content = {"_about": _ABOUT, **DEFAULTS, "_help": HELP}
    try:
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        tmp = path + ".tmp"
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(content, f, indent=2)
            f.write("\n")
        os.replace(tmp, path)
    except OSError as exc:
        logging.warning(f"Could not write the default Matilda config to {path} ({exc}); "
                        "built-in defaults are used.")
        return None
    logging.info(f"Wrote default Matilda config to {path}.")
    return path


def resolve_image_path(cfg: dict) -> str | None:
    """Web-plot folder: MATILDA_IMAGE_PATH env var, else the config value.

    'none' (any case) or null disables image saving (returns None).
    """
    value = os.environ.get("MATILDA_IMAGE_PATH", cfg.get("image_path"))
    if value is None or str(value).strip().lower() in ("", "none"):
        return None
    return value


def resolve_tiled_server(cfg: dict) -> str:
    """Tiled host: the configured name, or for 'auto' localhost when running on
    usaxscontrol itself (avoids proxy/firewall issues), the beamline host otherwise."""
    server = cfg.get("tiled_server") or "auto"
    if server.lower() != "auto":
        return server
    if socket.gethostname() in _BEAMLINE_HOSTS:
        return "localhost"
    return _BEAMLINE_TILED_HOST


# Loaded once at import (read-only); the daemon modules take their constants
# from here.
SITE_CONFIG: dict = load_config()


def main() -> int:
    """``python -m matilda.config``: create the file if missing, print the values."""
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    path = config_path()
    written = ensure_config_file(path)
    print(f"Matilda config file: {path}" + ("  (created with defaults)" if written else ""))
    for key, value in load_config(path).items():
        marker = "" if value == DEFAULTS[key] else "   (changed from default)"
        print(f"  {key:26s} = {value!r}{marker}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

"""Site configuration (matilda.config): defaults, file creation, validation."""

import json

import pytest

from matilda import config


def test_no_file_gives_defaults(tmp_path):
    cfg = config.load_config(str(tmp_path / "missing.json"))
    assert cfg == config.DEFAULTS
    assert not (tmp_path / "missing.json").exists()     # loading never writes


def test_ensure_writes_editable_defaults_once(tmp_path):
    path = tmp_path / ".matilda" / "config.json"
    assert config.ensure_config_file(str(path)) == str(path)
    data = json.loads(path.read_text())
    for key, value in config.DEFAULTS.items():
        assert data[key] == value
    assert set(data["_help"]) == set(config.DEFAULTS)
    assert config.load_config(str(path)) == config.DEFAULTS
    # an existing file is never touched
    path.write_text(json.dumps({"tiled_port": 9000}))
    assert config.ensure_config_file(str(path)) is None
    assert json.loads(path.read_text()) == {"tiled_port": 9000}


def test_user_values_merged_and_bad_ones_rejected(tmp_path, caplog):
    path = tmp_path / "config.json"
    path.write_text(json.dumps({
        "tiled_server": "myhost.example.org",
        "days_to_look_back": 3,
        "image_path": None,
        "tiled_port": "8000",          # wrong type -> default
        "days_to_look_back_blanks": True,   # bool is not an int here
        "conda_executable": None,      # not nullable -> default
        "no_such_key": 1,              # unknown -> ignored
        "_comment": "ignored",
    }))
    cfg = config.load_config(str(path))
    assert cfg["tiled_server"] == "myhost.example.org"
    assert cfg["days_to_look_back"] == 3
    assert cfg["image_path"] is None
    assert cfg["tiled_port"] == config.DEFAULTS["tiled_port"]
    assert cfg["days_to_look_back_blanks"] == config.DEFAULTS["days_to_look_back_blanks"]
    assert cfg["conda_executable"] == config.DEFAULTS["conda_executable"]
    assert "no_such_key" not in cfg
    assert "unknown key" in caplog.text


def test_broken_json_falls_back_and_is_kept(tmp_path):
    path = tmp_path / "config.json"
    path.write_text("{ not json")
    assert config.load_config(str(path)) == config.DEFAULTS
    assert config.ensure_config_file(str(path)) is None
    assert path.read_text() == "{ not json"


def test_env_var_selects_file(tmp_path, monkeypatch):
    monkeypatch.setenv("MATILDA_CONFIG", str(tmp_path / "x.json"))
    assert config.config_path() == str(tmp_path / "x.json")


@pytest.mark.parametrize("env,cfg_value,expected", [
    (None, "/web/", "/web/"),
    (None, None, None),
    ("none", "/web/", None),
    ("/other/", "/web/", "/other/"),
])
def test_image_path_env_precedence(monkeypatch, env, cfg_value, expected):
    if env is None:
        monkeypatch.delenv("MATILDA_IMAGE_PATH", raising=False)
    else:
        monkeypatch.setenv("MATILDA_IMAGE_PATH", env)
    assert config.resolve_image_path({"image_path": cfg_value}) == expected


def test_tiled_server_auto(monkeypatch):
    assert config.resolve_tiled_server({"tiled_server": "h.example"}) == "h.example"
    monkeypatch.setattr(config.socket, "gethostname", lambda: "usaxscontrol")
    assert config.resolve_tiled_server({"tiled_server": "auto"}) == "localhost"
    monkeypatch.setattr(config.socket, "gethostname", lambda: "laptop")
    assert config.resolve_tiled_server({"tiled_server": "auto"}) == \
        "usaxscontrol.xray.aps.anl.gov"


def test_daemon_constants_match_old_hardcoded_values(monkeypatch):
    """Without a config file the daemon must see exactly the old constants."""
    import matilda.matilda as m
    import matilda.readfromtiled as rt
    import matilda.convertTune as ct
    assert m.CONDA_EXECUTABLE == "/APSshare/miniconda/x86_64/bin/conda"
    assert m.PYNIKA_CONDA_ENV_PATH == "/home/beams/USAXS/.conda/envs/pynika"
    assert m.PYIRENA_CONDA_ENV_PATH == "/home/beams/USAXS/.conda/envs/pyirena"
    assert (m.NumberOfDaysToLookBack, m.NumberOfDaysToLookBackBlanks,
            m.NumberOfImagesInGraphs) == (1, 5, 10)
    assert (rt.port, rt.catalog, rt.TILED_TIMEOUT) == (8000, "usaxs_MongoDB", 10)
    assert (ct.NumberOfTunesToShow, ct.NumberOfDaysToLookBackTunes) == (5, 1)

# Site configuration

Matilda's daemon needs a few site-specific settings: where to write the live
web plots, where conda and the pynika / pyirena environments are, which Tiled
server and catalog to query, and how far back to look for scans. These live in
one JSON file:

```
~/.matilda/config.json
```

Set the `MATILDA_CONFIG` environment variable to use a different file
(for example a shared file for a service account).

The reduction GUI (`matilda-gui`) does not use any of these settings.

## How it behaves

| Situation | What Matilda does |
|---|---|
| No file | Uses the built-in defaults: the APS 12-ID-E values. The beamline behaves exactly as before the file existed. |
| Daemon start-up, no file | Creates `~/.matilda/` and writes the defaults to `config.json`, so you have a complete file to edit. If the home folder is read-only, it logs a warning and carries on with the defaults. |
| File exists | Reads it. An existing file is **never rewritten**. |
| A key is missing from the file | That key keeps its default. |
| Unknown key, or a value of the wrong type | Logged as a warning in `matilda.log`. That key keeps its default. |
| File is not valid JSON | Logged as an error. **All** keys keep their defaults, and the file is left as it is so you can fix it. |

The file is read once, when Matilda starts, so **restart the daemon after
editing it**.

To create the file and see the values in force without starting the daemon:

```bash
python -m matilda.config
```

This prints the file path, writes the defaults if the file was missing, and
marks every value that differs from the default.

To go back to the defaults, delete the file. The next daemon start (or
`python -m matilda.config`) writes a fresh one.

## Settings

| Key | Default | Meaning |
|---|---|---|
| `image_path` | `/home/joule/WEBUSAXS/www_live/` | Folder for the live JPEG web plots. `null` turns them off. |
| `conda_executable` | `/APSshare/miniconda/x86_64/bin/conda` | Conda used to run pynika and pyirena (`conda run`). |
| `pynika_conda_env_path` | `/home/beams/USAXS/.conda/envs/pynika` | pynika environment (auto-calibration). |
| `pyirena_conda_env_path` | `/home/beams/USAXS/.conda/envs/pyirena` | pyirena environment (auto-analysis, merging). |
| `tiled_server` | `auto` | Tiled host name. `auto` means `localhost` when running on usaxscontrol itself and `usaxscontrol.xray.aps.anl.gov` everywhere else. |
| `tiled_port` | `8000` | Tiled port. |
| `tiled_catalog` | `usaxs_MongoDB` | Tiled catalog that holds the USAXS runs. |
| `tiled_timeout_s` | `10` | Timeout for one Tiled request, in seconds. |
| `days_to_look_back` | `1` | How many days back to search for new scans. |
| `days_to_look_back_blanks` | `5` | How many days back to search for blanks. |
| `days_to_look_back_tunes` | `1` | How many days back to search for tune scans. |
| `images_in_graphs` | `10` | How many of the most recent scans each web plot shows. |
| `tunes_to_show` | `5` | Tune curves per tune plot. Keep this at 10 or fewer. |

The generated file also contains an `_about` note and a `_help` block that
describes each key. Matilda ignores every key that starts with `_`.

### Environment variables

Two variables that existed before the config file still work:

* `MATILDA_IMAGE_PATH` **overrides** `image_path`. The literal value `none`
  turns the web plots off.
* `MATILDA_LOG_DIR` sets where `matilda.log` goes. It is not in the config
  file, because the log has to exist before the config is read.

## Example: running off the beamline

```json
{
  "image_path": null,
  "tiled_server": "tiled.example.org",
  "tiled_catalog": "my_catalog",
  "conda_executable": "/opt/miniconda3/bin/conda",
  "pynika_conda_env_path": "/opt/miniconda3/envs/pynika",
  "pyirena_conda_env_path": "/opt/miniconda3/envs/pyirena"
}
```

Any key you leave out keeps its default.

## Not in this file

Per-folder switches stay next to the data, as before:

* `pyirena_config.json` turns on pyirena analysis for a folder.
* `merge_config.json` turns on USAXS+SAXS merging.

See [operations.md](operations.md).

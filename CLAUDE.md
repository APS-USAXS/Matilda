# CLAUDE.md

Orientation file for AI agents working in this repository. It tells you *where*
things are and *what rules apply*, not what every function does.

Matilda reduces USAXS / SAXS / WAXS data collected at APS beamline 12-ID-E:
two desktop GUIs for users, plus a headless daemon for beamline staff. The
daemon runs **in production on `usaxscontrol.xray.aps.anl.gov`** — a bug here
corrupts live experiment data during a user's beamtime. Treat changes to the
reduction path accordingly.

---

## 1. Commands

```bash
pip install -e ".[gui,dev]"   # dev install
pytest                        # full suite (tests/, test_*.py only)
ruff check matilda/           # lint
matilda                       # the polling daemon (beamline only)
matilda-gui                   # interactive data reduction
matilda-sample-plates         # sample-plate position editor
```

`tests/manual_test.py` is a hand-run integration script with hardcoded local
paths — it is deliberately excluded from collection. Don't wire it into CI.

---

## 2. Architecture

```
Tiled server (HTTP/REST)
        │  FindLastScanData() / FindLastBlankScan()
        ▼
matilda.py  ── 15-second polling loop
        ├── convertFlyscan.py   USAXS flyscan  → I(Q)
        ├── convertUSAXS.py     USAXS step-scan → I(Q)
        ├── convertSWAXS.py     SAXS/WAXS 2-D   → I(Q) via pyFAI
        ├── pynika   (subprocess, conda run)  triggered by AgBehenate/LaB6 names
        ├── pyirena  (subprocess, conda run)  fit_pyirena, merge_data
        ├── plotData.py         JPEG summary plots → web directory
        └── hdf5code.py         NXcanSAS read/write
```

| Path | Responsibility |
|---|---|
| `matilda/matilda.py` | Orchestrator; `main()` is the console-script entry point |
| `matilda/convert*.py` | Per-technique reduction to I(Q) |
| `matilda/desmearing*.py` | Slit desmearing (USAXS) |
| `matilda/hdf5code.py` | All NXcanSAS HDF5 read/write |
| `matilda/readfromtiled.py` | Tiled queries |
| `matilda/supportFunctions.py` | Shared helpers; re-exports HDF5 helpers (`F401` allowed) |
| `matilda/supportNikaFunctions.py` | Nika-derived 2-D routines |
| `matilda/plotData.py` | matplotlib JPEG plots (target for pyqtgraph replacement) |
| `matilda/gui/data_reduction/` | Interactive reduction GUI |
| `matilda/gui/sample_plate_setup.py` | Sample-plate editor |

Tunable module-level constants sit near the top of `matilda.py` (`imagePath`,
`NumberOfDaysToLookBack`, `CONDA_EXECUTABLE`, the conda env paths, and the
config filenames `pyirena_config.json` / `merge_config.json`). See
`docs/architecture.md` for the full table before changing any of them.

### Invariants — do not break these

1. **Deployment target is RHEL 8 with GLIBC 2.28.** `PySide6` is pinned
   `>=6.4,<6.8` because 6.8+ requires GLIBC 2.32. Do not raise that ceiling
   without checking the beamline machine first.
2. **PySide6, never PyQt6.** The two in one environment break Qt
   platform-plugin resolution (on macOS: "cocoa not found"). The `[gui]` extra
   installs PySide6 only.
3. **The daemon must stay headless.** `matplotlib` is a core dependency and is
   used from the daemon path; Qt and `pyepics` live in the `[gui]` extra so a
   server install can omit them.
4. pynika and pyirena are invoked as **subprocesses via `conda run`**, not
   imported. They have their own environments on purpose — do not "simplify"
   this into an import.
5. NXcanSAS output is consumed downstream by pyirena, DataReporter and Igor
   Irena. Changing field names or units is a breaking change for all three.

---

## 3. Conventions

**Scientific.** Q in Å⁻¹, intensity in cm⁻¹ (absolute) after calibration.
Single-letter physics names (`I`, `Q`, `l`) are idiomatic — `E741` is disabled
for that reason, don't rename them.

**Code.** Python ≥3.11. `ruff` with `select = ["E4","E7","E9","F","B006"]` —
correctness rules only; style rules are deliberately off for now.
`line-length = 200` is a **legacy accommodation**, annotated as such in
`pyproject.toml`. New code should still target 100 and the ceiling should come
down as legacy files get rewritten. Version comes from git tags via `hatch-vcs`.

**Testing.** New reduction math needs a test in `tests/`. `test_smoke_reduction.py`
is the end-to-end guard — keep it passing. Tests must not need Tiled or a display.

---

## 4. Where to look

| If you are… | Read |
|---|---|
| Getting oriented | `docs/architecture.md` (full module reference) |
| Working on desmearing | `docs/desmearing-methods.md`, `docs/desmearing-paper.md` |
| Working on the reduction GUI | `docs/matilda-gui.md` |
| Working on sample plates | `docs/sample-plate-setup.md` |
| Deploying or debugging the daemon | `docs/service.md`, `docs/operations.md`, `serv_matilda.sh`, `matilda_server.service` |
| Installing | `docs/installation.md` |
| Looking for planned work | `IMPROVEMENT_PLAN.md` |

`ControlData/`, `TestData/` and `IgorExample/` are fixtures and references, not
shipped code; `CodeFragments`, `IgorExample` and `docs` are excluded from ruff.

---

## 5. Maintaining this file

This is a map, not documentation. Update it when a module is added or removed
from §2, when an invariant changes, or when a command in §1 changes. Do not
duplicate `docs/architecture.md` here — link to it.

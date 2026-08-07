"""Guards for the two Qt rules in CLAUDE.md.

Neither test needs Qt installed, and neither needs a display:

1. Every GUI module imports Qt through ``matilda.gui._qt`` (invariant 2 —
   one import point, PySide6, no PyQt6 alongside it). Checked by scanning
   source, so it fails on a *new* direct import even where Qt is absent.
2. The daemon path stays headless (invariant 3) — no module on the reduction
   path imports a Qt binding, or anything under ``matilda.gui``.
"""

import ast
import os

import pytest

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_GUI_DIR = os.path.join(_REPO_ROOT, "matilda", "gui")
_SHIM = os.path.join(_GUI_DIR, "_qt.py")

_QT_BINDINGS = ("PySide6", "PyQt6", "PySide2", "PyQt5")


def _gui_sources():
    for dirpath, _dirnames, filenames in os.walk(_GUI_DIR):
        for name in sorted(filenames):
            if name.endswith(".py"):
                path = os.path.join(dirpath, name)
                if os.path.abspath(path) != os.path.abspath(_SHIM):
                    yield path


def _binding_imports(path):
    """Return the Qt binding modules imported directly by one source file."""
    with open(path, encoding="utf-8") as fh:
        tree = ast.parse(fh.read(), filename=path)

    hits = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            hits += [a.name for a in node.names
                     if a.name.split(".")[0] in _QT_BINDINGS]
        elif isinstance(node, ast.ImportFrom) and node.module:
            if node.module.split(".")[0] in _QT_BINDINGS:
                hits.append(node.module)
    return hits


@pytest.mark.parametrize("path", list(_gui_sources()),
                         ids=lambda p: os.path.relpath(p, _REPO_ROOT))
def test_gui_module_imports_qt_only_through_the_shim(path):
    hits = _binding_imports(path)
    assert not hits, (
        f"{os.path.relpath(path, _REPO_ROOT)} imports {', '.join(sorted(set(hits)))} "
        "directly. Import Qt names from matilda.gui._qt instead — see the "
        "'GUI framework' note in matilda/gui/__init__.py."
    )


def test_shim_is_the_one_place_that_names_a_binding():
    assert _binding_imports(_SHIM), "matilda/gui/_qt.py should import the Qt binding"


def _daemon_sources():
    """Top-level matilda/*.py — the daemon and reduction path, gui/ excluded."""
    pkg_dir = os.path.join(_REPO_ROOT, "matilda")
    for name in sorted(os.listdir(pkg_dir)):
        if name.endswith(".py"):
            yield os.path.join(pkg_dir, name)


@pytest.mark.parametrize("path", list(_daemon_sources()),
                         ids=lambda p: os.path.relpath(p, _REPO_ROOT))
def test_daemon_module_does_not_import_qt_or_the_gui_package(path):
    """A ``pip install matilda`` without the [gui] extra has no Qt at all.

    Checked on source rather than on ``sys.modules`` after an import: pyFAI
    opportunistically loads ``silx.gui`` (hence PySide6) from
    ``pyFAI.io.ponifile`` when a binding happens to be installed. That import
    is wrapped in its own try/except and degrades to ``GeometryModel = None``,
    so it does not make the daemon depend on Qt — but it does mean a runtime
    probe cannot tell a Matilda import from a third-party one.
    """
    hits = _binding_imports(path)
    assert not hits, (
        f"{os.path.relpath(path, _REPO_ROOT)} imports {', '.join(sorted(set(hits)))}. "
        "The daemon path must stay headless (CLAUDE.md invariant 3); Qt lives "
        "in the [gui] extra."
    )

    with open(path, encoding="utf-8") as fh:
        tree = ast.parse(fh.read(), filename=path)
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom) and node.module:
            assert not node.module.startswith("matilda.gui"), (
                f"{os.path.relpath(path, _REPO_ROOT)}:{node.lineno} imports "
                f"{node.module} — that pulls Qt into the daemon path."
            )
        elif isinstance(node, ast.Import):
            for alias in node.names:
                assert not alias.name.startswith("matilda.gui"), (
                    f"{os.path.relpath(path, _REPO_ROOT)}:{node.lineno} imports "
                    f"{alias.name} — that pulls Qt into the daemon path."
                )

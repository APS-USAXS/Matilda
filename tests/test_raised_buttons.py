"""
Contract for the vendored ``raised_buttons`` module (raised push buttons).

The file is a byte-identical copy of ``pyirena/gui/raised_buttons.py``; do not
edit it here. Refresh it with ``python <pyirena>/scripts/sync_raised_buttons.py
--write``. These checks need no display.
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

from matilda.gui import raised_buttons as rb

_HEX = re.compile(r"^#[0-9a-fA-F]{6}$")


def test_module_imports_no_gui_toolkit():
    tree = ast.parse(Path(rb.__file__).read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            names = ([a.name for a in node.names] if isinstance(node, ast.Import)
                     else [node.module or ""])
            for name in names:
                assert not re.match(r"(PySide|PyQt)", name), name


def test_every_role_has_face_and_text_colour():
    for name, (bg, fg, bold) in rb.DEFAULT_ROLES.items():
        assert _HEX.match(bg) and _HEX.match(fg), name
        assert isinstance(bold, bool), name


def test_install_is_idempotent():
    class FakeApp:
        qss = ""

        def styleSheet(self):
            return self.qss

        def setStyleSheet(self, s):
            self.qss = s

    app = FakeApp()
    rb.install(app)
    once = app.qss
    rb.install(app)
    assert app.qss == once and 'QPushButton[btnRole="stop"]' in once

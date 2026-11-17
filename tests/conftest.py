"""Shared pytest fixtures and environment setup for the Matilda test suite."""

import os
import sys

import pytest

# Make the repo importable when the package is not pip-installed.
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

# Never read (or create) the developer's ~/.matilda/config.json: point the
# site config at a path that does not exist, so every test sees the built-in
# defaults.  Must happen before any matilda module is imported.
import tempfile  # noqa: E402

os.environ["MATILDA_CONFIG"] = os.path.join(
    tempfile.gettempdir(), "matilda-tests-no-such-dir", "config.json")

# Headless matplotlib for any module that imports pyplot (plotData).
try:
    import matplotlib
    matplotlib.use("Agg")
except ImportError:
    pass

TESTDATA_DIR = os.path.join(REPO_ROOT, "TestData")


@pytest.fixture
def testdata_dir():
    """Path to the bundled TestData folder; skips the test if absent."""
    if not os.path.isdir(TESTDATA_DIR):
        pytest.skip("TestData folder not available")
    return TESTDATA_DIR

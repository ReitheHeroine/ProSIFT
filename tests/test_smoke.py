"""
Smoke tests for ProSIFT.

These exist so pytest has something to discover on first install. Replace
with real tests as modules are implemented. It is fine to keep a `test_smoke`
file around for "does the package import at all" checks, but the real
value comes from per-module tests and regression tests against reference data.
"""

import sys


def test_python_version():
    """Placeholder: confirm we are on a supported Python version."""
    assert sys.version_info >= (3, 10), (
        f"ProSIFT requires Python 3.10+, got {sys.version_info[:2]}"
    )


def test_imports_placeholder():
    """
    Placeholder for package import check.

    Once the prosift package exists, replace with:
        import prosift
        assert prosift.__version__
    """
    assert True

"""Sphinx configuration for the ocal documentation."""
from __future__ import annotations

import sys
from pathlib import Path

# Make the ``ocal`` package importable for autodoc (src layout).
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

# -- Project information -----------------------------------------------------
project = "ocal"
author = "Hiroyuki Matsui, Tomoharu Okada, Koki Ozawa"
copyright = "2026, Hiroyuki Matsui, Tomoharu Okada, Koki Ozawa"
release = "0.1.0"
version = "0.1.0"

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
]

templates_path = ["_templates"]
exclude_patterns: list[str] = ["_build"]
language = "en"

# -- autodoc / autosummary ---------------------------------------------------
autosummary_generate = True
autodoc_member_order = "bysource"
autodoc_typehints = "description"
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
}

# -- napoleon (NumPy-style docstrings) ---------------------------------------
napoleon_google_docstring = False
napoleon_numpy_docstring = True

# -- intersphinx -------------------------------------------------------------
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
}

# -- HTML output -------------------------------------------------------------
html_theme = "furo"
html_static_path = ["_static", "../assets"]

"""Sphinx configuration for SILEXlight documentation."""

# pylint: disable=redefined-builtin
import sys
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

project = 'SILEXlight'
copyright = '2026, SILEXlight contributors'
author = 'SILEXlight contributors'

version_file = ROOT / 'SILEXlight' / '__init__.py'
version_text = version_file.read_text(encoding='utf-8')
match = re.search(r"__version__\s*=\s*['\"]([^'\"]+)['\"]", version_text)
release = match.group(1) if match else '0.0.0'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
]

autosummary_generate = True
autodoc_member_order = 'bysource'
autodoc_typehints = 'description'

templates_path = ['_templates']
exclude_patterns = []

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_title = f'{project} {release}'


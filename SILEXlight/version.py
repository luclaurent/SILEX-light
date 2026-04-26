#!/usr/bin/env python
"""Version helper used by packaging scripts."""

import os


def get_version(rel_path):
    """Read ``__version__`` from the provided relative file path."""
    v = None
    vOk = False
    with open(os.path.join('./', rel_path), 'r') as f:
        for line in f:
            if line.lstrip().startswith('__version__'):
                v = line.split('=')[-1].strip().replace("'", "").replace('"', "")
                vOk = True
                break
    if not vOk:
        raise RuntimeError("Unable to find version string.")
    return v

if __name__ == '__main__':
    print(get_version('SILEXlight/__init__.py'))
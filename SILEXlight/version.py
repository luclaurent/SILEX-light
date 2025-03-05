#!/usr/bin/env python
import os
def get_version(rel_path):
    """ Get version from version.py."""
    v = None
    vOk = False
    with open(os.path.join('./',rel_path),'r') as f:
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
project('SILEXlight', 'c',
  version : run_command('SILEXlight/version.py', check: true).stdout().strip(),
  license: 'LGPL-3.0',
  meson_version: '>=0.64.0',
  default_options : ['warning_level=2', 'buildtype=release']
)

# Based on:
# - https://numpy.org/doc/stable/f2py/buildtools/meson.html
# - https://meson-python.readthedocs.io/en/stable/how-to-guides/first-project.html
# - https://github.com/jameskermode/f90wrap/blob/master/f90wrap/meson.build
# - https://github.com/nennigb/pypolsys/blob/meson-build/meson.build
# - https://github.com/aburrell/apexpy/blob/main/meson.build

add_languages('c')
add_languages('fortran')

# Python configuration
py = import('python').find_installation(pure: false)
message(py.full_path())
message(py.get_install_dir())

optf2py = []
flagFortranf2py = []
# manage debug
if get_option('buildtype') == 'debug'
message('Debug mode')
  add_project_arguments('-O0', 
      '-g3',
      '-fbacktrace',
      '-fcheck=all',
      '-ffpe-trap=zero,overflow,underflow',
      '-Wall',
      '-Wextra',
      '-Warray-temporaries',
      '-Wconversion',
      '-ffree-line-length-0', 
      language : 'fortran')
  optf2py += ['--debug-capi']
  # flagFortranf2py += ['-O0', 
  #     '-g3',
  #     '-fbacktrace',
  #     '-fcheck=all',
  #     '-ffpe-trap=zero,overflow,underflow',
  #     '-Wall',
  #     '-Wextra',
  #     '-Warray-temporaries',
  #     '-Wconversion',
  #     '-ffree-line-length-0']
endif

py_dep = py.dependency()

inc_root = include_directories('.')

# Get include from numpy and f2py
# based on https://numpy.org/doc/stable/f2py/buildtools/meson.html 
incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py)

# declare Fortran sources files (always given in relative path to SILEXlight folder)
src_files = {
  'silex_lib_tet4_fortran': ['silex_lib_tet4_fortran.f'],
  'silex_lib_tri3_fortran': ['silex_lib_tri3_fortran.f'],
}
#set folder path
src_path = 'SILEXlight'

# iterate over the sources
foreach name,sources : src_files
  message('Processing ' + name)

  # create outputs
  name_c_module = name + 'module.c'
  name_f2py_wrappers = []
  name_f2py_wrappers2 = []
  name_f2py_wrappers = name + '-f2pywrappers.f'
  # adapt the sources to the full path
  sources_full = []
  foreach src : sources
    sources_full += [src_path + '/' + src]
    message(src_path + '/' + src)
  endforeach

  # create a custom target for each source
  message('Creating custom target for ' + name)

  # build command
  commandf2py =[py, 
    '-m',
    'numpy.f2py', 
    # '-c',
    '@INPUT@', 
    '-m',
    name, 
    '--lower']
  if optf2py.length() > 0
    commandf2py += optf2py
  endif

  message('F2py command:')
  message( flagFortranf2py.length())
  message( optf2py.length())

  res_custom = custom_target(name,
    output: [name_c_module, 
            name_f2py_wrappers2,
            name_f2py_wrappers],
    input: sources_full,
    command : commandf2py
  )

  # list source files for Python's extension
  src_extension = sources_full
  src_extension += res_custom


  # Assumes that f2py was run before in "name" module
  py.extension_module(name, 
    src_extension,
    incdir_f2py / 'fortranobject.c',
    include_directories: inc_np,
    dependencies : [py_dep],
    install: true,
    # link_with: obj_dep,
    # link_language: 'fortran',
    # native: true,
    subdir: 'SILEXlight'  # Folder relative to site-packages to install to
    )

  message('Processing done: ' + name)
endforeach

# Declare the sources for Python
py.install_sources([
  'SILEXlight/__init__.py',
  'SILEXlight/silex_lib_gmsh.py',
  'SILEXlight/silex_lib_tet4_python.py',
  'SILEXlight/silex_lib_tri3_python.py',
  ],
  pure: false,
  preserve_path: false,
  subdir: 'SILEXlight'
  )

py.install_sources([
  'SILEXlight/tests/__init__.py',
  'SILEXlight/tests/test_silex_lib_gmsh.py',
  'SILEXlight/tests/test_silex_lib_tet4_fortran.py',
  'SILEXlight/tests/test_silex_lib_tet4_python.py',
  'SILEXlight/tests/test_silex_lib_tri3_fortran.py',
  'SILEXlight/tests/test_silex_lib_tri3_python.py',
  ],
  pure: false,
  preserve_path: false,
  subdir: 'SILEXlight/tests'
  )



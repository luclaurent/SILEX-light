# Create normal meson project.
project('SILEX-light', 
        'c',
        default_options : ['buildtype=release'])

# add Fortran language
add_languages('fortran')

# Import python3 meson module which is used to find the Python dependencies.
py3_mod = import('python')
# Locate the python executable.
py3 = py3_mod.find_installation('python',
    required: true,
)
py3_dep = py3.dependency()
#
message(py3.get_install_dir())
message(py3.path())
message(py3.get_path('purelib'))

# get includes of numpy
incdir_numpy = run_command(py3,
    ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
    check : true
    ).stdout().strip()

# get includes of f2py
incdir_f2py = run_command(py3,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

# add numpy's and f2py's includes to Meson includes
inc_np = include_directories(incdir_numpy,incdir_f2py)

srcpath = 'src/SILEXlight'
suf_module = 'module.c'
suf_f2pywrapper = '-f2pywrappers.f'

# build specific targets
module_name = 'silex_lib_tet4_fortran'
fortran_file = module_name + '.f'
silex_lib_tet4_source = custom_target(module_name + suf_module,
                            input : [ srcpath / fortran_file ], # .f so no F90 wrappers
                            output : [ module_name + suf_module, module_name + suf_f2pywrapper ],
                            command : [ py3, '-m', 'numpy.f2py', '@INPUT@',
                            '-m', module_name, '--lower']
                            )
# generate extension
py3.extension_module(
    module_name,
    srcpath / fortran_file,
    silex_lib_tet4_source,
    incdir_f2py / 'fortranobject.c',
    dependencies : [ py3_dep ],
    install: true,
    include_directories: inc_np,
    subdir: 'SILEXlight'
  )
# #
module_name = 'silex_lib_tri3_fortran'
fortran_file = module_name + '.f'
silex_lib_tri3_source = custom_target(module_name + suf_module,
                            input : [ srcpath / fortran_file ], # .f so no F90 wrappers
                            output : [ module_name + suf_module, module_name + suf_f2pywrapper ],
                            command : [ py3, '-m', 'numpy.f2py', '@INPUT@',
                            '-m', module_name, '--lower']
                            )
# generate extension
py3.extension_module(
    module_name,
    srcpath / fortran_file,
    silex_lib_tri3_source,
    incdir_f2py / 'fortranobject.c',
    dependencies : [ py3_dep ],
    install: true,
    include_directories: inc_np,
    subdir: 'SILEXlight'
  )

# #  Create the external C module using the python3 module helper function.
# py3.extension_module(
#     'silex_lib_tet4_fortran',
#     srcpath / 'silex_lib_tet4_fortran.f',
#     silex_lib_tet4_source,
#     incdir_f2py / 'fortranobject.c',
#     dependencies : [py3_dep],
#     install: true,
#     include_directories: inc_np,
#     subdir: 'SILEXlight'
#   )
# py3.extension_module(
#     'silex_lib_tri3_fortran',
#     'src/SILEXlight/silex_lib_tri3_fortran.f',
#     dependencies : [py3_dep],
#     install: true,
#     subdir: 'SILEXlight'
#   )

python_sources=[
    'src/SILEXlight/__init__.py',
    'src/SILEXlight/__version__.py',
    'src/SILEXlight/silex_lib_gmsh.py',
    'src/SILEXlight/silex_lib_tet4_python.py',
    'src/SILEXlight/silex_lib_tri3_python.py'
]
py3.install_sources(
    python_sources,
    pure: false,
    subdir: 'SILEXlight')

# Pathdir contains the dynamic library module.
pypathdir = meson.current_build_dir()


#   install_subdir(
#     'funmod',
#     install_dir: py3.get_install_dir()
# )


#   # Create a test script that runs a python script that uses
#   # the C external module.  Must populate PYTHONPATH to include
#   # a path that includes the mymath.so dynamic library.
# #   test('extmod',
# #     py3,
# #     args : files('mymath_test.py'),
# #     env : ['PYTHONPATH=' + pypathdir])

#   # Check we can apply a version constraint
#   # dependency('python3', version: '>=@0@'.format(py3_dep.version()))

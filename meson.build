# Create normal meson project.
project('SILEX-light', 
        'fortran',
        default_options : ['buildtype=release'])

# # Adding at project level causes many spurious -lgfortran flags.
# add_languages('fortran', native: false)
# ff = meson.get_compiler('fortran')
# if ff.has_argument('-Wno-conversion')
#   add_project_arguments('-Wno-conversion', language: 'fortran')
# endif



# Import python3 meson module which is used to find the Python dependencies.
py3_mod = import('python')
# Locate the python executable.
py3 = py3_mod.find_installation('python',
    required: true,
)
message(py3.get_install_dir())
message(py3.path())
message(py3.get_path('purelib'))
# # python_inst = py3.find_installation('python3')
# # py3insta = py3_mod.find_installation('python3')   
# # Create the Meson python3 dependency from the python3 module.
py3_dep = dependency('python3', required : false)

if py3_dep.found()
#  Create the external C module using the python3 module helper function.
py3.extension_module(
    'silex_lib_tet4_fortran',
    'src/SILEXlight/silex_lib_tet4_fortran.f',
    dependencies : [py3_dep],
    install: true,
    subdir: 'SILEXlight'
  )
py3.extension_module(
    'silex_lib_tri3_fortran',
    'src/SILEXlight/silex_lib_tri3_fortran.f',
    dependencies : [py3_dep],
    install: true,
    subdir: 'SILEXlight'
  )

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
else
  error('MESON_SKIP_TEST: Python3 libraries not found, skipping test.')
endif
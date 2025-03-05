import os
import sys
import sysconfig
import setuptools
import logging
from numpy.distutils.core import setup, Extension

logging.basicConfig(format='%(asctime)s - %(levelname)s -  %(message)s', level=logging.DEBUG)
logger = logging.getLogger()

# name of the package
__pkg_name__ = 'SILEXlight'
# default compilation argument
__default_compile_args__ = ['-O3',
                           '-ffast-math',
                           '-funroll-loops']  # "-fopenmp"

def convertToList(dataIn):
    if type(dataIn) is not list and dataIn:
        return [dataIn]
    else:
        return dataIn

def showMessageExtension(sourcefiles = [], 
                extensionName = None, 
                linkArgs = [],
                compileArgs = [],
                includeDirs = [],
                libraries = []):
    """ Show message about the extension building

    Args:
        sourcefiles (list): list of extension sources files (Fortran per instance). Defaults to [].
        extensionName (string): name of the extension inside the package. Defaults to [].
        linkArgs (list, optional): list of arguments for compile the sources with dedicated links. Defaults to [].
        compileArgs (list, optional): list of arguments for specific options for compilation (per instance -O3). Defaults to [].
        includeDirs (list, optional): list of directories that must be used for includes and .mod (Fortran's modules). Defaults to [].
        libraries (list, optional): list of specific libraries to be use to compile and build with the current setup.py. Defaults to [].
    """
    logger.info('Build extension {}'.format(__pkg_name__ + '.' + extensionName))
    txtMessage = ' >> source files: ' + ', '.join(sourcefiles)
    logger.debug(txtMessage)
    txtMessage = ' >> linking arguments: '
    if not linkArgs:
        txtMessage = txtMessage + 'None'
    else:
        txtMessage = txtMessage + ', '.join(linkArgs)
    logger.debug(txtMessage)
    txtMessage = ' >> compilation arguments: '
    if not compileArgs:
        txtMessage = txtMessage + 'None'
    else:
        txtMessage = txtMessage + ', '.join(compileArgs)
    logger.debug(txtMessage)
    txtMessage = ' >> include directories: '
    if not includeDirs:
        txtMessage = txtMessage + 'None'
    else:
        txtMessage = txtMessage + ', '.join(includeDirs)
    logger.debug(txtMessage)
    txtMessage = ' >> libraries for dependency: '
    if not libraries:
        txtMessage = txtMessage + 'None'
    else:
        txtMessage = txtMessage + ', '.join(libraries)
    logger.debug(txtMessage)
    pass

def showMessageLibrary(sourcefiles = [], libName = None):
    """ Show message about the extension building

    Args:
        sourcefiles (list): list of libraries sources files (Fortran per instance). Defaults to [].
        libName (string): name of the library. Defaults to None.
    """
    logger.info('Build library {}'.format(libName))
    txtMessage = ' >> source files: ' + ', '.join(sourcefiles)
    logger.debug(txtMessage)    
    pass

def lib_wrapper(sourcefiles = [], 
                libName = None):
    """Build library tuple

    Args:
        sourcefiles (list): list of libraries sources files (Fortran per instance). Defaults to [].
        libName (string): name of the library. Defaults to None.
    """
    # format source files
    sourcefiles = convertToList(sourcefiles)
    # show message
    showMessageLibrary(sourcefiles, libName)
    #
    return (libName, dict(sources=sourcefiles))

def ext_wrapper(sourcefiles = [], 
                extensionName = None, 
                linkArgs = [],
                compileArgs = [],
                includeDirs = [],
                libraries = []):
    """ Build extension for distutils

    Args:
        sourcefiles (list): list of extension sources files (Fortran per instance). Defaults to [].
        extensionName (string): name of the extension inside the package. Defaults to None.
        linkArgs (list, optional): list of arguments for compile the sources with dedicated links. Defaults to [].
        compileArgs (list, optional): list of arguments for specific options for compilation (per instance -O3). Defaults to [].
        includeDirs (list, optional): list of directories that must be used for includes and .mod (Fortran's modules). Defaults to [].
        libraries (list, optional): list of specific libraries to be use to compile and build with the current setup.py. Defaults to [].

    Returns:
        Extension: built extension
    """
    # format data
    sourcefiles = convertToList(sourcefiles)
    linkArgs = convertToList(linkArgs)
    compileArgs = convertToList(compileArgs)
    includeDirs = convertToList(includeDirs)
    libraries = convertToList(libraries)
    # show message
    showMessageExtension(sourcefiles, extensionName, linkArgs, compileArgs, includeDirs, libraries)
    
    ##
    return Extension(name=__pkg_name__ + '.' + extensionName,
                     sources=sourcefiles,
                     extra_compile_args=compileArgs,
                     extra_link_args=linkArgs,
                     include_dirs=includeDirs,
                     libraries=libraries)

def getTempDir():
    """Get the temporary build directory

    Returns:
        path: path of the temporary dir used to build fortran module (use to locate '.mod' files)
    """    
    plat_specifier = ".{}-{}".format(
        sysconfig.get_platform(), sys.version[0:3])
    return os.path.join('build','temp','.',plat_specifier)
    

# Load README
with open("README.md", "r") as fh:
    long_description = fh.read()
    
# Load version
def version():
    """ Get version from version.py."""
    v = None
    with open(os.path.join('./',__pkg_name__,'__init__.py')) as f:
        for line in f:
            if line.lstrip().startswith('__version__'):
                v = line.split('=')[-1].strip().replace("'", "").replace('"', "")
                break
        return v
    

#list of modules
modules = [
    ext_wrapper(sourcefiles=os.path.join(__pkg_name__,'silex_lib_tri3_fortran.f'),
                extensionName='silex_lib_tri3_fortran',
                linkArgs=[],
                compileArgs=__default_compile_args__,
                includeDirs=[],
                libraries=[]
    ),
    ext_wrapper(sourcefiles=os.path.join(__pkg_name__,'silex_lib_tet4_fortran.f'),
                extensionName='silex_lib_tet4_fortran',
                linkArgs=[],
                compileArgs=__default_compile_args__,
                includeDirs=[],
                libraries=[]
    )
    ]

# Fill setuptools
setup(
    name="SILEXlight",
    version=version(),
    long_description=long_description,
    long_description_content_type="text/markdown",
    # Use find_packages() to automatically discover all packages and subpackages
    packages=setuptools.find_packages(),
    include_package_data=True,
    # build f90 module
    ext_modules=modules
)
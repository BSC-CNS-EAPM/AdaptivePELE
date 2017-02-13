import numpy
from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension


ext_modules = [
        Extension("src.atomset.atomset", ["src/atomset/atomset.pyx"], include_dirs = ["src", "src/atomset"]),
        Extension("src.atomset.SymmetryContactMapEvaluator", ["src/atomset/SymmetryContactMapEvaluator.pyx"], include_dirs = ["src","src/atomset"]),
        Extension("src.atomset.RMSDCalculator", ["src/atomset/RMSDCalculator.pyx"], include_dirs = ["src", "src/atomset"])
            ]
setup(
    package_data={ "src/atomset": ['*.pxd'] },
    ext_modules = cythonize(ext_modules),  # accepts a glob pattern
    include_dirs=[numpy.get_include()]
)
# Run the following line to compile atomset package
# python setup.py build_ext --inplace

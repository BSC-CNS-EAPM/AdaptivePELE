import numpy
from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension


ext_modules = [
        Extension("atomset/atomset", ["atomset/atomset.pyx"], include_dirs = [".", "atomset"]),
        Extension("atomset/SymmetryContactMapEvaluator", ["atomset/SymmetryContactMapEvaluator.pyx"], include_dirs = [".","atomset"]),
        Extension("atomset/RMSDCalculator", ["atomset/RMSDCalculator.pyx"], include_dirs = [".", "atomset"])
            ]
setup(
    ext_modules = cythonize(ext_modules),  # accepts a glob pattern
    include_dirs=[numpy.get_include()]

)
# Run the following line to compile atomset package
# python setup.py build_ext --inplace


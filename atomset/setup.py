import numpy
from distutils.core import setup
from Cython.Build import cythonize
from distutils.extension import Extension


ext_modules = [
        Extension("atomset", ["atomset.pyx"]),
        Extension("SymmetryContactMapEvaluator", ["SymmetryContactMapEvaluator.pyx"]),
        Extension("RMSDCalculator", ["RMSDCalculator.pyx"])
            ]
setup(
    ext_modules = cythonize(ext_modules),  # accepts a glob pattern
    include_dirs=[numpy.get_include()]

)

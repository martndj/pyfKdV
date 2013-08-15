from distutils.core import setup 
from distutils.extension import Extension 
from Cython.Distutils import build_ext 
import numpy

npy_include_dir = numpy.get_include()

kdvLib=["kdvLib.so",]

ext_modules = [Extension("fKdV", ["fKdV.pyx"], 
                         include_dirs = [npy_include_dir],
                         extra_objects=kdvLib,
                         libraries=["fftw3",]
                         )]  

setup(name = 'fKdV',
      cmdclass = {'build_ext': build_ext},
      ext_modules = ext_modules)


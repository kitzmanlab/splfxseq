from setuptools import setup
#from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

import pysam
import numpy
import glob
import os.path as op

setup(
    name="splfxseq",

    packages=['spliceprocess','splanl'],

    include_dirs = [numpy.get_include()]+pysam.get_include(),

    entry_points = {
        'console_scripts': [ 'mpsa_cml = mpsa_pipe.mpsa_cml:main',
                             'gnomad_cml = mpsa_pipe.gnomad_cml:main',
                             'clinvar_cml = mpsa_pipe.clinvar_cml:main',
                             'splai_cml = mpsa_pipe.splai_cml:main',
                             'splai_dnv_cml = mpsa_pipe.splai_dnv_cml:main',
                           ]
    }
)

# ============================================================
# Copyright (c) 2012, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by Joel Bernier <bernier2@llnl.gov> and others.
# LLNL-CODE-529294.
# All rights reserved.
#
# This file is part of HEXRD. For details on dowloading the source,
# see the file COPYING.
#
# Please also see the file LICENSE.
#
# This program is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License (as published by the Free Software
# Foundation) version 2.1 dated February 1999.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this program (see file LICENSE); if not, write to
# the Free Software Foundation, Inc., 59 Temple Place, Suite 330,
# Boston, MA 02111-1307 USA or visit <http://www.gnu.org/licenses/>.
# ============================================================

import os
import sys

import numpy
np_include_dir = os.path.join(numpy.get_include(), 'numpy')

try:
    from setuptools import setup, Extension
    from setuptools import Command
    using_setuptools = True
except ImportError:
    if sys.platform.startswith('win'):
        # Require setuptools on Windows for installing the hexrd.exe entry point
        raise
    from distutils.core import setup, Extension
    from distutils.cmd import Command
    using_setuptools = False


import versioneer


versioneer.VCS = 'git'
versioneer.versionfile_source = 'hexrd/_version.py'
versioneer.versionfile_build = 'hexrd/_version.py'
versioneer.tag_prefix = 'v' # tags are like v1.2.0
versioneer.parentdir_prefix = 'hexrd-' # dirname like 'myproject-1.2.0'


class test(Command):

    """Run the test suite."""

    description = "Run the test suite"

    user_options = [('verbosity', 'v', 'set test report verbosity')]

    def initialize_options(self):
        self.verbosity = 0

    def finalize_options(self):
        try:
            self.verbosity = int(self.verbosity)
        except ValueError:
            raise ValueError('Verbosity must be an integer.')

    def run(self):
        import unittest
        suite = unittest.TestLoader().discover('hexrd')
        unittest.TextTestRunner(verbosity=self.verbosity+1).run(suite)


# for SgLite
srclist = [
    'sgglobal.c', 'sgcb.c', 'sgcharmx.c', 'sgfile.c', 'sggen.c', 'sghall.c',
    'sghkl.c', 'sgltr.c', 'sgmath.c', 'sgmetric.c', 'sgnorm.c', 'sgprop.c',
    'sgss.c', 'sgstr.c', 'sgsymbols.c', 'sgtidy.c', 'sgtype.c', 'sgutil.c',
    'runtests.c', 'sglitemodule.c'
    ]
srclist = [os.path.join('hexrd/sglite', f) for f in srclist]
sglite_mod = Extension(
    'hexrd.xrd.sglite',
    sources=srclist,
    define_macros=[('PythonTypes', 1)]
    )


# for transforms
srclist = ['transforms_CAPI.c', 'transforms_CFUNC.c']
srclist = [os.path.join('hexrd/transforms', f) for f in srclist]
transforms_mod = Extension(
    'hexrd.xrd._transforms_CAPI',
    sources=srclist,
    include_dirs=[np_include_dir]
    )


# all modules
ext_modules = [sglite_mod, transforms_mod]


packages = []
for dirpath, dirnames, filenames in os.walk('hexrd'):
    if '__init__.py' in filenames:
        packages.append('.'.join(dirpath.split(os.sep)))
    else:
        del(dirnames[:])


kwds = {'scripts': []}
if sys.platform.startswith('win') and using_setuptools:
    kwds['entry_points'] = {'console_scripts':
                                ["hexrd = hexrd.cli.main:main"]
                           }
else:
    kwds['scripts'].append('scripts/hexrd')
if ('bdist_wininst' in sys.argv) or ('bdist_msi' in sys.argv):
    kwds['scripts'].append('scripts/hexrd_win_post_install.py')


package_data = [
    'COPYING',
    'LICENSE',
    'wx/hexrd.png',
    'qt/*.ui',
    'data/materials.cfg',
    'data/all_materials.cfg',
    ]


data_files = [
    'share/example_config.yml',
    'share/calibrate_from_single_crystal.ipynb'
    ]


cmdclass = versioneer.get_cmdclass()
cmdclass['test'] = test

setup(
    name = 'hexrd',
    version = versioneer.get_version(),
    author = 'The HEXRD Development Team',
    author_email = 'praxes@googlegroups.com',
    description = 'hexrd diffraction data analysis',
    long_description = open('README.md').read(),
    license = 'LGPLv2',
    url = 'http://hexrd.readthedocs.org',
    classifiers = [
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Lesser General Public License v2 (LGPLv2)',
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering',
        ],
    ext_modules = ext_modules,
    packages = packages,
    install_requires = (
        'matplotlib',
        'numpy',
        'pyyaml',
        'scikit-learn',
        'scipy',
        'wxpython',
        ),
    package_data = {'hexrd': package_data},
    data_files = [('share/hexrd', data_files)],
    cmdclass = cmdclass,
    **kwds
    )

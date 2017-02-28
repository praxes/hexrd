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
from __future__ import print_function

import glob
import os
import sys

import numpy
np_include_dir = os.path.join(numpy.get_include(), 'numpy')

from setuptools import Command, Extension, find_packages, setup
from setuptools.command.build_ext import build_ext

import versioneer

cmdclass = versioneer.get_cmdclass()

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

cmdclass['test'] = test



_EXTENSION_FLAG_OVERRIDES = {
    'hexrd.xrd._transforms_CAPI': {
        'extra_compile_args': {
            'msvc': []
        }
    }
}
# handle overriding of options for specific compilers
class build_ext_override_keys( build_ext ):
    def build_extensions(self):
        c = self.compiler.compiler_type
        for e in self.extensions:
            overrides = _EXTENSION_FLAG_OVERRIDES.get(e._full_name)
            if overrides:
                for var, val_dict in overrides.iteritems():
                    val = val_dict.get(c)
                    if val is not None:
                        setattr(e, var, val)

        build_ext.build_extensions(self)

cmdclass['build_ext'] = build_ext_override_keys


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


srclist = ['transforms_CAPI.c', 'transforms_CFUNC.c']
srclist = [os.path.join('hexrd/transforms', f) for f in srclist]
transforms_mod = Extension(
    'hexrd.xrd._transforms_CAPI',
    sources=srclist,
    include_dirs=[np_include_dir],
    extra_compile_args=['-std=gnu99'],
    )

ext_modules = [sglite_mod, transforms_mod]

# use entry_points, not scripts:
entry_points = {
    'console_scripts': ["hexrd = hexrd.cli.main:main"]
    }

# only defining scripts so bdist_wininst can make an entry in the start menu
scripts = []
if ('bdist_wininst' in sys.argv) or ('bdist_msi' in sys.argv):
    scripts.append('scripts/hexrd_win_post_install.py')

data_files = [
    ('share/hexrd', glob.glob('share/*')),
    ]

package_data = {
    'hexrd': [
        'COPYING',
        'LICENSE',
        'data/*.cfg',
        'qt/resources/*.ui',
        'qt/resources/*.png',
        'wx/*.png',
        ]
    }

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
    packages = find_packages(),
    entry_points = entry_points,
    scripts = scripts,
    data_files = data_files,
    package_data = package_data,
    cmdclass = cmdclass,
    )

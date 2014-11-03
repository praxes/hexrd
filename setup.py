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

from distutils.core import setup, Extension
from distutils.cmd import Command
import os
import sys

import numpy
np_include_dir = os.path.join(numpy.get_include(), 'numpy')


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


scripts = []
if sys.platform.startswith('win'):
    # scripts calling multiprocessing must be importable
    import shutil
    shutil.copy('scripts/hexrd', 'scripts/hexrd_app.py')
    scripts.append('scripts/hexrd_app.py')
else:
    scripts.append('scripts/hexrd')
if ('bdist_wininst' in sys.argv) or ('bdist_msi' in sys.argv):
    scripts.append('scripts/hexrd_win_post_install.py')


# release.py contains version, authors, license, url, keywords, etc.
repo_root = os.path.dirname(os.path.abspath(__file__))
execfile(os.path.join(repo_root, 'hexrd','release.py'), globals())


package_data = [
    'COPYING',
    'LICENSE',
    'wx/hexrd.png',
    'data/materials.cfg',
    'data/all_materials.cfg',
    ]


data_files = [
    'share/example_config.yml',
    'share/calibrate_from_single_crystal.ipynb'
    ]


description = "hexrd diffraction data analysis"


setup(
    name = name,
    version = __version__,
    author = author,
    author_email = author_email,
    description = description,
    license = license,
    ext_modules = ext_modules,
    packages = packages,
    requires = (
        'python (>=2.7)',
        'numpy (>=1.4.0)',
        'scipy (>=0.7.0)',
        'wxpython (>= 2.8)',
        ),
    scripts = scripts,
    package_data = {'hexrd': package_data},
    data_files = [('share/hexrd', data_files)],
    cmdclass = {'test': test}
    )

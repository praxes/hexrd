# ============================================================
# Copyright (c) 2007-2012, Lawrence Livermore National Security, LLC. 
# Produced at the Lawrence Livermore National Laboratory. 
# Written by Joel Bernier <bernier2@llnl.gov> and others. 
# LLNL-CODE-529294. 
# All rights reserved.
# 
# This file is part of HEXRD. For details, see https://github.com/joelvbernier/hexrd.
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
import os

srclist = ['sgglobal.c','sgcb.c','sgcharmx.c','sgfile.c',
           'sggen.c','sghall.c','sghkl.c','sgltr.c','sgmath.c','sgmetric.c',
           'sgnorm.c','sgprop.c','sgss.c','sgstr.c','sgsymbols.c',
           'sgtidy.c','sgtype.c','sgutil.c','runtests.c','sglitemodule.c']
srclist = [os.path.join('hexrd/sglite', f) for f in srclist]

sglite_mod = Extension('hexrd.XRD.sglite', sources=srclist,
                   define_macros = [('PythonTypes', 1)])

ext_modules = [sglite_mod]

packages = []
for dirpath, dirnames, filenames in os.walk('.'):
    if '__init__.py' in filenames:
        packages.append('.'.join(dirpath.split(os.sep)))
    else:
        del(dirnames[:])

setup(
    author = 'Joel Bernier, et al.',
    author_email = 'bernier2@llnl.gov',
    description = 'High energy diffraction microscopy',
    ext_modules = ext_modules,
    name = 'hexrd',
    packages = packages,
    requires = (
        'python (>=2.6)',
        'numpy (>=1.4.0)',
        'scipy (>=0.7.0)',
        'wxpython (>= 2.8)',
        ),
)


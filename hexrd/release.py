# -*- coding: utf-8 -*-
"""Release data for the HEXRD project."""

#-----------------------------------------------------------------------------
#  Copyright (c) 2014, HEXRD Development Team
#  Copyright (c) 2012, Joel Bernier <>
#  Copyright (c) 2012, Lawrence Livermore National Security, LLC.
#
#  Distributed under the terms of the LGPL License, version 2.1.
#
#  The full license is in the file hexrd/LICENSE
#-----------------------------------------------------------------------------

# Name of the package for release purposes.  This is the name which labels
# the tarballs and RPMs made by distutils, so it's best to lowercase it.
name = 'hexrd'

# HEXRD version information.  An empty _version_extra corresponds to a full
# release.  'dev' as a _version_extra string means this is a development
# version
_version_major = 0
_version_minor = 2
_version_patch = 1
_version_extra = ''
#_version_extra = '.dev' # uncomment this during development
# _version_extra = 'rc1' # uncomment this for release candidates

# Construct full version string from these.
_ver = [_version_major, _version_minor, _version_patch]

__version__ = '.'.join(map(str, _ver))
if _version_extra:
    __version__ = __version__ + _version_extra

version_info = (_version_major, _version_minor, _version_patch, _version_extra)

description = "HEXRD: High energy x-ray diffraction"

long_description = \
"""
HEXRD provides a collection of resources for analysis of x-ray diffraction
data, especially high-energy x-ray diffraction. HEXRD is comprised of a
library and API for writing scripts, a command line interface, and an
interactive graphical user interface.
"""

license = 'LGPL'

authors = {'Joel' : ('Joel Bernier',''),
           'Don'    : ('Don Boyce',''),
           'Darren'   : ('Darren Dale','dsdale24@gmail.com'),
           }

author = 'The HEXRD Development Team'

author_email = 'praxes@googlegroups.com'

url = 'http://hexrd.readthedocs.org'

download_url = 'https://pypi.python.org/pypi/hexrd'

platforms = ['Linux','Mac OSX','Windows XP/Vista/7/8']

keywords = ['Synchrotron','X-ray diffraction']

classifiers = [
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: GNU Lesser General Public License v2 (LGPLv2)',
    'Programming Language :: Python',
    'Programming Language :: Python :: 2',
    'Programming Language :: Python :: 2.7',
    'Topic :: Scientific/Engineering',
    ]

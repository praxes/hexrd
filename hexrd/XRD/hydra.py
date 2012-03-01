#! /usr/bin/env python
# ============================================================
# Copyright (c) 2012, Lawrence Livermore National Security, LLC. 
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
#
"""Hydra detector tools

This is just a first pass at laying out a class for the hydra reader. 
Needs much more development.
"""
#
# ---------------------------------------------------CLASS:  Hydra
#
class Hydra(object):
    """Hydra image processing"""
    def __init__(self):
        """Constructor for Hydra."""
        #
        #  These arrays need four entries
        #
        self.readers = 4*[None]
	self.images  = 4*[None]
        #
        return
    #
    # ============================== API
    #
    #                     ========== Properties
    #
    #                     ========== Public Methods
    #
    def loadImages(self):
        """Load the four hydra images"""
        print 'loading images'
        reader1 = self.readers[0]
        aggMode = reader1.aggModeOp
        nrFrame = reader1.getNumberOfFrames() # number of reader frames
        if aggMode:
            rdFrames = nrFrame
        else:
            rdFrames = 1
            pass
        
        for i in range(4):
            ri = self.readers[i].makeReader()
            self.images[i] = ri.read(nframes= rdFrames, 
                                     sumImg = aggMode)
            print 'done %d' % i
            pass

        return
    #
    pass  # end class
#
# -----------------------------------------------END CLASS:  Hydra


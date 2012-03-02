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
"""Configuration for XRD application
"""
import platform

import wx

import matplotlib

matplotlib.use('WXAgg')
#  
#  * PLATFORM 
#
sysName   = platform.system()
onWindows = (sysName.lower() == 'windows')
onLinux   = (sysName.lower() == 'linux')
#
print 'system:  ', sysName
#
# ---------------------------------------------------CLASS:  WindowParameters
#
class FilesAndPaths:
    """File and path names"""
    DATA_ROOT = ''
    pass

class WindowParameters:
    """WindowParameters"""
    #
    #  COLOR SPECS
    #
    #  Ways to Set Color:
    #
    #  win.SetBackgroundColour(wxColour(0,0,255))
    #  win.SetBackgroundColour('BLUE')
    #  win.SetBackgroundColour('#0000FF')
    #  win.SetBackgroundColour((0,0,255))
    #
    #  * Containers
    #
    BG_COLOR_FRAME  = (224, 255, 224)
    BG_COLOR_PANEL  = (240, 240, 255) # light purple
    BG_COLOR_PANEL1 = 'ORANGE'
    #
    #  * Titlebars
    #
    #  - New naming convention on color variables; others
    #    are set below for backward compatibility
    #
    BG_COLOR_TITLEBAR_FRAME  = (64,  255,  64)
    BG_COLOR_TITLEBAR_PAGE   = (128, 255, 128)
    BG_COLOR_TITLEBAR_PANEL  = (192, 192,  64)
    BG_COLOR_TITLEBAR_PANEL1 = (192, 192, 128)
    #
    #
    #  * Miscellaneous
    #
    CANVAS_BG_COLOR = 'white'
    CANVAS_BG_COLOR = BG_COLOR_PANEL
    #
    #  MARGINS & PADDING
    #
    PAGE_SUBPANEL_PADDING = 10
    #
    #  ======================================== Old Variables
    #
    #  These are here until replaced in code.
    #
    TITLEBAR_BG_COLOR_FRAME  = BG_COLOR_TITLEBAR_FRAME
    TITLEBAR_BG_COLOR_PANEL  = BG_COLOR_TITLEBAR_PANEL
    TITLEBAR_BG_COLOR_PANEL1 = BG_COLOR_TITLEBAR_PANEL1
    #
    #  *-* new naming convention
    #
    BG_COLOR_FRAME_TITLEBAR  = BG_COLOR_TITLEBAR_FRAME 
    BG_COLOR_PAGE_TITLEBAR   = BG_COLOR_TITLEBAR_PAGE  
    BG_COLOR_PANEL_TITLEBAR  = BG_COLOR_TITLEBAR_PANEL 
    BG_COLOR_PANEL1_TITLEBAR = BG_COLOR_TITLEBAR_PANEL1
    #
    pass  # end class
#
# -----------------------------------------------END CLASS:  WindowParameters

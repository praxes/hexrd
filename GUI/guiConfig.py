#! /usr/bin/env python
#
"""Configuration for XRD application
"""
import platform

import wx

import matplotlib;  matplotlib.use('WXAgg')
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

#! /usr/bin/env python
#
#  $Id: mainApp.py 907 2011-07-22 22:27:02Z boyce6 $
#
"""Main application file
"""
import os, sys

import wx

import guiConfig
from xrdMainFrame import xrdMainFrame
#
#  mdef Modules
#
from XRD import detector as detectorModule

from XRD.Experiment import loadExp, ImageModes
#
# ---------------------------------------------------CLASS:  xrdApp
#
class xrdApp(wx.PySimpleApp):
    """xrdApp"""
    def __init__(self, *args):
	"""Constructor for xrdApp"""
	#
	wx.PySimpleApp.__init__(self)
        #
        # No command args for now, due to mac build issue (64bit, argv emulation)
        #
        f = ''
        #if len(args) == 0:
        #    f = ''
        #else:
        #    f = args[0]
        #    pass

        self.__makeData(f)
	#
	return

    def __makeData(self, inpFile):
        """Keep globals for easy access from subwindows"""
	#
        #  Data
        #
        #  * XRD Workspace
        #
        self.ws = loadExp(inpFile)
        #
        #  * Image Information
        #
        self.imgMode  = ImageModes.SINGLE_FRAME
        self.imgCal   = None
        self.imgSweep = None
         
        return

    def __getNotebook(self):
        """Return the notebook"""
        topWin = self.GetTopWindow()

        return topWin.nBook
    #
    # ============================== API
    #
    @property
    def imgFrame(self):
        """Image frame according to image mode"""
        if self.imgMode == ImageModes.SINGLE_FRAME:
            return self.imgCal
        else:
            return self.imgSweep
    
    def getCanvas(self):
        """Return the canvas panel"""
        
        return self.GetTopWindow().canvasPanel

    def updateFromExp(self):
        """Update GUI display"""
        self.GetTopWindow().updateFromExp()

        return
    
    pass # end class
#
# -----------------------------------------------END CLASS:  xrdApp

if __name__ == '__main__':
    #
    #  Run program stand-alone.
    #
    app = xrdApp(*sys.argv[1:])

    #if len(sys.argv) == 1:
    #    app = xrdApp()
    #else:
    #    app = xrdApp(*sys.argv[1:])
    #    pass
    #
    #  The main window cannot be imported until after the app 
    #  is instantiated due to the wx.ColourDatabase() call.
    #
    mFrame = xrdMainFrame(None, wx.NewId())
    app.SetTopWindow(mFrame)

    app.MainLoop()

    pass

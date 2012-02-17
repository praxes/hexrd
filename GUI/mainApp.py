#! /usr/bin/env python
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

	self.mFrame = None
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
    app.mFrame = xrdMainFrame(None, wx.NewId())
    app.SetTopWindow(app.mFrame)

    #if len(sys.argv) == 1:
    #    app = xrdApp()
    #else:
    #    app = xrdApp(*sys.argv[1:])
    #    pass
    #
    #  The main window cannot be imported until after the app 
    #  is instantiated due to the wx.ColourDatabase() call.
    #
    #
    # Splash screen.
    #
    splashFile = 'hexrd.png'
    splashDir = os.path.dirname(__file__)
    print 'splash:  ', os.path.join(splashDir, splashFile)
    splashImage = wx.Bitmap(os.path.join(splashDir, splashFile))
    #
    wx.SplashScreen(splashImage, wx.SPLASH_CENTRE_ON_PARENT|wx.SPLASH_TIMEOUT,
                    1000, app.mFrame)
    #
    # Main frame
    #
    app.mFrame.Show(True)
    #
    # GUI main loop
    #
    app.MainLoop()

    pass

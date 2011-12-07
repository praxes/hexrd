#! /usr/bin/env python
#
#  $Id$
#
"""Panel for selecting image reader
"""
import wx

from guiConfig    import WindowParameters as WP
from guiUtilities import makeTitleBar

from geReader import geReaderPanel

DET_GE      = 'ge'
DET_CHOICES = [DET_GE]
#
#  * mode choices for image reader
#
[MODE_CAL, MODE_SPOTS] = range(2)
#
# ---------------------------------------------------CLASS:  readerPanel
#
class readerPanel(wx.Panel):
    """readerPanel """
    def __init__(self, parent, id, **kwargs):
	"""Constructor for readerPanel."""
	#
	wx.Panel.__init__(self, parent, id, **kwargs)
        self.SetBackgroundColour(WP.BG_COLOR_PANEL)
	#
        #  Data
        #

        #
	#  Window Objects.
	#
        self.__makeObjects()
	#
	#  Bindings.
	#
	self.__makeBindings()
	#
	#  Sizing.
	#
	self.__makeSizers()
	#
	self.SetAutoLayout(True)
        self.SetSizerAndFit(self.sizer)
	#
	return
    #
    # ============================== Internal Methods
    #
    def __makeObjects(self):
        """Add interactors"""

        self.tbarSizer = makeTitleBar(self, ' Reader Panel ')
        #
        #  Add detector choice
        #
        self.det_cho = wx.Choice(self, wx.NewId(), 
                                 choices=DET_CHOICES)
        self.rdr_pan = geReaderPanel(self, wx.NewId())

        return

    def __makeBindings(self):
        """Bind interactors"""
	self.Bind(wx.EVT_CHOICE, self.OnDetChoice, self.det_cho)

        return

    def __makeSizers(self):
	"""Lay out the interactors"""
	
	self.sizer = wx.BoxSizer(wx.VERTICAL)
	#self.sizer.Add(self.titlebar, 0, wx.EXPAND|wx.ALIGN_CENTER)
	self.sizer.Add(self.tbarSizer, 0, wx.EXPAND|wx.ALIGN_CENTER)
	self.sizer.Add(self.det_cho,   0, wx.ALIGN_RIGHT)
	self.sizer.Add(self.rdr_pan,   1, wx.EXPAND|wx.ALIGN_RIGHT)

	return
    #
    # ============================== API
    #
    #                     ========== *** Access Methods
    #
    def updateFromExp(self):
        """Reset interactors for new exp"""
        self.rdr_pan.updateFromExp()

        return
    #
    #                     ========== *** Event Callbacks
    #
    def OnDetChoice(self, e):
        """On choice of detector"""
        val = self.det_cho.GetStringSelection()
        if not val == 'ge':
            msg = 'detector type "%s" not implemented: resetting' % val
            wx.MessageBox(msg)
            self.det_cho.SetStringSelection('ge')
            pass

        return
    
    pass # end class
#
# -----------------------------------------------END CLASS:  readerPanel

#! /usr/bin/env python
#
"""Panel for index
"""
import wx

from guiConfig    import WindowParameters as WP
from guiUtilities import makeTitleBar, callJoel
#
# ---------------------------------------------------CLASS:  indexPanel
#
class indexPanel(wx.Panel):
    """indexPanel """
    def __init__(self, parent, id, **kwargs):
	"""Constructor for indexPanel."""
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

        self.sz_titlebar = makeTitleBar(self, 'Indexing')
        self.hpage = callJoel(self)

        return

    def __makeBindings(self):
        """Bind interactors"""
        return

    def __makeSizers(self):
	"""Lay out the interactors"""
	
	self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.sizer.Add(self.sz_titlebar, 0, wx.EXPAND|wx.ALIGN_CENTER)
	self.sizer.Add(self.hpage,    1, wx.EXPAND|wx.ALIGN_CENTER)

	return
    #
    # ============================== API
    #
    def updateFromExp(self):
        """Update page"""
        return

    pass # end class
#
# -----------------------------------------------END CLASS:  indexPanel

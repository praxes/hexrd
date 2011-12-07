#! /usr/bin/env python
#
"""Notebook for XRD main window
"""
#
#  Notebook Styles and Events
#
"""
wxNB_TOP 	Place tabs on the top side.
wxNB_LEFT 	Place tabs on the left side.
wxNB_RIGHT 	Place tabs on the right side.
wxNB_BOTTOM 	Place tabs under instead of above the notebook pages.
wxNB_FIXEDWIDTH 	(Windows only) All tabs will have same width.
wxNB_MULTILINE 	(Windows only) There can be several rows of tabs.
wxNB_NOPAGETHEME 	(Windows only) Display a solid colour on notebook pages, and not a gradient, which can reduce performance.
wxNB_FLAT 	(Windows CE only) Show tabs in a flat style. 


EVT_NOTEBOOK_PAGE_CHANGED
EVT_NOTEBOOK_PAGE_CHANGING
"""
import wx

from guiConfig import WindowParameters as WP

from materialsPanel import matPanel
from readerPanel    import readerPanel

from detectorPanel import detectorPanel
from spotsPanel    import spotsPanel
from indexPanel    import indexPanel
from grainPanel    import grainPanel
#
# ---------------------------------------------------CLASS:  xrdNoteBook
#
class xrdNoteBook(wx.Notebook):
    """xrdNoteBook """
    def __init__(self, parent, id, **kwargs):
	"""Constructor for xrdNoteBook."""
	#
	wx.Notebook.__init__(self, parent, id, **kwargs)
        #
        #  Add pages to notebook, maintaining a dictionary
        #  of pages by page title.
        #
        self.pageDict = dict()
        #
        title = 'Materials'
        panel = matPanel(self, wx.NewId())
        self.materialsPanel = panel
        self.AddPage(panel, title)
        self.pageDict[title] = panel
        #
        title = 'Reader'
        self.readerPanel = readerPanel(self, wx.NewId())
        self.AddPage(self.readerPanel, title)
        self.pageDict[title] = self.readerPanel
        #
        title = 'Detector'
        self.detectorPanel = detectorPanel(self, wx.NewId())
        self.AddPage(self.detectorPanel, title)
        self.pageDict[title] = self.detectorPanel
        #
        title = 'Spots'
        self.spotsPanel = spotsPanel(self, wx.NewId())
        self.AddPage(self.spotsPanel, title)
        self.pageDict[title] = self.spotsPanel
        #
        #
        self.AddPage(indexPanel(self, wx.NewId()), 
                     'Indexing')
        self.AddPage(grainPanel(self, wx.NewId()), 
                     'Grains')
        #
        #  Make sure page is updated on page change.
        #
        self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, self.OnPageChange, self)
        
	#
	return
    #
    # ============================== Internal Methods
    #
    #
    # ============================== API
    #
    #                     ========== *** Access Methods
    #
    def updateFromExp(self):
        """Update all subwindows"""
        #for i in range(self.GetPageCount()):
        for i in range(3):
            self.GetPage(i).updateFromExp()
            pass

        return

    def OnPageChange(self, e):
        """Notebook page changed:  run pages update"""
        #print self.GetSelection()
        #print e.Selection, e.OldSelection
        #print dir(e)

        sel = e.Selection
        cp = self.GetPage(sel)
        cp.updateFromExp()
        e.Skip()

        return

    # These are out of date
    def getPage_geReader(self):  return self.GetPage(0)
    def getPage_Detector(self):  return self.GetPage(1)
    def getPage_Spot(self):      return self.GetPage(2)
    def getPage_Indexing(self):  return self.GetPage(3)
    def getPage_Grains(self):    return self.GetPage(4)

    def getPageByName(self, n):  return self.pageDict[n]
    #
    #                     ========== *** Event Callbacks
    #

    pass # end class
#
# -----------------------------------------------END CLASS:  xrdNoteBook

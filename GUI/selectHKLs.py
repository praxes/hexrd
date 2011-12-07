#! /usr/bin/env python
#
#  Copyright-Info-Goes-Here
#
"""Window used for selecting HKLs.
"""
import sys
import math

import wx
#import wx.lib.mixins.listctrl  as  listMixins

from guiConfig import WindowParameters as WP
#
# ---------------------------------------------------CLASS:  selectHKLsPanel.py
#
class selectHKLsPanel(wx.Panel):
    """selectHKLsPanel.py """
    def __init__(self, parent, id, mat, **kwargs):
	"""Constructor for selectHKLsPanel.py

        mat - crystallographic material
"""
	#
	wx.Panel.__init__(self, parent, id, **kwargs)
	#
        #  Data
        #
        self.mat = mat
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

        self.__makeTitleBar('Select HKLs')
        self.listctrl =  self.__makeListCtrl()

        return

    def __makeTitleBar(self, t):
        """Add titlebar"""
	self.titlebar = wx.StaticText(self, -1, t, 
					 style=wx.ALIGN_CENTER|wx.SIMPLE_BORDER) 
	self.titlebar.SetBackgroundColour(WP.TITLEBAR_BG_COLOR_PANEL)
        myToolTip = r"""
PANEL FOR ...
"""
        self.titlebar.SetToolTipString(myToolTip)

        return

    def __makeListCtrl(self):
	"""Make the list control"""
        #
	LStyle = wx.LC_REPORT
	#
	listctrl = wx.ListView(self, wx.NewId(), style=LStyle)
	listctrl.InsertColumn(0, 'HKL')
	listctrl.InsertColumn(1, 'd-spacing')
	listctrl.InsertColumn(2, '2-theta (deg)')
        listctrl.SetColumnWidth(0, 200)
        listctrl.SetColumnWidth(1, wx.LIST_AUTOSIZE_USEHEADER)
        listctrl.SetColumnWidth(2, wx.LIST_AUTOSIZE_USEHEADER)
        #
        #  Now add HKLs
        #
        pdat = self.mat.planeData
        hkls = pdat.hklDataList # property
        excl = pdat.exclusions

        print 'exclusions:\n', excl
        print len(hkls)
        print hkls[0]
        
        n = pdat.getNhklRef()
        for i in range(n):
            hklData = hkls[i]
            hkl = hklData['hkl']
            hklStr = '(%d, %d, %d)' % (hkl[0], hkl[1], hkl[2])
            index = listctrl.InsertStringItem(sys.maxint, hklStr)
            dspace = '%.6g' % hklData['dSpacings']
            tth    = hklData['tTheta'] * (180/math.pi)
            tTheta = '%.6g' % tth
            listctrl.SetStringItem(index, 1, dspace)
            listctrl.SetStringItem(index, 2, tTheta)
            #
            #  Show exclusions by background color
            #
            if excl[i]:
                listctrl.SetItemBackgroundColour(index, 'grey')
                pass
            pass
        #
        #  Save certain data
        #
        self.nhkls   = n
        self.hkls    = hkls
        self.exclude = excl

	return listctrl

    def __makeBindings(self):
        """Bind interactors"""
        #self.Bind(wx.EVT_LISTBOX, self.OnSelection, wx.EVT_LIST_ITEM_SELECTED)

        return

    def __makeSizers(self):
	"""Lay out the interactors"""
	
	self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.sizer.Add(self.titlebar, 0, wx.EXPAND|wx.ALIGN_CENTER)
	self.sizer.Add(self.listctrl, 1, wx.EXPAND|wx.ALIGN_CENTER)

	return
    #
    # ============================== API
    #
    #                     ========== *** Access Methods
    #
    def getExclusions(self):
        e = self.exclude
        for i in range(self.nhkls):
            e[i] = not self.listctrl.IsSelected(i)
            pass

        return e
    #
    #                     ========== *** Event Callbacks
    #
    def OnSelection(self, e):
        """List control selection has changed"""
        
        return

    
    pass # end class
#
# -----------------------------------------------END CLASS:  selectHKLsPanel.py
#
# ---------------------------------------------------CLASS:  selectHKLsDialog
#
class selectHKLsDialog(wx.Dialog):
    """selectHKLsDialog """
    def __init__(self, parent, id, mat, **kwargs):
	"""Constructor for selectHKLsDialog"""
	#
        myCaption = 'HKL Selector'
        myStyle   = wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER
	wx.Dialog.__init__(self, parent, id, 
                           style=myStyle, title=myCaption)
	#
	#  Data Objects.
	#
	
	#
	#  Windows.
	#
	self.titlebar = wx.StaticText(self, -1, 'selectHKLsDialog', 
				      style=wx.ALIGN_CENTER|wx.SIMPLE_BORDER) 
        self.dataPanel = selectHKLsPanel(self, wx.NewId(), mat)
        self.dataPanel.SetMinSize((400,400))
	#
	#  Bindings.
	#
	self._makeBindings()
	#
	#  Sizing.
	#
	self._makeSizers()
	#
	self.SetAutoLayout(True)
        self.SetSizerAndFit(self.sizer)
	#
	#  Events.
	#
	
	#
	return
    #
    # ============================== Internal Methods
    #
    def _makeBindings(self):
	"""Bind interactors to functions"""

	return

    def _makeSizers(self):
	"""Lay out windows"""
	self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.sizer.Add(self.titlebar,  0, wx.EXPAND|wx.ALIGN_CENTER)
	self.sizer.Add(self.dataPanel, 1, wx.EXPAND|wx.ALIGN_CENTER)
        #
        #  Add buttons
        #
        bs = self.CreateButtonSizer(wx.OK|wx.CANCEL)
        self.sizer.Add(bs, 0, wx.ALIGN_CENTER)
        #
        self.sizer.Show(self.titlebar, False)

	return
    #
    # ============================== API
    #
    #                     ========== *** Access Methods
    #
    def getExclusions(self):
        return self.dataPanel.getExclusions()
    #
    #                     ========== *** Event Callbacks
    #
    pass # end class
#
# -----------------------------------------------END CLASS:  selectHKLsDialog
#

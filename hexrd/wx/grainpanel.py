#! /usr/bin/env python
# ============================================================
# Copyright (c) 2012, Lawrence Livermore National Security, LLC. 
# Produced at the Lawrence Livermore National Laboratory. 
# Written by Joel Bernier <bernier2@llnl.gov> and others. 
# LLNL-CODE-529294. 
# All rights reserved.
# 
# This file is part of HEXRD. For details on dowloading the source,
# see the file COPYING.
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
"""Panel for grain
"""
import wx

from hexrd.GUI.guiConfig import WindowParameters as WP
from hexrd.GUI.guiUtilities import makeTitleBar
#
# ---------------------------------------------------CLASS:  GrainPanel
#
class GrainPanel(wx.Panel):
    """GrainPanel """
    def __init__(self, parent, id, **kwargs):
	"""Constructor for GrainPanel."""
	#
	wx.Panel.__init__(self, parent, id, **kwargs)
        self.SetBackgroundColour(WP.BG_COLOR_PANEL)
	#
        #  Data
        #

        #
	#  Window Objects.
	#
        self.__make_objects()
	#
	#  Bindings.
	#
	self.__make_bindings()
	#
	#  Sizing.
	#
	self.__make_sizers()
	#
	self.SetAutoLayout(True)
        self.SetSizerAndFit(self.sizer)
	#
	return
    #
    # ============================== Internal Methods
    #
    def __make_objects(self):
        """Add interactors"""

        self.sz_titlebar = makeTitleBar(self, 'Grains')
        self.glist_pan = GrainListSubPanel(self, wx.NewId())
        self.refine_pan = RefinementSubPanel(self, wx.NewId())

        return

    def __make_bindings(self):
        """Bind interactors"""
        return

    def __make_sizers(self):
	"""Lay out the interactors"""
	
	self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.sizer.Add(self.sz_titlebar, 0, wx.EXPAND|wx.ALIGN_CENTER)
	self.sizer.Add(self.glist_pan,   1, wx.EXPAND|wx.ALIGN_CENTER)
	self.sizer.Add(self.refine_pan,  1, wx.EXPAND|wx.ALIGN_CENTER)

	return
    #
    # ============================== API
    #
    def updateFromExp(self):
        """Update page"""
        return

    
    pass # end class
#
# -----------------------------------------------END CLASS:  GrainPanel
# ---------------------------------------------------CLASS:  GrainListSubPanel
#
class GrainListSubPanel(wx.Panel):
    """GrainListSubPanel """
    def __init__(self, parent, id, **kwargs):
	"""Constructor for GrainListSubPanel."""
	#
	wx.Panel.__init__(self, parent, id, **kwargs)
	#
        #  Data
        #

        #
	#  Window Objects.
	#
        self.__make_objects()
	#
	#  Bindings.
	#
	self.__make_bindings()
	#
	#  Sizing.
	#
	self.__make_sizers()
	#
	self.SetAutoLayout(True)
        self.SetSizerAndFit(self.sizer)
	#
	return
    #
    # ============================== Internal Methods
    #
    def __make_objects(self):
        """Add interactors"""

        self.tbarSizer = makeTitleBar(self, 'Grain List',
                                      color=WP.BG_COLOR_PANEL1_TITLEBAR)
        self.listctrl = self.__make_listctrl()
        self.ginfo_but = wx.Button(self, wx.NewId(), 'show info')

        return
    
    def __make_listctrl(self):
	"""Make the list control"""
        #
	LStyle = wx.LC_REPORT
	#
	listctrl = wx.ListView(self, wx.NewId(), style=LStyle)
	listctrl.InsertColumn(0, ' ID ')
	listctrl.InsertColumn(1, 'completeness')
        listctrl.SetColumnWidth(0, wx.LIST_AUTOSIZE_USEHEADER)
        listctrl.SetColumnWidth(1, wx.LIST_AUTOSIZE_USEHEADER)

        return listctrl

    def __make_bindings(self):
        """Bind interactors"""
        self.Bind(wx.EVT_BUTTON, self.OnTable, self.ginfo_but)
        
        return

    def __make_sizers(self):
	"""Lay out the interactors"""
	
	self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.sizer.Add(self.tbarSizer, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.sizer.Add(self.listctrl)
        self.sizer.Add(self.ginfo_but, 0, wx.ALIGN_CENTER)

	return
    #
    # ============================== API
    #
    #                     ========== *** Access Methods
    #

    #
    #                     ========== *** Event Callbacks
    #
    def OnTable(self, e):
        """Raise table with grain info"""
        dlg = GrainTableDlg(self, wx.NewId())
        dlg.ShowModal()
        
        return
    
    pass # end class
#
# -----------------------------------------------END CLASS:  GrainListSubPanel
# ---------------------------------------------------CLASS:  RefinementSubPanel
#
class RefinementSubPanel(wx.Panel):
    """RefinementSubPanel """
    def __init__(self, parent, id, **kwargs):
	"""Constructor for RefinementSubPanel."""
	#
	wx.Panel.__init__(self, parent, id, **kwargs)
	#
        #  Data
        #

        #
	#  Window Objects.
	#
        self.__make_objects()
	#
	#  Bindings.
	#
	self.__make_bindings()
	#
	#  Sizing.
	#
	self.__make_sizers()
	#
	self.SetAutoLayout(True)
        self.SetSizerAndFit(self.sizer)
	#
	return
    #
    # ============================== Internal Methods
    #
    def __make_objects(self):
        """Add interactors"""

        self.tbarSizer = makeTitleBar(self, 'Grain Refinement',
                                        color=WP.BG_COLOR_PANEL1_TITLEBAR)
        self.refine_but  = wx.Button(self, wx.NewId(), 'refine')

        return

    def __make_bindings(self):
        """Bind interactors"""
        return

    def __make_sizers(self):
	"""Lay out the interactors"""
	
	self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.sizer.Add(self.tbarSizer, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.sizer.Add(self.refine_but, 0, wx.ALIGN_CENTER|wx.TOP, 5)

	return
    #
    # ============================== API
    #
    #                     ========== *** Access Methods
    #

    #
    #                     ========== *** Event Callbacks
    #
    
    pass # end class
#
# -----------------------------------------------END CLASS:  RefinementSubPanel
# ---------------------------------------------------CLASS:  GrainTablePanel
#
class GrainTablePanel(wx.Panel):
    """GrainTablePanel """
    def __init__(self, parent, id, **kwargs):
	"""Constructor for GrainTablePanel."""
	#
	wx.Panel.__init__(self, parent, id, **kwargs)
	#
        #  Data
        #

        #
	#  Window Objects.
	#
        self.__make_objects()
	#
	#  Bindings.
	#
	self.__make_bindings()
	#
	#  Sizing.
	#
	self.__make_sizers()
	#
	self.SetAutoLayout(True)
        self.SetSizerAndFit(self.sizer)
	#
	return
    #
    # ============================== Internal Methods
    #
    def __make_objects(self):
        """Add interactors"""

        self.tbarSizer = makeTitleBar(self, 'Grain Table')
        self.listctrl = self.__make_listctrl()
        return

    def __make_listctrl(self):
	"""Make the list control"""
        #
	LStyle = wx.LC_REPORT
	#
	listctrl = wx.ListView(self, wx.NewId(), style=LStyle)
	listctrl.InsertColumn(0, 'spot ID')
	listctrl.InsertColumn(1, '  hkl  ')
        listctrl.SetColumnWidth(0, wx.LIST_AUTOSIZE_USEHEADER)
        listctrl.SetColumnWidth(1, wx.LIST_AUTOSIZE_USEHEADER)

        return listctrl

    def __make_bindings(self):
        """Bind interactors"""
        return

    def __make_sizers(self):
	"""Lay out the interactors"""
	
	self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.sizer.Add(self.tbarSizer, 0, wx.EXPAND|wx.ALIGN_CENTER)
	self.sizer.Add(self.listctrl, 1, wx.EXPAND|wx.ALIGN_CENTER)

	return
    #
    # ============================== API
    #
    #                     ========== *** Access Methods
    #

    #
    #                     ========== *** Event Callbacks
    #
    
    pass # end class
#
# -----------------------------------------------END CLASS:  GrainTablePanel
# ---------------------------------------------------CLASS:  GrainTableDlg
#
class GrainTableDlg(wx.Dialog):
    """GrainTableDlg """
    def __init__(self, parent, id, **kwargs):
	"""Constructor for GrainTableDlg"""
	#
        myStyle = wx.RESIZE_BORDER|wx.DEFAULT_DIALOG_STYLE
	wx.Dialog.__init__(self, parent, id, style=myStyle)
	#
	#  Data Objects.
	#
	
	#
	#  Windows.
	#
	self.titlebar = wx.StaticText(self, -1, 'Grain Table Dialog', 
				      style=wx.ALIGN_CENTER|wx.SIMPLE_BORDER)
        self.panel = GrainTablePanel(self, wx.NewId())
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
	self.sizer.Add(self.titlebar, 0, wx.EXPAND|wx.ALIGN_CENTER)
	self.sizer.Add(self.panel, 1, wx.EXPAND|wx.ALIGN_CENTER)

	return
    #
    # ============================== API
    #
    #                     ========== *** Access Methods
    #

    #
    #                     ========== *** Event Callbacks
    #
    

    pass # end class
#
# -----------------------------------------------END CLASS:  GrainTableDlg
#

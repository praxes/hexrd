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
"""Panel for index
"""
import wx

from hexrd.GUI.guiConfig    import WindowParameters as WP
from hexrd.GUI.guiUtilities import makeTitleBar, callJoel
#
# ---------------------------------------------------CLASS:  indexPanel
#
class indexPanel(wx.Panel):
    """indexPanel """
    #
    # Class data
    #
    INDEX_CHOICES = ['Fiber Search', 'GrainSpotter']
    INDEX_CHOICE_IDS = [IND_FIBER, IND_GSPOT] = range(2)
    
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
        self.ChooseMethod(None)
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
        self.method_cho = wx.Choice(self, wx.NewId(), choices=self.INDEX_CHOICES)
        self.run_but  = wx.Button(self, wx.NewId(), 'Run Indexer')

        self.fiber_pan = FiberSearchPanel(self, wx.NewId())
        self.gspot_pan = GrainSpotterPanel(self, wx.NewId())
        
        return

    def __makeBindings(self):
        """Bind interactors"""
	self.Bind(wx.EVT_CHOICE, self.ChooseMethod, self.method_cho)
        return

    def __makeSizers(self):
	"""Lay out the interactors"""
        # Top sizer for method and run
	self.topsizer = wx.BoxSizer(wx.HORIZONTAL)
        self.topsizer.Add(self.method_cho, 0)
        self.topsizer.Add(self.run_but,  0, wx.LEFT, 5)
    
	self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.sizer.Add(self.sz_titlebar, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.sizer.Add(self.topsizer, 0)
        self.sizer.Add(self.fiber_pan, 1, wx.EXPAND|wx.ALIGN_CENTER)
        self.sizer.Add(self.gspot_pan, 1, wx.EXPAND|wx.ALIGN_CENTER)

	return
    #
    # ============================== API
    #
    #                     ========== *** Access Methods
    #
    def updateFromExp(self):
        """Update page"""
        return
    #
    #                     ========== *** Event Callbacks
    #
    def ChooseMethod(self, e):
        """Choose method button has been pressed"""
        if self.method_cho.GetSelection() == self.IND_FIBER:
            self.fiber_pan.Enable()
            self.gspot_pan.Disable()
        else:
            self.fiber_pan.Disable()
            self.gspot_pan.Enable()
            
        return
    
    pass # end class
#
# -----------------------------------------------END CLASS:  indexPanel
# ---------------------------------------------------CLASS:  FiberSearchPanel
#
class FiberSearchPanel(wx.Panel):
    """Handles fiber search options """
    def __init__(self, parent, id, **kwargs):
	"""Constructor for FiberSearchPanel"""
	#
	wx.Panel.__init__(self, parent, id, **kwargs)
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

        self.tbarSizer = makeTitleBar(self, 'Fiber Search Options', 
                                      color=WP.BG_COLOR_PANEL1_TITLEBAR)

        return

    def __makeBindings(self):
        """Bind interactors"""
        return

    def __makeSizers(self):
	"""Lay out the interactors"""
	
	self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.sizer.Add(self.tbarSizer, 0, wx.EXPAND|wx.ALIGN_CENTER)

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
# -----------------------------------------------END CLASS:  FiberSearchPanel
# ---------------------------------------------------CLASS:  GrainSpotterPanel
#
class GrainSpotterPanel(wx.Panel):
    """Handles grain spotter options """
    def __init__(self, parent, id, **kwargs):
	"""Constructor for GrainSpotterPanel"""
	#
	wx.Panel.__init__(self, parent, id, **kwargs)
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

        self.tbarSizer = makeTitleBar(self, 'Grain Spotter Options', 
                                      color=WP.BG_COLOR_PANEL1_TITLEBAR)

        return

    def __makeBindings(self):
        """Bind interactors"""
        return

    def __makeSizers(self):
	"""Lay out the interactors"""
	
	self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.sizer.Add(self.tbarSizer, 0, wx.EXPAND|wx.ALIGN_CENTER)

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
# -----------------------------------------------END CLASS:  GrainSpotterPanel

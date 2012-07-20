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
        self.method_cho = wx.Choice(self, wx.NewId(), choices=['Fiber Search', 'GrainSpotter'])
        self.run_but  = wx.Button(self, wx.NewId(), 'Run Indexer')
        
        return

    def __makeBindings(self):
        """Bind interactors"""
	#self.Bind(wx.EVT_CHOICE, self.OnChoice, self.choice)
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

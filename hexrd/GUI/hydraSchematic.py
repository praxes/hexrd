#! /usr/bin/env python
# ============================================================
# Copyright (c) 2012, Lawrence Livermore National Security, LLC. 
# Produced at the Lawrence Livermore National Laboratory. 
# Written by Joel Bernier <bernier2@llnl.gov> and others. 
# LLNL-CODE-529294. 
# All rights reserved.
# 
# This file is part of HEXRD. For details, see https://github.com/joelvbernier/hexrd.
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
"""Window showing layout of detectors
"""
import wx

from hexrd.GUI.guiConfig import WindowParameters as WP
#
# ---------------------------------------------------CLASS:  HydraSchematicPanel
#
class HydraSchematicPanel(wx.Panel):
    """HydraSchematicPanel """
    def __init__(self, parent, id, **kwargs):
	"""Constructor for HydraSchematicPanel."""
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

        self.__makeTitleBar('Hydra Schematic')

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

    def __makeBindings(self):
        """Bind interactors"""
        return

    def __makeSizers(self):
	"""Lay out the interactors"""
	
	self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.sizer.Add(self.titlebar, 0, wx.EXPAND|wx.ALIGN_CENTER)

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
# -----------------------------------------------END CLASS:  HydraSchematicPanel

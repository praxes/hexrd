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
"""Panel for selecting image reader
"""
import wx

from hexrd.GUI.guiConfig    import WindowParameters as WP
from hexrd.GUI.guiUtilities import makeTitleBar

from hexrd.GUI.geReader import geReaderPanel

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

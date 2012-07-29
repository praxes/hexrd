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
"""Panel for spots
"""
import wx

from hexrd.GUI.guiConfig    import WindowParameters as WP
from hexrd.GUI.guiUtilities import makeTitleBar

from hexrd.XRD.crystallography    import processWavelength
#
# ---------------------------------------------------CLASS:  spotsPanel
#
class spotsPanel(wx.Panel):
    """spotsPanel """
    def __init__(self, parent, id, **kwargs):
	"""Constructor for spotsPanel."""
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

        self.tbarSizer = makeTitleBar(self, 'Spots')

        app = wx.GetApp(); exp = app.ws

        # Booleans

        self.disc_box  = wx.CheckBox(self, wx.NewId(), 'Discard at bounds')
        self.bbox_box  = wx.CheckBox(self, wx.NewId(), 'Keep in bounding box')
        self.pado_box  = wx.CheckBox(self, wx.NewId(), 'Pad Omega')
        self.pads_box  = wx.CheckBox(self, wx.NewId(), 'Pad Spots')


        # Threshold

        self.thresh_lab = wx.StaticText(
            self, wx.NewId(), 
            'Threshold', style=wx.ALIGN_CENTER)
        self.thresh_txt = wx.TextCtrl(
            self, wx.NewId(), value='', 
            style=wx.RAISED_BORDER)

        # Min PX

        self.minpx_lab = wx.StaticText(
            self, wx.NewId(), 
            'Min PX', style=wx.ALIGN_CENTER)
        self.minpx_txt = wx.TextCtrl(
            self, wx.NewId(), value='', 
            style=wx.RAISED_BORDER)

        self.run  = wx.Button(self, wx.NewId(), 'Find Spots')        

        return

    def __makeBindings(self):
        """Bind interactors"""
        self.Bind(wx.EVT_TEXT_ENTER, self.OnThreshold,  self.thresh_txt)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnMinPX,      self.minpx_txt)

        self.Bind(wx.EVT_BUTTON, self.OnRun, self.run)

        return

    def __makeSizers(self):
	"""Lay out the interactors"""
        padtop = 10
        
        # ========== Checkboxes
        
        self.cbsizer = wx.BoxSizer(wx.VERTICAL)  
        self.cbsizer.Add(self.disc_box, 0, wx.ALIGN_LEFT)
        self.cbsizer.Add(self.bbox_box, 0, wx.ALIGN_LEFT)
        self.cbsizer.Add(self.pado_box, 0, wx.ALIGN_LEFT)
        self.cbsizer.Add(self.pads_box, 0, wx.ALIGN_LEFT)
        
        # ========== Valued options

        nrow = 0; ncol = 2; padx = 5; pady = 5
        self.fgSizer = wx.FlexGridSizer(nrow, ncol, padx, pady) 
        self.fgSizer.AddGrowableCol(1, 1)
        #  threshold
        self.fgSizer.Add(self.thresh_lab,       0, wx.ALIGN_RIGHT)
        self.fgSizer.Add(self.thresh_txt,       0, wx.ALIGN_RIGHT)
        #  min PX
        self.fgSizer.Add(self.minpx_lab,       0, wx.ALIGN_RIGHT)
        self.fgSizer.Add(self.minpx_txt,       0, wx.ALIGN_RIGHT)

        # ========== Main Sizer

	self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.sizer.Add(self.tbarSizer, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.sizer.Add(self.cbsizer,   0, wx.ALIGN_RIGHT|wx.TOP, padtop)
        self.sizer.Add(self.fgSizer,   1, wx.ALIGN_RIGHT|wx.TOP, padtop)
        self.sizer.Add(self.run,       0, wx.ALIGN_RIGHT)

	return
    #
    # ============================== API
    #
    #                     ========== *** Access Methods
    #
    def updateFromExp(self):
        """Update all subwindows"""
        app  = wx.GetApp(); exp  = app.ws

        self.minpx_txt.SetValue(str(0))
        self.thresh_txt.SetValue(str(0))

        #AngstromTimesKev = processWavelength(1.0)
        #ang = exp.phase0.planeData.wavelength
        #kev = AngstromTimesKev/ang
        #self.beam_txt.SetValue(str(kev))

        #self.mat_cho.SetStringSelection(exp.phase0.name)

        return
    #
    #                     ========== *** Event Callbacks
    #
    def OnMinPX(self, evt):
        """Callback for minpx_txt control"""
        app  = wx.GetApp(); exp  = app.ws

        return

    def OnThreshold(self, evt):
        """Callback for beam_txt control"""

        return

    def OnBeamEnergy(self, evt):
        """Callback for beam_txt control"""

        return

    def OnChooseMat(self, evt):
        """Choose sweep material"""
        app = wx.GetApp(); exp = app.ws

        matName = self.mat_cho.GetStringSelection()
        #exp.phase0 = exp.matDict[matName]
        #exp.phase0.planeData.tThMax = exp.detector.getTThMax()

        return

    def OnRun(self, evt):
        """Callback for run"""
        exp = wx.GetApp().ws
        #
        # Fill in spot options from the form
        #
        opts = exp.spotOpts
        #
        opts.thresh = int(self.thresh_txt.GetValue())
        opts.minpx = int(self.minpx_txt.GetValue())
        #
        exp.findSpots()
        
        return
    
    pass # end class
#
# -----------------------------------------------END CLASS:  spotsPanel

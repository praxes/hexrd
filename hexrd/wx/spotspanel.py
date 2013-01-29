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

from hexrd.wx.guiconfig    import WindowParameters as WP
from hexrd.wx.guiutil import makeTitleBar

from hexrd.xrd.crystallography    import processWavelength
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

        # Material Choice

        self.mat_lab = wx.StaticText(self, -1, 'material',
                                     style=wx.ALIGN_CENTER|wx.SIMPLE_BORDER)
        self.mat_cho = wx.Choice(self, wx.NewId(), choices=exp.matNames)

        # Beam Energy

        self.beam_lab = wx.StaticText(
            self, wx.NewId(),
            'Beam Energy (keV)', style=wx.ALIGN_CENTER)
        self.beam_txt = wx.TextCtrl(
            self, wx.NewId(), value='',
            style=wx.RAISED_BORDER)

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
        self.Bind(wx.EVT_CHOICE,     self.OnChooseMat,  self.mat_cho)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnBeamEnergy, self.beam_txt)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnThreshold,  self.thresh_txt)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnMinPX,      self.minpx_txt)

        self.Bind(wx.EVT_BUTTON, self.OnRun, self.run)

        return

    def __makeSizers(self):
        """Lay out the interactors"""
        nrow = 3; ncol = 2; padx = 5; pady = 5
        self.fgSizer = wx.FlexGridSizer(nrow, ncol, padx, pady)
        self.fgSizer.AddGrowableCol(1, 1)
        #  1. material selector
        self.fgSizer.Add(self.mat_lab,       0, wx.ALIGN_RIGHT)
        self.fgSizer.Add(self.mat_cho,       0, wx.ALIGN_RIGHT)
        #  2. beam energy
        self.fgSizer.Add(self.beam_lab,       0, wx.ALIGN_RIGHT)
        self.fgSizer.Add(self.beam_txt,       0, wx.ALIGN_RIGHT)
        #  3. threshold
        self.fgSizer.Add(self.thresh_lab,       0, wx.ALIGN_RIGHT)
        self.fgSizer.Add(self.thresh_txt,       0, wx.ALIGN_RIGHT)
        #  4. min PX
        self.fgSizer.Add(self.minpx_lab,       0, wx.ALIGN_RIGHT)
        self.fgSizer.Add(self.minpx_txt,       0, wx.ALIGN_RIGHT)
        #
        #  ========== Main Sizer
        #
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.tbarSizer, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.sizer.Add(self.fgSizer,   1, wx.ALIGN_RIGHT)
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
        #wx.GetApp().ws.findSpots()
        return

    pass # end class
#
# -----------------------------------------------END CLASS:  spotsPanel

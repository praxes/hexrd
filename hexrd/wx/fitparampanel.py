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
"""Panel for fit parameters

Includes EVT_FIT_PARAM as new event.
"""
import wx

from hexrd.wx.guiconfig import WindowParameters as WP
#
#  Module Data
#
EVT_FIT_PARAM_T = wx.NewEventType()
EVT_FIT_PARAM   = wx.PyEventBinder(EVT_FIT_PARAM_T, 1)
#
# ---------------------------------------------------CLASS:  fitParamPanel
#
class fitParamPanel(wx.Panel):
    """fitParamPanel """
    def __init__(self, parent, id, fParams, **kwargs):
        """Constructor for fitParamPanel.

        NOTE that fitParams is argument.
"""
        #
        wx.Panel.__init__(self, parent, id, **kwargs)
        #
        #  Data
        #
        self.fParams = fParams
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

        #self.__makeTitleBar('fitParamPanel')

        # Now add several items for each parameter

        self.nColumns = 2 # 2 items per parameter
        self.rowDict  = dict()

        for p in self.fParams:
            name  = p.getProp('name')
            valu  = p.getProp('value')
            cbox  = wx.CheckBox(self, wx.NewId(), name)
            spin  = wx.SpinCtrl(self, wx.NewId(), str(valu), initial=50, name=name)
            self.rowDict[name] = [cbox, spin]
            pass

        return

    def __makeTitleBar(self, t):
        """Add titlebar"""
        self.titlebar = wx.StaticText(self, -1, t,
                                      style=wx.ALIGN_CENTER|wx.SIMPLE_BORDER)
        self.titlebar.SetBackgroundColour(WP.TITLEBAR_BG_COLOR)
        myToolTip = r"""
        PANEL FOR managing data for fit parameters
        """
        self.titlebar.SetToolTip(myToolTip)

        return

    def __makeBindings(self):
        """Bind interactors"""
        for p in self.fParams:
            name = p.getProp('name')
            row  = self.rowDict[name]
            spin = row[1]
            self.Bind(wx.EVT_SPINCTRL, self.UpdateParamValue, spin)

            pass

        return

    def __makeSizers(self):
        """Lay out the interactors"""
        nrow = self.fParams.getNumParam()
        ncol = self.nColumns
        padx = pady =5 # px

        self.fgsizer = wx.FlexGridSizer(nrow, ncol, padx, pady) # m x n, paddings
        self.fgsizer.AddGrowableCol(0, 1)
        self.fgsizer.AddGrowableCol(1, 1)
        for p in self.fParams:
            name = p.getProp('name')
            row  = self.rowDict[name]
            self.fgsizer.Add(row[0], 1, wx.EXPAND)
            self.fgsizer.Add(row[1], 1, wx.EXPAND)
            pass

        self.sizer = self.fgsizer

        return
    #
    # ============================== API
    #
    def UpdateParamValue(self, e):
        """Update parameter value from spin control"""
        s = e.GetEventObject()
        n = s.GetName()
        p = self.fParams.getParam(n)
        i = s.GetValue()

        vmin = p.getProp('min')
        vmax = p.getProp('max')
        v    =  vmin + (i/100.0)*(vmax - vmin)

        # note:  SetValueString is wxPython (SetValue is overloaded in wx widgets)
        s.SetValueString(str(v))

        fpEvent = fitParamEvent(EVT_FIT_PARAM_T, self.GetId())
        self.GetEventHandler().ProcessEvent(fpEvent)

        return

    pass # end class
#
# -----------------------------------------------END CLASS:  fitParamPanel
# ---------------------------------------------------CLASS:  fitParamEvent
#
class fitParamEvent(wx.PyCommandEvent):
    """fitParamEvent

    Indicates that the fit parameters have been updated.
"""
    def __init__(self, evtType, id):
        """Constructor for fitParamEvent"""
        wx.PyCommandEvent.__init__(self, evtType, id)
        self.myVal = None

        return
    #
    # ============================== API
    #

    #
    pass  # end class
#
# -----------------------------------------------END CLASS:  fitParamEvent

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
"""Panel for detector calibration

STATUS

* Needs to work with the MAR detector

TO DO

* Panel for modifying plane-data (in progress)

"""
import os

import wx
import numpy
#
#  XRD Modules
#
from hexrd.xrd            import detector
from hexrd.xrd.Material   import *
from hexrd.xrd.Experiment import FitModes

from hexrd.xrd.crystallography import dUnit as WAVELENGTH_UNIT
from hexrd.xrd.crystallography import processWavelength
#
#  wx Modules
#
from hexrd.wx.guiConfig    import WindowParameters as WP
from hexrd.wx.guiUtilities import ResetChoice, AddSpacer, EmptyWindow, makeTitleBar
from hexrd.wx.FloatControl import *
from hexrd.wx.LogWindows   import logWindow
from hexrd.wx.ringSubPanel import ringPanel
from hexrd.wx.selectHKLs   import selectHKLsDialog as hklsDlg #TBR

# AngstromTimesKev = 12.39854 # from APS site (lose digits accuracy this way)
AngstromTimesKev = processWavelength(1.0)
#
# ---------------------------------------------------CLASS:  plotPanel
#
class detectorPanel(wx.Panel):
    """plotPanel """
    def __init__(self, parent, id, **kwargs):
	"""Constructor for plotPanel."""
	#
	wx.Panel.__init__(self, parent, id, **kwargs)
        self.SetBackgroundColour(WP.BG_COLOR_PANEL)
	#
        #  Data
        #
        #  *** NONE YET ***
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
        #
        exp = wx.GetApp().ws
        #
        #  Titlebar
        #
        self.tbarSizer = makeTitleBar(self, 'Detector Calibration')
        #
        #  Material Selection
        #
        self.mats_lab = wx.StaticText(self, wx.NewId(),
                                        'Active Material', style=wx.ALIGN_CENTER)
        self.mats_cho = wx.Choice(self, wx.NewId(),
                                  choices=[m.name for m in exp.matList])
        #
        #  Rings panel
        #
        self.ring_pan = ringPanel(self, wx.NewId())
        #
        #  II.  Geometry
        #
        #  Will have a checkbox and a spin control for each parameter.
        #
        self.detLabSizer = makeTitleBar(self, 'Detector Parameters',
                                          color=WP.TITLEBAR_BG_COLOR_PANEL1)
        deltaBoxTip = "Use this box to set the increment for the spinner to the left"
        app = wx.GetApp()
        det = app.ws.detector

        name = 'x Center'
        self.cbox_xc  = wx.CheckBox(self, wx.NewId(), name)
        self.float_xc = FloatControl(self, wx.NewId())
        self.float_xc.SetValue(det.xc)
        self.float_xc.SetDelta(0.5*det.pixelPitch)

        name = 'y Center'
        self.cbox_yc  = wx.CheckBox(self, wx.NewId(), name)
        self.float_yc = FloatControl(self, wx.NewId())
        self.float_yc.SetValue(det.yc)
        self.float_yc.SetDelta(0.5*det.pixelPitch)

        name = 'Distance'
        self.cbox_D  = wx.CheckBox(self, wx.NewId(), name)
        self.float_D = FloatControl(self, wx.NewId())
        self.float_D.SetValue(det.workDist)
        self.float_D.SetDelta(10*det.pixelPitch)

        name = 'x Tilt'
        self.cbox_xt  = wx.CheckBox(self, wx.NewId(), name)
        self.float_xt = FloatControl(self, wx.NewId())
        self.float_xt.SetValue(det.xTilt)

        name = 'y Tilt'
        self.cbox_yt  = wx.CheckBox(self, wx.NewId(), name)
        self.float_yt = FloatControl(self, wx.NewId())
        self.float_yt.SetValue(det.yTilt)

        name = 'z Tilt'
        self.cbox_zt  = wx.CheckBox(self, wx.NewId(), name)
        self.float_zt = FloatControl(self, wx.NewId())
        self.float_zt.SetValue(det.yTilt)
        #
        #  Distortion parameters
        #
        #  *) NOTE THAT THESE ARE SPECIFIC FOR THE GE
        #  *) must break these out into a subpanel, as the
        #     number (if any at all) will change for each
        #     detector type.
        name = 'p0'
        self.cbox_d1  = wx.CheckBox(self, wx.NewId(), name)
        self.float_d1 = FloatControl(self, wx.NewId())
        self.float_d1.SetValue(det.dparms[0])

        name = 'p1'
        self.cbox_d2  = wx.CheckBox(self, wx.NewId(), name)
        self.float_d2 = FloatControl(self, wx.NewId())
        self.float_d2.SetValue(det.dparms[1])

        name = 'p2'
        self.cbox_d3  = wx.CheckBox(self, wx.NewId(), name)
        self.float_d3 = FloatControl(self, wx.NewId())
        self.float_d3.SetValue(det.dparms[2])

        name = 'n0'
        self.cbox_d4  = wx.CheckBox(self, wx.NewId(), name)
        self.float_d4 = FloatControl(self, wx.NewId())
        self.float_d4.SetValue(det.dparms[3])

        name = 'n1'
        self.cbox_d5  = wx.CheckBox(self, wx.NewId(), name)
        self.float_d5 = FloatControl(self, wx.NewId())
        self.float_d5.SetValue(det.dparms[4])

        name = 'n2'
        self.cbox_d6  = wx.CheckBox(self, wx.NewId(), name)
        self.float_d6 = FloatControl(self, wx.NewId())
        self.float_d6.SetValue(det.dparms[5])
        #
        #  Fitting method
        #
        self.fitLabelSizer = makeTitleBar(self, 'Fitting Method',
                                          color=WP.TITLEBAR_BG_COLOR_PANEL1)
        self.fitDir_rb = wx.RadioButton(self, wx.NewId(), 'Direct Fit',
                                       style=wx.RB_GROUP)
        self.fitBin_rb = wx.RadioButton(self, wx.NewId(), 'Binned Fit')
        #
        #  III. Caking
        #
        self.numEta_lab = wx.StaticText(self, wx.NewId(),
                                        'Azimuthal bins',
                                         style=wx.ALIGN_RIGHT)
        self.numRho_lab = wx.StaticText(self, wx.NewId(),
                                        'Radial bins per ring',
                                         style=wx.ALIGN_RIGHT)
        self.numEta_spn = wx.SpinCtrl(self, wx.NewId(), min=12, initial=36)
        self.numRho_spn = wx.SpinCtrl(self, wx.NewId(), min=10, initial=20)
        #
        #  Fit button with options (at some point)
        #
        self.runFit_but  = wx.Button(self, wx.NewId(), 'Run Fit')

        return

    def __makeBindings(self):
        """Bind interactors"""
        # calibrant section
        self.Bind(wx.EVT_CHOICE,     self.OnMatChoice,   self.mats_cho)

        self.Bind(wx.EVT_SPINCTRL,   self.OnNumRho,  self.numRho_spn)
        self.Bind(wx.EVT_SPINCTRL,   self.OnNumEta,  self.numEta_spn)

        # detector section
        self.Bind(EVT_FLOAT_CTRL, self.OnFloatXC, self.float_xc)
        self.Bind(EVT_FLOAT_CTRL, self.OnFloatYC, self.float_yc)
        self.Bind(EVT_FLOAT_CTRL, self.OnFloatD,  self.float_D)
        self.Bind(EVT_FLOAT_CTRL, self.OnFloatXT, self.float_xt)
        self.Bind(EVT_FLOAT_CTRL, self.OnFloatYT, self.float_yt)
        self.Bind(EVT_FLOAT_CTRL, self.OnFloatZT, self.float_zt)

        self.Bind(EVT_FLOAT_CTRL, self.OnFloatd1, self.float_d1)
        self.Bind(EVT_FLOAT_CTRL, self.OnFloatd2, self.float_d2)
        self.Bind(EVT_FLOAT_CTRL, self.OnFloatd3, self.float_d3)
        self.Bind(EVT_FLOAT_CTRL, self.OnFloatd4, self.float_d4)
        self.Bind(EVT_FLOAT_CTRL, self.OnFloatd5, self.float_d5)
        self.Bind(EVT_FLOAT_CTRL, self.OnFloatd6, self.float_d6)

        # checkboxes
        self.Bind(wx.EVT_CHECKBOX,   self.OnCheck_xc, self.cbox_xc)
        self.Bind(wx.EVT_CHECKBOX,   self.OnCheck_yc, self.cbox_yc)
        self.Bind(wx.EVT_CHECKBOX,   self.OnCheck_D,  self.cbox_D)
        self.Bind(wx.EVT_CHECKBOX,   self.OnCheck_xt, self.cbox_xt)
        self.Bind(wx.EVT_CHECKBOX,   self.OnCheck_yt, self.cbox_yt)
        self.Bind(wx.EVT_CHECKBOX,   self.OnCheck_zt, self.cbox_zt)

        self.Bind(wx.EVT_CHECKBOX,   self.OnCheck_d1, self.cbox_d1)
        self.Bind(wx.EVT_CHECKBOX,   self.OnCheck_d2, self.cbox_d2)
        self.Bind(wx.EVT_CHECKBOX,   self.OnCheck_d3, self.cbox_d3)
        self.Bind(wx.EVT_CHECKBOX,   self.OnCheck_d4, self.cbox_d4)
        self.Bind(wx.EVT_CHECKBOX,   self.OnCheck_d5, self.cbox_d5)
        self.Bind(wx.EVT_CHECKBOX,   self.OnCheck_d6, self.cbox_d6)

        # fitting section
        self.Bind(wx.EVT_RADIOBUTTON, self.OnFitDir, self.fitDir_rb)
        self.Bind(wx.EVT_RADIOBUTTON, self.OnFitBin, self.fitBin_rb)

        self.Bind(wx.EVT_BUTTON, self.OnRunFit, self.runFit_but)

        return

    def __makeSizers(self):
	"""Lay out the interactors"""
        #
        #  Material sizer
        #
	matSizer = wx.BoxSizer(wx.HORIZONTAL)
        matSizer.Add(self.mats_lab, 0, wx.ALIGN_RIGHT)
        matSizer.Add(self.mats_cho, 0, wx.ALIGN_LEFT|wx.LEFT, 5)
        #
        #  Geometry sizer
        #
        nrow = 10; ncol = 2; padx = 5; pady = 5
        self.geoSizer = wx.FlexGridSizer(nrow, ncol, padx, pady)
        self.geoSizer.AddGrowableCol(0, 1)
        self.geoSizer.AddGrowableCol(1, 1)
        #  * x-center
        self.geoSizer.Add(self.cbox_xc,  1, wx.EXPAND)
        self.geoSizer.Add(self.float_xc, 1, wx.EXPAND)
        #  * y-center
        self.geoSizer.Add(self.cbox_yc,  1, wx.EXPAND)
        self.geoSizer.Add(self.float_yc, 1, wx.EXPAND)
        #  * distance
        self.geoSizer.Add(self.cbox_D,  1, wx.EXPAND)
        self.geoSizer.Add(self.float_D, 1, wx.EXPAND)
        #  * x-tilt
        self.geoSizer.Add(self.cbox_xt,  1, wx.EXPAND)
        self.geoSizer.Add(self.float_xt, 1, wx.EXPAND)
        #  * y-tilt
        self.geoSizer.Add(self.cbox_yt,  1, wx.EXPAND)
        self.geoSizer.Add(self.float_yt, 1, wx.EXPAND)
        #  * z-tilt
        self.geoSizer.Add(self.cbox_zt,  1, wx.EXPAND)
        self.geoSizer.Add(self.float_zt, 1, wx.EXPAND)

        #  *** distortion parameters ***
        self.geoSizer.Add( self.cbox_d1,  1, wx.EXPAND)
        self.geoSizer.Add(self.float_d1, 1, wx.EXPAND)
        #
        self.geoSizer.Add( self.cbox_d2,  1, wx.EXPAND)
        self.geoSizer.Add(self.float_d2, 1, wx.EXPAND)
        #
        self.geoSizer.Add( self.cbox_d3,  1, wx.EXPAND)
        self.geoSizer.Add(self.float_d3, 1, wx.EXPAND)
        #
        self.geoSizer.Add( self.cbox_d4,  1, wx.EXPAND)
        self.geoSizer.Add(self.float_d4, 1, wx.EXPAND)
        #
        self.geoSizer.Add( self.cbox_d5,  1, wx.EXPAND)
        self.geoSizer.Add(self.float_d5, 1, wx.EXPAND)
        #
        self.geoSizer.Add( self.cbox_d6,  1, wx.EXPAND)
        self.geoSizer.Add(self.float_d6, 1, wx.EXPAND)
        #
        #  Radio buttons
        #
        self.rbSizer = wx.BoxSizer(wx.VERTICAL)
        self.rbSizer.Add(self.fitDir_rb, 0)
        self.rbSizer.Add(self.fitBin_rb, 0)
        #
        #  Caking sizer
        #
        nrow = 0; ncol = 2; padx = pady = 5
	self.cakeSizer = wx.FlexGridSizer(nrow, ncol, padx, pady)
        #
	self.cakeSizer.AddGrowableCol(0, 1)
        self.cakeSizer.Add(self.numEta_lab, 1, wx.ALIGN_RIGHT)
        self.cakeSizer.Add(self.numEta_spn, 0, wx.ALIGN_LEFT)
        self.cakeSizer.Add(self.numRho_lab, 1, wx.ALIGN_RIGHT)
        self.cakeSizer.Add(self.numRho_spn, 0, wx.ALIGN_LEFT)
        #
        #  Fit sizer
        #
        self.fitSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.fitSizer.Add(self.rbSizer,   1, wx.ALIGN_LEFT)
        self.fitSizer.Add(self.cakeSizer, 0)
        #
        #  Main Sizer: two subpanels
        #
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.tbarSizer,  0, wx.TOP | wx.EXPAND, 5)
        self.sizer.Add(matSizer,  0, wx.TOP | wx.EXPAND, 5)
        self.sizer.Add(self.ring_pan,   0, wx.TOP | wx.EXPAND, 5)
        self.sizer.Add(self.detLabSizer, 1, wx.TOP | wx.EXPAND, 5)
        self.sizer.Add(self.geoSizer,   0, wx.TOP | wx.EXPAND, 5)
        self.sizer.Add(self.fitLabelSizer, 1, wx.TOP | wx.EXPAND, 5)
        self.sizer.Add(self.fitSizer,   0, wx.TOP | wx.EXPAND, 5)
        AddSpacer(self, self.sizer, WP.BG_COLOR_TITLEBAR_PANEL1)
        self.sizer.Add(self.runFit_but, 0, wx.ALIGN_RIGHT|wx.RIGHT, 5)

	return

    def __initExclusions(self):
        """Nothing yet"""
        return

    def __showCbox(self, cb, fl, val):
        """Set status on checkbox"""
        cb.SetValue(val)
        fl.Enable(val)

        return

    def __showGeoParams(self):
        """Display current values of geometric parameters"""
        det = wx.GetApp().ws.detector

        self.float_xc.SetValue(det.xc)
        self.float_yc.SetValue(det.yc)
        self.float_D.SetValue(det.workDist)
        self.float_xt.SetValue(det.xTilt)
        self.float_yt.SetValue(det.yTilt)
        self.float_zt.SetValue(det.zTilt)

        self.float_d1.SetValue(det.dparms[0])
        self.float_d2.SetValue(det.dparms[1])
        self.float_d3.SetValue(det.dparms[2])
        self.float_d4.SetValue(det.dparms[3])
        self.float_d5.SetValue(det.dparms[4])
        self.float_d6.SetValue(det.dparms[5])

        self.__showCbox(self.cbox_xc, self.float_xc, det.refineFlags[0])
        self.__showCbox(self.cbox_yc, self.float_yc, det.refineFlags[1])
        self.__showCbox(self.cbox_D,  self.float_D,  det.refineFlags[2])
        self.__showCbox(self.cbox_xt, self.float_xt, det.refineFlags[3])
        self.__showCbox(self.cbox_yt, self.float_yt, det.refineFlags[4])
        self.__showCbox(self.cbox_zt, self.float_zt, det.refineFlags[5])
        self.__showCbox(self.cbox_d1, self.float_d1, det.refineFlags[6])
        self.__showCbox(self.cbox_d2, self.float_d2, det.refineFlags[7])
        self.__showCbox(self.cbox_d3, self.float_d3, det.refineFlags[8])
        self.__showCbox(self.cbox_d4, self.float_d4, det.refineFlags[9])
        self.__showCbox(self.cbox_d5, self.float_d5, det.refineFlags[10])
        self.__showCbox(self.cbox_d6, self.float_d6, det.refineFlags[11])

        return

    def __showFitOpts(self):
        """Set option interactors"""
        exp = wx.GetApp().ws
        cin = exp.calInput

        self.numRho_spn.SetValue(cin.numRho)
        self.numEta_spn.SetValue(cin.numEta)

        if cin.fitType == FitModes.DIRECT:
            self.fitDir_rb.SetValue(True)
            self.numRho_spn.Disable()
            self.numEta_spn.Disable()
        else:
            self.fitBin_rb.SetValue(True)
            self.numRho_spn.Enable()
            self.numEta_spn.Enable()
            pass

        return
    #
    # ============================== API
    #
    #                     ========== *** Access Methods
    #
    def updateFromExp(self):
        """Update all subwindows"""
        app = wx.GetApp()
        exp = app.ws
        mat = exp.activeMaterial

        ResetChoice(self.mats_cho, exp.matNames, mat.name)

        self.__showGeoParams()
        self.ring_pan.updateFromExp()
        self.__showFitOpts()

        app.getCanvas().update()

        return
    #
    #                     ========== *** Event Callbacks
    #
    def OnMatChoice(self, e):
        """Select new material"""
        exp = wx.GetApp().ws

        sel = self.mats_cho.GetSelection()
        if sel >= 0:
            exp.activeMaterial = sel
            self.updateFromExp()
            pass

        return

    def OnRunFit(self, evt):
        """Callback for runFit_but"""
        #
        #  Get workspace.
        #
        exp = wx.GetApp().ws

        try:
            action = {
                'exec': exp.calibrate,
                'args': (),
                'kwargs': dict()
                }
            logwin = logWindow(self, wx.NewId(), action, 'Fitting Log')
            logwin.ShowModal()

        except Exception as e:
            wx.MessageBox(str(e))
            pass

        self.Refresh()
        self.updateFromExp()
        #
        return

    def OnFitDir(self, e):
        """Direct fit chosen"""
        exp = wx.GetApp().ws
        exp.calInput.fitType = FitModes.DIRECT

        self.numRho_spn.Disable()
        self.numEta_spn.Disable()
        return

    def OnFitBin(self, e):
        """Binned fit chosen"""
        exp = wx.GetApp().ws
        exp.calInput.fitType = FitModes.MULTIRING

        self.numRho_spn.Enable()
        self.numEta_spn.Enable()
        return
    #
    #  Detector Parameters
    #
    def OnFloatXC(self, evt):
        """Callback for float_xc choice"""
        try:
            a = wx.GetApp()
            a.ws.detector.xc = evt.floatValue
            a.getCanvas().update()

        except Exception as e:
            msg = 'Failed to set x-center value: \n%s' % str(e)
            wx.MessageBox(msg)
            pass

        return

    def OnFloatYC(self, evt):
        """Callback for float_yc choice"""
        try:
            a = wx.GetApp()
            a.ws.detector.yc = evt.floatValue
            a.getCanvas().update()

        except Exception as e:
            msg = 'Failed to set y-center value: \n%s' % str(e)
            wx.MessageBox(msg)
            pass

        return

    def OnFloatD(self, evt):
        """Callback for float_d choice"""
        try:
            a = wx.GetApp()
            a.ws.detector.workDist = evt.floatValue
            a.getCanvas().update()

        except Exception as e:
            msg = 'Failed to set detector distance value: \n%s' % str(e)
            wx.MessageBox(msg)
            pass

        return

    def OnFloatXT(self, evt):
        """Callback for float_xt choice"""
        try:
            a = wx.GetApp()
            a.ws.detector.xTilt = evt.floatValue
            a.getCanvas().update()

        except Exception as e:
            msg = 'Failed to set x-tilt value: \n%s' % str(e)
            wx.MessageBox(msg)
            pass

        return

    def OnFloatYT(self, evt):
        """Callback for float_yt choice"""
        try:
            a = wx.GetApp()
            a.ws.detector.yTilt = evt.floatValue
            a.getCanvas().update()

        except Exception as e:
            msg = 'Failed to set y-tilt value: \n%s' % str(e)
            wx.MessageBox(msg)
            pass

        return

    def OnFloatZT(self, evt):
        """Callback for float_yt choice"""
        try:
            a = wx.GetApp()
            a.ws.detector.zTilt = evt.floatValue
            a.getCanvas().update()

        except Exception as e:
            msg = 'Failed to set z-tilt value: \n%s' % str(e)
            wx.MessageBox(msg)
            pass

        return

    def OnFloatd1(self, evt):
        """Callback for float_d1 choice"""
        try:
            a = wx.GetApp()
            a.ws.detector.dparms[0] = evt.floatValue
            a.getCanvas().update()

        except Exception as e:
            msg = 'Failed to set distortion parameter 1: \n%s' % str(e)
            wx.MessageBox(msg)
            pass

        return

    def OnFloatd2(self, evt):
        """Callback for float_d2 choice"""
        try:
            a = wx.GetApp()
            a.ws.detector.dparms[1] = evt.floatValue
            a.getCanvas().update()

        except Exception as e:
            msg = 'Failed to set distortion parameter 1: \n%s' % str(e)
            wx.MessageBox(msg)
            pass

        return

    def OnFloatd3(self, evt):
        """Callback for float_d3 choice"""
        try:
            a = wx.GetApp()
            a.ws.detector.dparms[2] = evt.floatValue
            a.getCanvas().update()

        except Exception as e:
            msg = 'Failed to set distortion parameter 1: \n%s' % str(e)
            wx.MessageBox(msg)
            pass

        return

    def OnFloatd4(self, evt):
        """Callback for float_d4 choice"""
        try:
            a = wx.GetApp()
            a.ws.detector.dparms[3] = evt.floatValue
            a.getCanvas().update()

        except Exception as e:
            msg = 'Failed to set distortion parameter 1: \n%s' % str(e)
            wx.MessageBox(msg)
            pass

        return

    def OnFloatd5(self, evt):
        """Callback for float_d5 choice"""
        try:
            a = wx.GetApp()
            a.ws.detector.dparms[4] = evt.floatValue
            a.getCanvas().update()

        except Exception as e:
            msg = 'Failed to set distortion parameter 5: \n%s' % str(e)
            wx.MessageBox(msg)
            pass

        return

    def OnFloatd6(self, evt):
        """Callback for float_d6 choice"""
        try:
            a = wx.GetApp()
            a.ws.detector.dparms[5] = evt.floatValue
            a.getCanvas().update()

        except Exception as e:
            msg = 'Failed to set distortion parameter 6: \n%s' % str(e)
            wx.MessageBox(msg)
            pass

        return


    def OnNumRho(self, e):
        """Number of rho bins has changed"""
        exp = wx.GetApp().ws
        exp.calInput.numRho = self.numRho_spn.GetValue()
        print 'changing numrho:', exp.calInput.numRho

        return

    def OnNumEta(self, e):
        """Number of eta bins has changed"""
        exp = wx.GetApp().ws
        exp.calInput.numEta = self.numEta_spn.GetValue()

        return

    #
    #  Checkboxes
    #
    def OnCheck_xc(self, e):
        """xc box is checked"""
        ind  = 0
        fc   = self.float_xc

        exp  = wx.GetApp().ws
        stat = e.IsChecked()

        exp.refineFlags[ind] = stat
        fc.Enable(stat)

        #print 'xc:  ', stat

        return

    def OnCheck_yc(self, e):
        """xc box is checked"""
        fc   = self.float_yc
        ind  = 1

        exp  = wx.GetApp().ws
        stat = e.IsChecked()

        exp.refineFlags[ind] = stat
        fc.Enable(stat)

        return

    def OnCheck_D(self, e):
        """xc box is checked"""
        fc   = self.float_D
        ind  = 2

        exp  = wx.GetApp().ws
        stat = e.IsChecked()

        exp.refineFlags[ind] = stat
        fc.Enable(stat)

        return

    def OnCheck_xt(self, e):
        """xc box is checked"""
        ind  = 3
        fc   = self.float_xt

        exp  = wx.GetApp().ws
        stat = e.IsChecked()

        exp.refineFlags[ind] = stat
        fc.Enable(stat)

        return

    def OnCheck_yt(self, e):
        """xc box is checked"""
        ind  = 4
        fc   = self.float_yt

        exp  = wx.GetApp().ws
        stat = e.IsChecked()

        exp.refineFlags[ind] = stat
        fc.Enable(stat)

        return

    def OnCheck_zt(self, e):
        """xc box is checked"""
        ind  = 5
        fc   = self.float_zt

        exp  = wx.GetApp().ws
        stat = e.IsChecked()

        exp.refineFlags[ind] = stat
        fc.Enable(stat)

        return

    def OnCheck_d1(self, e):
        """xc box is checked"""
        ind  = 6
        fc   = self.float_d1

        exp  = wx.GetApp().ws
        stat = e.IsChecked()

        exp.refineFlags[ind] = stat
        fc.Enable(stat)

        return

    def OnCheck_d2(self, e):
        """xc box is checked"""
        ind  = 7
        fc   = self.float_d2

        exp  = wx.GetApp().ws
        stat = e.IsChecked()

        exp.refineFlags[ind] = stat
        fc.Enable(stat)

        return

    def OnCheck_d3(self, e):
        """xc box is checked"""
        ind  = 8
        fc   = self.float_d3

        exp  = wx.GetApp().ws
        stat = e.IsChecked()

        exp.refineFlags[ind] = stat
        fc.Enable(stat)

        return

    def OnCheck_d4(self, e):
        """xc box is checked"""
        ind  = 9
        fc   = self.float_d4

        exp  = wx.GetApp().ws
        stat = e.IsChecked()

        exp.refineFlags[ind] = stat
        fc.Enable(stat)

        return

    def OnCheck_d5(self, e):
        """xc box is checked"""
        ind  = 10
        fc   = self.float_d5

        exp  = wx.GetApp().ws
        stat = e.IsChecked()

        exp.refineFlags[ind] = stat
        fc.Enable(stat)

        return

    def OnCheck_d6(self, e):
        """xc box is checked"""
        ind  = 11
        fc   = self.float_d6

        exp  = wx.GetApp().ws
        stat = e.IsChecked()

        exp.refineFlags[ind] = stat
        fc.Enable(stat)

        return
    pass # end class
#
# -----------------------------------------------END CLASS:  plotPanel

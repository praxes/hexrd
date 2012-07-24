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
"""Module to present caking options
"""
import wx

import numpy

from hexrd.xrd import detector
from hexrd.xrd.Experiment  import PolarRebinOpts as prOpts
from hexrd.xrd.xrdUtils    import CollapseOmeEta

from hexrd.wx.canvasUtilities import *
from hexrd.wx.guiConfig    import WindowParameters as WP
from hexrd.wx.guiUtilities import makeTitleBar
from hexrd.wx.ringSubPanel import ringPanel
from hexrd.wx.selectHKLs   import selectHKLsDialog as hklsDlg
from hexrd.wx.LogWindows   import logWindow
from hexrd.wx.FloatControl import *
from hexrd.wx.cakingCanvas import cakeDisplay
#
#  Module Data
#
piBy180 = numpy.pi / 180.0
#
# ---------------------------------------------------CLASS:  cakingPanel
#
class cakingPanel(wx.Panel):
    """cakingPanel """
    def __init__(self, parent, id, **kwargs):
	"""Constructor for cakingPanel."""
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
        self.update(showImg=True)
	#
	return
    #
    # ============================== Internal Methods
    #
    def __makeObjects(self):
        """Add interactors"""

        self.tbarSizer = makeTitleBar(self, 'Polar Rebinning')
        #
        self.std_pan = standardOptsPanel(self, wx.NewId())
        self.mrb_pan = multiringOptsPanel(self, wx.NewId())
        self.sph_pan = sphericalOptsPanel(self, wx.NewId())
        #
        #  Method
        #
        self.method_cho = wx.Choice(self, wx.NewId(), choices=prOpts.cakeMethods)
        self.method_cho.SetSelection(1)
        #
        #  Run
        #
        self.run_but   = wx.Button(self, wx.NewId(), 'Run Polar Rebin')
        #
        #  Canvas for figures moved to separate window:  see cakingCanvas.py
        #
        return

    def __makeBindings(self):
        """Bind interactors"""
        #
	self.Bind(wx.EVT_CHOICE, self.OnChooseMethod, self.method_cho)
        self.Bind(wx.EVT_BUTTON, self.OnRun,    self.run_but)

        return

    def __makeSizers(self):
	"""Lay out the interactors"""
	#
        #  Control sizer
        #
        self.ctrlSizer = wx.BoxSizer(wx.VERTICAL)
        self.ctrlSizer.Add(self.method_cho, 0, wx.ALIGN_RIGHT|wx.BOTTOM, 15)
        self.ctrlSizer.Add(self.run_but,    0, wx.ALIGN_RIGHT|wx.BOTTOM, 5)
        #
        #  Options sizer
        #
        self.optSizer = wx.BoxSizer(wx.HORIZONTAL)
        #
        self.optSizer.Add(self.std_pan, 0, wx.ALIGN_RIGHT)
        self.optSizer.Add(self.mrb_pan, 0, wx.ALIGN_RIGHT)
        self.optSizer.Add(self.sph_pan, 0, wx.ALIGN_RIGHT)
        #
        self.optSizer.Show(self.mrb_pan, False)
        self.optSizer.Show(self.sph_pan, False)
        #
        #  Mode sizer
        #
        sep = wx.StaticLine(self, -1, style=wx.LI_VERTICAL, size=(5,5))
        sep.SetBackgroundColour('black')
        #
        self.modeSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.modeSizer.Add(self.ctrlSizer, 0, wx.ALIGN_RIGHT)
        self.modeSizer.Add(sep, 0,
                           wx.EXPAND|wx.ALIGN_CENTER|wx.LEFT|wx.RIGHT, 10)
        self.modeSizer.Add(self.optSizer,  1, wx.EXPAND|wx.ALIGN_CENTER)
        #
        #  Main Sizer
        #
	self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.sizer.Add(self.tbarSizer, 0, wx.EXPAND|wx.ALIGN_CENTER|wx.BOTTOM, 5)
	self.sizer.Add(self.modeSizer, 0, wx.EXPAND|wx.ALIGN_CENTER)
        #
        self.Layout()
        #  Initialize method choice

        self.OnChooseMethod(None)

	return

    def __cake_img(self):
        """Execute rebinning for IMG case
"""
        exp = wx.GetApp().ws
        det = exp.detector

        action = {
            'exec'  : self.__makeImgInfo,
            'args'  : (),
            'kwargs': dict()
            }
        logwin = logWindow(self, wx.NewId(), action, 'Standard Polar Rebinning')
        logwin.ShowModal()

        return
        #
        #  ==============================
        #  Show rebin output on new window
        #  ==============================
        #
        #  Now draw
        #
        self.axes   = self.figure.gca()
        self.axes.set_aspect('equal')

        intensity = self.img_info['intensity']

        self.axes.images = []
        # show new image
        self.axes.imshow(intensity, origin='upper',
                         interpolation='nearest',
                         cmap=self.cmPanel.cmap,
                         vmin=self.cmPanel.cmin_val,
                         vmax=self.cmPanel.cmax_val,
                         aspect='auto')
        #self.axes.set_aspect('equal')
        self.axes.set_autoscale_on(False)

        self.canvas.draw()

        #
        pass

    def __cake_rng(self):
        """Multiring cake method"""
        exp = wx.GetApp().ws
        det = exp.detector

        action = {
            'exec'  : self.__makeMRB,
            'args'  : (),
            'kwargs': dict()
            }
        logwin = logWindow(self, wx.NewId(), action, 'Multiring Binning')
        logwin.ShowModal()

        # ====================

        cCan = cakeDisplay(self, wx.NewId(), prOpts.CAKE_RNG, self.mrb)

        return

    def __makeImgInfo(self, log=None):
        """routine for calling polarRebin"""
        # GET data from interactors

        exp = wx.GetApp().ws
        det = exp.detector

        etaMin = self.std_pan.emin_spn.GetValue()*piBy180
        etaMax = self.std_pan.emax_spn.GetValue()*piBy180
        rhoMin = self.std_pan.rmin_spn.GetValue()
        rhoMax = self.std_pan.rmax_spn.GetValue()
        numEta = self.std_pan.enum_spn.GetValue()
        numRho = self.std_pan.rnum_spn.GetValue()

        correct = self.std_pan.corr_cbox.GetValue()

        kwa = {
            'etaRange' : [etaMin, etaMax],
            'numEta'   : numEta,
            'rhoRange' : [rhoMin, rhoMax],
            'numRho'   : numRho,
            'corrected': correct
            }

        self.img_info = det.polarRebin(exp.activeImage, log=log, **kwa)

        log.write('\ndone', done=True)

        return

    def __makeMRB(self, log=None):
        """routine for calling MultiRingBinned"""
        exp = wx.GetApp().ws
        det = exp.detector
        pdat = exp.activeMaterial.planeData
        img  = exp.activeImage

        nrho = self.mrb_pan.numRho_spn.GetValue()
        neta = self.mrb_pan.numEta_spn.GetValue()
        etaMin = self.mrb_pan.emin*piBy180
        etaMax = self.mrb_pan.emax*piBy180

        cakeArgs = {'verbose'  : True,
                    'numEta'   : neta,
                    'corrected': True,
                    'etaRange' : [etaMin, etaMax]}

        self.mrb = detector.MultiRingBinned(det, pdat, img,
                                            targetNRho=nrho,
                                            polarRebinKWArgs=cakeArgs,
                                            log=log)

        log.write('\ndone', done=True)

        return

    def __cake_sph(self):
        """Build eta-omega plots"""
        exp = wx.GetApp().ws

        reader = exp.activeReader.makeReader() # need to load actual reader
        pdata  = exp.activeMaterial.planeData
        hklIDs = pdata.getHKLID(pdata.hkls.T.tolist())
        det    = exp.detector

        nlump  = self.sph_pan.lump_spn.GetValue()
        nbins  = self.sph_pan.bins_spn.GetValue()
        thresh = self.sph_pan.thresh_spn.GetValue()
        kwargs = dict(
            nframesLump=nlump,
            nEtaBins=nbins,
            threshold=thresh,
            debug=True)

        omeEta = CollapseOmeEta(reader, pdata, hklIDs, det, **kwargs)

        print omeEta
        print dir(omeEta)

        # Show some pole figures
        style = {'marker':'d',    'ls':'None', 'mec':'y',
                 'mfc'   :'None', 'ms':10.,    'mew':2}

        res    = 100
        cmap   = self.cmPanel.cmap
        pt     = numpy.r_[1.0, 0.0, 0.0]
        pList  = [(pt, style)] # to do
        iHKL   = 0
        kwargs = dict(
            pfig=res,
            cmap=cmap,
            iData=iHKL,
            pointLists=pList,
            rangeVV=(0., 150.),
            drawColorbar=False,
            pfigDisplayKWArgs={'doY90Rot':False}
            )
        pfig = omeEta.display(**kwargs)

        print pfig, dir(pfig)

        return
    #
    # ============================== API
    #
    #                     ========== *** Access Methods
    #
    def update(self, showImg=False):
        """Update figure"""
        #  off for now ...

        return

    #
    #                     ========== *** Event Callbacks
    #
    def OnChooseMethod(self, e):
        """Binning method has changed"""
        #self.sizer.Show(self.canvas, False)

        s = self.method_cho.GetStringSelection()
        print 'rebin method:  ', s
        if s == prOpts.CAKE_IMG:
            self.optSizer.Show(self.mrb_pan, False)
            self.optSizer.Show(self.sph_pan, False)
            self.optSizer.Show(self.std_pan)

        elif s == prOpts.CAKE_RNG:
            self.optSizer.Show(self.std_pan, False)
            self.optSizer.Show(self.sph_pan, False)
            self.optSizer.Show(self.mrb_pan)

        elif s == prOpts.CAKE_SPH:
            self.optSizer.Show(self.std_pan, False)
            self.optSizer.Show(self.mrb_pan, False)
            self.optSizer.Show(self.sph_pan)
            pass

        self.sizer.Layout()

        self.update(showImg=True)

        #self.sizer.Show(self.canvas, True)

        self.sizer.Layout()

        return

    def OnRun(self, e):
        """Run caking"""

        s = self.method_cho.GetStringSelection()

        if s == prOpts.CAKE_IMG:
            self.__cake_img()
        elif s == prOpts.CAKE_RNG:
            self.__cake_rng()
        elif s == prOpts.CAKE_SPH:
            self.__cake_sph()
            pass

        return

    pass # end class
#
# -----------------------------------------------END CLASS:  cakingPanel
# ---------------------------------------------------CLASS:  cakingDialog
#
class cakingDialog(wx.Dialog):
    """cakingDialog """
    def __init__(self, parent, id, **kwargs):
	"""Constructor for cakingDialog"""
	#
	wx.Dialog.__init__(self, parent, id, title="Polar Rebinning",
                           style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
	#
	#  Data Objects.
	#

	#
	#  Windows.
	#
	self.titlebar = wx.StaticText(self, -1, 'cakingDialog',
				      style=wx.ALIGN_CENTER|wx.SIMPLE_BORDER)
        self.dataPanel = cakingPanel(self, wx.NewId())
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
        self.sizer.Show(self.titlebar, False)

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
# -----------------------------------------------END CLASS:  cakingDialog
# ---------------------------------------------------CLASS:  standardOptsPanel
#
class standardOptsPanel(wx.Panel):
    """standardOptsPanel """
    def __init__(self, parent, id, **kwargs):
	"""Constructor for standardOptsPanel."""
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
        self.update()
	#
	return
    #
    # ============================== Internal Methods
    #
    def __makeObjects(self):
        """Add interactors"""

        self.tbarSizer = makeTitleBar(self, 'Options for Standard Rebin',
                                      color=WP.BG_COLOR_TITLEBAR_PANEL1)
        #
        #  Labels
        #
        self.min_lab = wx.StaticText(self, wx.NewId(), 'min', style=wx.ALIGN_CENTER)
        self.max_lab = wx.StaticText(self, wx.NewId(), 'max', style=wx.ALIGN_CENTER)
        self.num_lab = wx.StaticText(self, wx.NewId(), 'num', style=wx.ALIGN_CENTER)
        self.rho_lab = wx.StaticText(self, wx.NewId(), 'rho', style=wx.ALIGN_CENTER)
        self.eta_lab = wx.StaticText(self, wx.NewId(), 'eta', style=wx.ALIGN_CENTER)
        #
        #  Rho
        #
        self.rmin_spn  = wx.SpinCtrl(self, wx.NewId(), min=1, max=10000, value=str(100))
        self.rmax_spn  = wx.SpinCtrl(self, wx.NewId(), min=1, max=10000, value=str(1000))
        self.rnum_spn  = wx.SpinCtrl(self, wx.NewId(), min=1, max=10000, value=str(500))
        #
        #  Eta
        #
        self.emin_spn  = wx.SpinCtrl(self, wx.NewId(), min=-360, max=360, value=str(0))
        self.emax_spn  = wx.SpinCtrl(self, wx.NewId(), min=-360, max=360, value=str(360))
        self.enum_spn  = wx.SpinCtrl(self, wx.NewId(), min=1,    max=360, value=str(36))
        #
        #  Other options
        #
        self.corr_cbox = wx.CheckBox(self, wx.NewId(), 'corrected')
        self.npdv_lab  = wx.StaticText(self, wx.NewId(), 'pixel divisions',
                                      style=wx.ALIGN_CENTER)
        self.npdv_spn  = wx.SpinCtrl(self, wx.NewId(), min=1, max=10, initial=1)
        self.frame_lab = wx.StaticText(self, wx.NewId(), 'frame',
                                       style=wx.ALIGN_CENTER)
        self.frame_cho = wx.Choice(self, wx.NewId(), choices=['frame 1'])

        return

    def __makeBindings(self):
        """Bind interactors"""#  SpinCtrl Events:
        self.Bind(wx.EVT_SPINCTRL, self.OnRminSpn, self.rmin_spn)
        self.Bind(wx.EVT_SPINCTRL, self.OnRminSpn, self.rmax_spn)
        self.Bind(wx.EVT_SPINCTRL, self.OnRminSpn, self.emin_spn)
        self.Bind(wx.EVT_SPINCTRL, self.OnRminSpn, self.emax_spn)
        self.Bind(wx.EVT_SPINCTRL, self.OnRminSpn, self.enum_spn)
        self.Bind(wx.EVT_CHECKBOX, self.OnRminSpn, self.corr_cbox)

        return

    def __makeSizers(self):
	"""Lay out the interactors"""
        #
        #  Rho and Eta subdivisions
        #
        nrow = 0; ncol = 4; padx = pady = 5
	self.mmSizer = wx.FlexGridSizer(nrow, ncol, padx, pady)
        # top row (header)
        self.mmSizer.AddSpacer(1)
        self.mmSizer.Add(self.min_lab, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.mmSizer.Add(self.max_lab, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.mmSizer.Add(self.num_lab, 0, wx.EXPAND|wx.ALIGN_CENTER)
        # rho row
        self.mmSizer.Add(self.rho_lab, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.mmSizer.Add(self.rmin_spn, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.mmSizer.Add(self.rmax_spn, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.mmSizer.Add(self.rnum_spn, 0, wx.EXPAND|wx.ALIGN_CENTER)
        # eta row
        self.mmSizer.Add(self.eta_lab, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.mmSizer.Add(self.emin_spn, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.mmSizer.Add(self.emax_spn, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.mmSizer.Add(self.enum_spn, 0, wx.EXPAND|wx.ALIGN_CENTER)
        #
        #  Other options
        #
        nrow = 0; ncol = 2; padx = pady = 5
	self.optSizer = wx.FlexGridSizer(nrow, ncol, padx, pady)
        #
        self.optSizer.AddSpacer(1)
        self.optSizer.Add(self.corr_cbox, 0, wx.EXPAND|wx.ALIGN_CENTER)
        #
        self.optSizer.Add(self.npdv_spn, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.optSizer.Add(self.npdv_lab, 0, wx.EXPAND|wx.ALIGN_CENTER)
        #
        self.optSizer.Add(self.frame_cho, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.optSizer.Add(self.frame_lab, 0, wx.EXPAND|wx.ALIGN_CENTER)
        #
        self.hsizer = wx.BoxSizer(wx.HORIZONTAL)
        sep = wx.StaticLine(self, -1, style=wx.LI_VERTICAL, size=(5,5))
        sep.SetBackgroundColour('black')
	self.hsizer.Add(self.optSizer,  0, wx.EXPAND|wx.ALIGN_CENTER)
        self.hsizer.Add(sep, 0, wx.EXPAND|wx.ALIGN_CENTER|wx.LEFT|wx.RIGHT, 5)
	self.hsizer.Add(self.mmSizer,   0, wx.EXPAND|wx.ALIGN_CENTER)
        #
	self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.sizer.Add(self.tbarSizer, 0, wx.EXPAND|wx.ALIGN_CENTER)
	self.sizer.Add(self.hsizer,    1, wx.EXPAND|wx.ALIGN_CENTER)
        #
        #  For mac, refresh the windows to display initial value
        #
        self.rmin_spn.Refresh()
        self.Refresh()

	return
    #
    # ============================== API
    #
    #                     ========== *** Access Methods
    #
    def update(self):
        """Update the image """
        app    = wx.GetApp()
        exp    = app.ws
        canvas = app.getCanvas()

        etaMin = self.emin_spn.GetValue()*piBy180
        etaMax = self.emax_spn.GetValue()*piBy180
        rhoMin = self.rmin_spn.GetValue()
        rhoMax = self.rmax_spn.GetValue()
        numEta = self.enum_spn.GetValue()

        correct = self.corr_cbox.GetValue()

        kwa = {
            'etaRange' : [etaMin, etaMax],
            'numEta'   : numEta,
            'rhoRange' : [rhoMin, rhoMax],
            'corrected': correct
            }
        xyRings = exp.detector.getPRBOverlay(kwa)

        opts = {'color': 'r'}
        canvas.clearLines()
        canvas.addXYplot(xyRings, opts=opts)
        return

    #
    #                     ========== *** Event Callbacks
    #
    def OnRminSpn(self, e):
        """Rho min has changed"""
        self.update()

        return
#OnRmaxSpn
#OnEminSpn
#OnEmaxSpn
#OnEnumSpn


    pass # end class
#
# -----------------------------------------------END CLASS:  standardOptsPanel
# ---------------------------------------------------CLASS:  multiringOptsPanel
#
class multiringOptsPanel(wx.Panel):
    """multiringOptsPanel """
    def __init__(self, parent, id, **kwargs):
	"""Constructor for multiringOptsPanel."""
	#
	wx.Panel.__init__(self, parent, id, **kwargs)
	#
        #  Data
        #
        self.emin = 0
        self.emax = 360
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
        self.update()
	#
	return
    #
    # ============================== Internal Methods
    #
    def __makeObjects(self):
        """Add interactors"""

        self.tbarSizer = makeTitleBar(self, 'Options for Multiring Rebin')

        self.ring_pan = ringPanel(self, wx.NewId())

        self.emin_lab = wx.StaticText(self, wx.NewId(),
                                      'Eta min',
                                      style=wx.ALIGN_RIGHT)
        self.emax_lab = wx.StaticText(self, wx.NewId(),
                                      'Eta max',
                                      style=wx.ALIGN_RIGHT)
        self.emin_txt = wx.TextCtrl(self, wx.NewId(),
                                    value='0',
                                    style=wx.RAISED_BORDER | wx.TE_PROCESS_ENTER)
        self.emax_txt = wx.TextCtrl(self, wx.NewId(),
                                    value='360',
                                    style=wx.RAISED_BORDER | wx.TE_PROCESS_ENTER)
        self.numEta_lab = wx.StaticText(self, wx.NewId(),
                                        'Number of Eta Bins',
                                         style=wx.ALIGN_RIGHT)
        self.numRho_lab = wx.StaticText(self, wx.NewId(),
                                        'Rho Bins Per Ring',
                                         style=wx.ALIGN_RIGHT)
        self.numEta_spn = wx.SpinCtrl(self, wx.NewId(), min=1, value=str(36) )
        self.numRho_spn = wx.SpinCtrl(self, wx.NewId(), min=1, value=str(20) )

        return

    def __makeBindings(self):
        """Bind interactors"""
        self.Bind(wx.EVT_TEXT_ENTER, self.OnEMinTxt, self.emin_txt)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnEMaxTxt, self.emax_txt)
        self.Bind(wx.EVT_SPINCTRL, self.OnUpdate, self.numEta_spn)
        self.Bind(wx.EVT_SPINCTRL, self.OnUpdate, self.numRho_spn)

        rp = self.ring_pan

        self.Bind(wx.EVT_TEXT_ENTER, self.OnUpdate, rp.waveAng_txt)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnUpdate, rp.waveKEV_txt)

	self.Bind(wx.EVT_CHOICE,     self.OnUpdate, rp.width_cho)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnUpdate, rp.width_txt)

        self.Bind(wx.EVT_BUTTON,     self.OnUpdate, rp.hkl_but)

        return

    def __makeSizers(self):
	"""Lay out the interactors"""
        #
        #  Caking sizer
        #
        nrow = 0; ncol = 4; padx = pady = 5
	self.cakeSizer = wx.FlexGridSizer(nrow, ncol, padx, pady)
        #
	self.cakeSizer.AddGrowableCol(0, 1)
        self.cakeSizer.Add(self.numEta_lab, 1, wx.ALIGN_RIGHT)
        self.cakeSizer.Add(self.numEta_spn, 0, wx.ALIGN_LEFT)
        self.cakeSizer.Add(self.numRho_lab, 1, wx.ALIGN_RIGHT)
        self.cakeSizer.Add(self.numRho_spn, 0, wx.ALIGN_LEFT)
        #
        self.cakeSizer.Add(self.emin_lab, 1, wx.ALIGN_RIGHT)
        self.cakeSizer.Add(self.emin_txt, 0, wx.ALIGN_LEFT)
        self.cakeSizer.Add(self.emax_lab, 1, wx.ALIGN_RIGHT)
        self.cakeSizer.Add(self.emax_txt, 0, wx.ALIGN_LEFT)

	self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.sizer.Add(self.tbarSizer, 0, wx.EXPAND|wx.ALIGN_CENTER)
	self.sizer.Add(self.cakeSizer,  0, wx.EXPAND|wx.ALIGN_CENTER)
	self.sizer.Add(self.ring_pan,  0, wx.EXPAND|wx.ALIGN_CENTER)

	return
    #
    # ============================== API
    #
    #                     ========== *** Access Methods
    #
    def update(self):
        """Update the image """
        app    = wx.GetApp()
        exp    = app.ws
        pdat   = exp.activeMaterial.planeData
        canvas = app.getCanvas()

        xyRings = exp.detector.getRings(pdat, ranges=True)

        opts = {'color': 'r'}
        canvas.clearLines()
        canvas.addXYplot(xyRings, opts=opts)

        return
    #
    #                     ========== *** Event Callbacks
    #
    def OnUpdate(self, e):
        """Update image"""
        self.update()

        return
    def OnEMinTxt(self, e):
        """Callback for emin_txt choice"""
        self.emin = float(self.emin_txt.GetValue())
        print 'setting Eta min: ', self.emin
        self.update()

        return
    def OnEMaxTxt(self, e):
        """Callback for emin_txt choice"""
        self.emax = float(self.emax_txt.GetValue())
        print 'setting Eta min: ', self.emax
        self.update()

        return

    pass # end class
#
# -----------------------------------------------END CLASS:  multiringOptsPanel
# ---------------------------------------------------CLASS:  sphericalOptsPanel
#
class sphericalOptsPanel(wx.Panel):
    """sphericalOptsPanel """
    def __init__(self, parent, id, **kwargs):
	"""Constructor for sphericalOptsPanel."""
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

        self.tbarSizer = makeTitleBar(self, 'Options for Spherical Rebin')
        #
        #  Integer inputs
        #
        self.lump_lab = wx.StaticText(self, wx.NewId(), 'lump', style=wx.ALIGN_CENTER)
        self.lump_spn  = wx.SpinCtrl(self, wx.NewId(), min=1, max=1000, initial=1)

        self.bins_lab = wx.StaticText(self, wx.NewId(), 'bins', style=wx.ALIGN_CENTER)
        self.bins_spn  = wx.SpinCtrl(self, wx.NewId(), min=1, max=10000, initial=600)

        self.thresh_lab = wx.StaticText(self, wx.NewId(), 'threshold',
                                        style=wx.ALIGN_CENTER)
        self.thresh_spn  = wx.SpinCtrl(self, wx.NewId(), min=1, max=1000, initial=20)
        #
        #  Material and HKLs selector
        #
        exp = wx.GetApp().ws
        self.matl_cho = wx.Choice(self, wx.NewId(), choices=exp.matNames)
        self.matl_cho.SetSelection(0)
        self.read_cho = wx.Choice(self, wx.NewId(), choices=exp.readerNames)
        self.read_cho.SetSelection(0)
        self.hkls_but = wx.Button(self, wx.NewId(), 'Select HKL')
        #
        #  Angle/axis
        #
        name = 'angle'
        self.angle_lab = wx.StaticText(self, wx.NewId(), name, style=wx.ALIGN_CENTER)
        self.angle_flt = FloatControl(self, wx.NewId())
        self.angle_flt.SetValue(1.0)

        name = 'axis x'
        self.axis1_lab = wx.StaticText(self, wx.NewId(), name, style=wx.ALIGN_CENTER)
        self.axis1_flt = FloatControl(self, wx.NewId())
        self.axis1_flt.SetValue(1.0)

        name = 'axis y'
        self.axis2_lab = wx.StaticText(self, wx.NewId(), name, style=wx.ALIGN_CENTER)
        self.axis2_flt = FloatControl(self, wx.NewId())
        self.axis2_flt.SetValue(1.0)

        name = 'axis z'
        self.axis3_lab = wx.StaticText(self, wx.NewId(), name, style=wx.ALIGN_CENTER)
        self.axis3_flt = FloatControl(self, wx.NewId())
        self.axis3_flt.SetValue(1.0)


        return

    def __makeBindings(self):
        """Bind interactors"""
        self.Bind(wx.EVT_CHOICE, self.OnChooseMatl, self.matl_cho)
        self.Bind(wx.EVT_CHOICE, self.OnChooseReader, self.read_cho)
        self.Bind(wx.EVT_BUTTON, self.OnSelectHKLs, self.hkls_but)

        return

    def __makeSizers(self):
	"""Lay out the interactors"""
        nrow = 0; ncol = 2; padx = pady = 5
        self.sz_spn = wx.FlexGridSizer(nrow, ncol, padx, pady)
        vpad = 6
        self.sz_spn.Add(self.lump_lab, 0, wx.ALIGN_RIGHT|wx.RIGHT, 10)
        self.sz_spn.Add(self.lump_spn, 0, wx.ALIGN_RIGHT|wx.BOTTOM, vpad)
        self.sz_spn.Add(self.bins_lab, 0, wx.ALIGN_RIGHT|wx.RIGHT, 10)
        self.sz_spn.Add(self.bins_spn, 0, wx.ALIGN_RIGHT|wx.BOTTOM, vpad)
        self.sz_spn.Add(self.thresh_lab, 0, wx.ALIGN_RIGHT|wx.RIGHT, 10)
        self.sz_spn.Add(self.thresh_spn, 0, wx.ALIGN_RIGHT|wx.BOTTOM, vpad)

	self.sz_cho = wx.BoxSizer(wx.VERTICAL)
        self.sz_cho.Add(self.matl_cho,  0, wx.EXPAND|wx.ALIGN_CENTER|wx.TOP, 5)
        self.sz_cho.Add(self.read_cho,  0, wx.EXPAND|wx.ALIGN_CENTER|wx.TOP, 5)
        self.sz_cho.Add(self.hkls_but,  0, wx.ALIGN_CENTER|wx.TOP, 5)

        nrow = 0; ncol = 2; padx = pady = 5
	self.sz_orient = wx.FlexGridSizer(nrow, ncol, padx, pady)
	self.sz_orient.AddGrowableCol(1, 1)
        self.sz_orient.Add(self.angle_lab, 0, wx.EXPAND|wx.ALIGN_RIGHT|wx.RIGHT, 10)
        self.sz_orient.Add(self.angle_flt, 1, wx.EXPAND|wx.ALIGN_LEFT)
        self.sz_orient.Add(self.axis1_lab, 0, wx.EXPAND|wx.ALIGN_RIGHT|wx.RIGHT, 10)
        self.sz_orient.Add(self.axis1_flt, 1, wx.EXPAND|wx.ALIGN_LEFT)
        self.sz_orient.Add(self.axis2_lab, 0, wx.EXPAND|wx.ALIGN_RIGHT|wx.RIGHT, 10)
        self.sz_orient.Add(self.axis2_flt, 1, wx.EXPAND|wx.ALIGN_LEFT)
        self.sz_orient.Add(self.axis3_lab, 0, wx.EXPAND|wx.ALIGN_RIGHT|wx.RIGHT, 10)
        self.sz_orient.Add(self.axis3_flt, 1, wx.EXPAND|wx.ALIGN_LEFT)

        self.sz_opts = wx.BoxSizer(wx.HORIZONTAL)
        self.sz_opts.Add(self.sz_spn,    0,
                         wx.EXPAND|wx.ALIGN_CENTER|wx.RIGHT, 10)
        self.sz_opts.Add(self.sz_cho,    0,
                         wx.EXPAND|wx.ALIGN_CENTER|wx.RIGHT, 10)
        self.sz_opts.Add(self.sz_orient, 1, wx.EXPAND|wx.ALIGN_CENTER)

	self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.sizer.Add(self.tbarSizer, 0, wx.EXPAND|wx.ALIGN_CENTER)
	self.sizer.Add(self.sz_opts,   1, wx.EXPAND|wx.ALIGN_CENTER)

	return
    #
    # ============================== API
    #
    #                     ========== *** Access Methods
    #

    #
    #                     ========== *** Event Callbacks
    #
    def OnSelectHKLs(self, evt):
        """Callback for hkl_but"""

        app = wx.GetApp()
        exp = app.ws
        mat = exp.activeMaterial

        dlg = hklsDlg(self, wx.NewId(), mat)

        if dlg.ShowModal() == wx.ID_OK:
            mat.planeData.exclusions = dlg.getExclusions()
            pass

        return

    def OnChooseMatl(self, e):
        """Change material"""
        app = wx.GetApp()
        exp = app.ws
        exp.activeMaterial = self.matl_cho.GetStringSelection()

        return

    def OnChooseReader(self, e):
        """Change reader"""
        app = wx.GetApp()
        exp = app.ws
        exp.activeReader = self.matl_cho.GetStringSelection()

        return

    pass # end class
#
# -----------------------------------------------END CLASS:  sphericalOptsPanel

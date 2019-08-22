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
import copy

import wx

from hexrd.wx.guiconfig    import WindowParameters as WP
from hexrd.wx.guiutil      import makeTitleBar, EmptyWindow
from hexrd.wx.selecthkls   import selectHKLsDialog as HklsDlg

import hexrd.xrd.xrdbase as xrdbase

if xrdbase.haveMultiProc:
    ncpus_DFLT = xrdbase.multiprocessing.cpu_count()
else:
    ncpus_DFLT = 1

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
        exp = wx.GetApp().ws
        ind_opts = exp.index_opts

        self.sz_titlebar = makeTitleBar(self, 'Indexing')
        self.method_cho = wx.Choice(self, wx.NewId(),
                                    choices=ind_opts.INDEX_CHOICES)
        self.method_cho.SetSelection(ind_opts.IND_FIBER)
        self.run_but  = wx.Button(self, wx.NewId(), 'Run Indexer')

        self.fiber_pan = FiberSearchPanel(self, wx.NewId())
        self.gspot_pan = GrainSpotterPanel(self, wx.NewId())

        return

    def __makeBindings(self):
        """Bind interactors"""
        self.Bind(wx.EVT_CHOICE, self.ChooseMethod, self.method_cho)
        self.Bind(wx.EVT_BUTTON, self.OnRunIndex, self.run_but)

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
    def OnRunIndex(self, e):
        """Run Indexer"""
        # Set indexer options and run
        exp = wx.GetApp().ws
        iopts = exp.index_opts

        if iopts.index_method == iopts.IND_FIBER:
            self.fiber_pan.set_index_opts()
            exp.run_indexer()

        print 'index results: ', exp.fitRMats
        return

    def ChooseMethod(self, e):
        """Choose method button has been pressed"""
        exp = wx.GetApp().ws
        ind_opts = exp.index_opts

        ind_opts.index_method = self.method_cho.GetSelection()

        if self.method_cho.GetSelection() == ind_opts.IND_FIBER:
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
        exp = wx.GetApp().ws
        iopts = exp.index_opts

        self.tbarSizer = makeTitleBar(self, 'Fiber Search Options',
                                      color=WP.BG_COLOR_PANEL1_TITLEBAR)

        # checkboxes

        self.friedel_cbox = wx.CheckBox(self, wx.NewId(), 'Friedel Only')
        self.friedel_cbox.SetValue(iopts.friedelOnly)
        self.claims_cbox = wx.CheckBox(self, wx.NewId(), 'Preserve Claims')
        self.claims_cbox.SetValue(iopts.preserveClaims)
        self.refine_cbox  = wx.CheckBox(self, wx.NewId(), 'Do Refinement')
        self.refine_cbox.SetValue(iopts.doRefinement)
        self.multi_cbox  = wx.CheckBox(self, wx.NewId(), 'Use Multiprocessing')
        self.multi_cbox.SetValue(iopts.doMultiProc)

        # value boxes

        self.etol_lab = wx.StaticText(self, wx.NewId(), 'Eta Tolerance',
                                      style=wx.ALIGN_RIGHT)
        self.etol_txt = wx.TextCtrl(self, wx.NewId(), value=str(iopts.etaTol),
                                    style=wx.RAISED_BORDER)

        self.otol_lab = wx.StaticText(self, wx.NewId(), 'Omega Tolerance',
                                      style=wx.ALIGN_RIGHT)
        self.otol_txt = wx.TextCtrl(self, wx.NewId(), value=str(iopts.omeTol),
                                    style=wx.RAISED_BORDER)

        self.steps_lab = wx.StaticText(self, wx.NewId(), 'Number of Steps',
                                      style=wx.ALIGN_RIGHT)
        self.steps_spn = wx.SpinCtrl(self, wx.NewId(),
                                     min=36, max=36000, initial=iopts.nsteps)

        label = 'Minimum Completeness'
        self.comp_lab = wx.StaticText(self, wx.NewId(), label,
                                      style=wx.ALIGN_RIGHT)
        self.comp_txt = wx.TextCtrl(self, wx.NewId(),
                                    value=str(iopts.minCompleteness),
                                    style=wx.RAISED_BORDER)

        label = 'Minimum Fraction Claimed'
        self.claim_lab = wx.StaticText(self, wx.NewId(), label,
                                      style=wx.ALIGN_RIGHT)
        self.claim_txt = wx.TextCtrl(self, wx.NewId(),
                                     value=str(iopts.minPctClaimed),
                                     style=wx.RAISED_BORDER)

        label = 'Number of CPUs'
        self.ncpus_lab = wx.StaticText(self, wx.NewId(), label,
                                      style=wx.ALIGN_RIGHT)
        self.ncpus_spn = wx.SpinCtrl(self, wx.NewId(),
                                     min=1, max=ncpus_DFLT, initial=ncpus_DFLT)

        label = 'Quit After This Many'
        self.qafter_lab = wx.StaticText(self, wx.NewId(), label,
                                      style=wx.ALIGN_RIGHT)
        self.qafter_spn = wx.SpinCtrl(self, wx.NewId(),
                                     min=0, max=100000, initial=0)

        self.hkls_but = wx.Button(self, wx.NewId(), 'HKLs')

        return

    def __makeBindings(self):
        """Bind interactors"""
        self.Bind(wx.EVT_BUTTON, self.OnRunHKLs, self.hkls_but)
        return

    def __makeSizers(self):
        """Lay out the interactors"""
        # checkbox sizer

        self.cbsizer = wx.BoxSizer(wx.VERTICAL) # checkboxes
        self.cbsizer.Add(self.friedel_cbox, 0, wx.ALIGN_LEFT)
        self.cbsizer.Add(self.claims_cbox, 0, wx.ALIGN_LEFT)
        self.cbsizer.Add(self.refine_cbox, 0, wx.ALIGN_LEFT)
        self.cbsizer.Add(self.multi_cbox, 0, wx.ALIGN_LEFT)

        # value controls

        nrow = 0; ncol = 2; padx = pady = 5
        self.valsizer = wx.FlexGridSizer(nrow, ncol, padx, pady)
        self.valsizer.AddGrowableCol(1, 1)
        #self.valsizer.AddGrowableRow(num, proportion)
        self.valsizer.Add(self.etol_lab, 0, wx.EXPAND|wx.ALIGN_RIGHT)
        self.valsizer.Add(self.etol_txt, 0, wx.EXPAND|wx.ALIGN_LEFT)
        self.valsizer.Add(self.otol_lab, 0, wx.EXPAND|wx.ALIGN_RIGHT)
        self.valsizer.Add(self.otol_txt, 0, wx.EXPAND|wx.ALIGN_LEFT)
        self.valsizer.Add(self.steps_lab, 0, wx.EXPAND|wx.ALIGN_RIGHT)
        self.valsizer.Add(self.steps_spn, 0, wx.EXPAND|wx.ALIGN_LEFT)
        self.valsizer.Add(self.comp_lab, 0, wx.EXPAND|wx.ALIGN_RIGHT)
        self.valsizer.Add(self.comp_txt, 0, wx.EXPAND|wx.ALIGN_LEFT)
        self.valsizer.Add(self.claim_lab, 0, wx.EXPAND|wx.ALIGN_RIGHT)
        self.valsizer.Add(self.claim_txt, 0, wx.EXPAND|wx.ALIGN_LEFT)
        self.valsizer.Add(self.ncpus_lab, 0, wx.EXPAND|wx.ALIGN_RIGHT)
        self.valsizer.Add(self.ncpus_spn, 0, wx.EXPAND|wx.ALIGN_LEFT)
        self.valsizer.Add(self.qafter_lab, 0, wx.EXPAND|wx.ALIGN_RIGHT)
        self.valsizer.Add(self.qafter_spn, 0, wx.EXPAND|wx.ALIGN_LEFT)
        self.valsizer.Add(EmptyWindow(self), 0)
        self.valsizer.Add(self.hkls_but, 0, wx.EXPAND|wx.ALIGN_LEFT)


        # Main Sizer

        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.tbarSizer, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.sizer.Add(self.cbsizer, 0, wx.TOP, 10)
        self.sizer.Add(self.valsizer, 0, wx.TOP, 10)

        return
    #
    # ============================== API
    #
    #                     ========== *** Access Methods
    #
    def set_index_opts(self):
        """Set index options from interactors"""
        exp = wx.GetApp().ws
        iopts = exp.index_opts

        # iopts.fsHKLs = tuple(exp.activeMaterial.planeData.hkls)
        iopts.preserveClaims = self.claims_cbox.GetValue()
        iopts.friedelOnly = self.friedel_cbox.GetValue()
        iopts.doRefinement = self.refine_cbox.GetValue()
        iopts.doMultiProc = self.multi_cbox.GetValue()
        iopts.etaTol = float(self.etol_txt.GetValue())
        iopts.omeTol=float(self.otol_txt.GetValue())
        iopts.minCompleteness=float(self.comp_txt.GetValue())
        iopts.minPctClaimed=float(self.claim_txt.GetValue())
        iopts.nsteps = self.steps_spn.GetValue()

        ncpusVal = self.ncpus_spn.GetValue()
        iopts.nCPUs = ncpusVal if ncpusVal else None

        qafterVal = self.qafter_spn.GetValue()
        iopts.quitAfter = qafterVal if qafterVal else None

        iopts.dspTol=None # no interactor yet

        return
    #
    #                     ========== *** Event Callbacks
    #
    def OnRunHKLs(self, evt):
        """Select HKLs to use for indexing"""
        exp = wx.GetApp().ws
        iopts = exp.index_opts
        pd = exp.activeMaterial.planeData

        hkls_dlg = HklsDlg(self, wx.NewId(), exp.activeMaterial)

        if hkls_dlg.ShowModal() == wx.ID_OK:
            # pd.exclusions = hkls_dlg.getExclusions()
            iopts.fsHKLs = hkls_dlg.getTuples()
            print 'hkls: ', iopts.fsHKLs

        return

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

        self.pfit_cbox = wx.CheckBox(self, wx.NewId(), 'Position Fit')

        label = 'Minimum Completeness'
        self.comp_lab = wx.StaticText(self, wx.NewId(), label,
                                      style=wx.ALIGN_CENTER)
        self.comp_txt = wx.TextCtrl(self, wx.NewId(), value='0.5',
                                    style=wx.RAISED_BORDER)

        label = 'Minimum Fraction of G-Vectors'
        self.fracG_lab = wx.StaticText(self, wx.NewId(), label,
                                      style=wx.ALIGN_CENTER)
        self.fracG_txt = wx.TextCtrl(self, wx.NewId(), value='0.5',
                                    style=wx.RAISED_BORDER)

        label = 'Sigmas'
        self.sigmas_lab = wx.StaticText(self, wx.NewId(), label,
                                      style=wx.ALIGN_CENTER)
        self.sigmas_txt = wx.TextCtrl(self, wx.NewId(), value='2.0',
                                    style=wx.RAISED_BORDER)

        self.trials_lab = wx.StaticText(self, wx.NewId(), 'Number of Trials',
                                      style=wx.ALIGN_CENTER)
        self.trials_spn = wx.SpinCtrl(self, wx.NewId(),
                                     min=1000, max=1000000, initial=100000)


        return

    def __makeBindings(self):
        """Bind interactors"""
        return

    def __makeSizers(self):
        """Lay out the interactors"""
        # value controls

        nrow = 0; ncol = 2; padx = pady = 5
        self.valsizer = wx.FlexGridSizer(nrow, ncol, padx, pady)
        self.valsizer.AddGrowableCol(1, 1)
        #self.valsizer.AddGrowableRow(num, proportion)
        self.valsizer.Add(self.comp_lab, 0, wx.EXPAND|wx.ALIGN_RIGHT)
        self.valsizer.Add(self.comp_txt, 0, wx.EXPAND|wx.ALIGN_LEFT)
        self.valsizer.Add(self.fracG_lab, 0, wx.EXPAND|wx.ALIGN_RIGHT)
        self.valsizer.Add(self.fracG_txt, 0, wx.EXPAND|wx.ALIGN_LEFT)
        self.valsizer.Add(self.sigmas_lab, 0, wx.EXPAND|wx.ALIGN_RIGHT)
        self.valsizer.Add(self.sigmas_txt, 0, wx.EXPAND|wx.ALIGN_LEFT)
        self.valsizer.Add(self.trials_lab, 0, wx.EXPAND|wx.ALIGN_RIGHT)
        self.valsizer.Add(self.trials_spn, 0, wx.EXPAND|wx.ALIGN_LEFT)

        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.tbarSizer, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.sizer.Add(self.pfit_cbox, 0, wx.TOP, 10)
        self.sizer.Add(self.valsizer, 0, wx.TOP, 10)

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

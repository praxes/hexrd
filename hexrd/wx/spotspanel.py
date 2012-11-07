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
"""Panel for spot finding"""

import wx
import numpy

from hexrd.wx.guiconfig    import WindowParameters as WP
from hexrd.wx.guiutil      import makeTitleBar, EmptyWindow, ResetChoice
from hexrd.wx.logwindows   import logWindow
from hexrd.wx.selecthkls   import selectHKLsDialog as HklsDlg 

from hexrd.xrd.crystallography import processWavelength

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
        self.tbar_raw = makeTitleBar(self, 'Raw Spots',
                                     color=WP.BG_COLOR_TITLEBAR_PANEL1)
        self.tbar_ind = makeTitleBar(self, 'Spots for Indexing',
                                     color=WP.BG_COLOR_TITLEBAR_PANEL1)

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
            self, wx.NewId(), value='500', 
            style=wx.RAISED_BORDER)

        # Min PX

        self.minpx_lab = wx.StaticText(
            self, wx.NewId(),
            'Min PX', style=wx.ALIGN_CENTER)
        self.minpx_txt = wx.TextCtrl(
            self, wx.NewId(), value='4', 
            style=wx.RAISED_BORDER)

        # Spots info

        self.aread_lab = wx.StaticText(
            self, wx.NewId(), 
            'Active Reader', style=wx.ALIGN_RIGHT)
        self.aread_cho = wx.Choice(self, wx.NewId(), choices=['reader list'])
        
        self.rdr_lab = wx.StaticText(
            self, wx.NewId(), 
            'Used Readers', style=wx.ALIGN_RIGHT)
        self.rdr_lbx =  wx.ListBox(self, wx.NewId(), choices = ['r1', 'r2'])

        self.nspot_lab = wx.StaticText(
            self, wx.NewId(), 
            'Number of Spots', style=wx.ALIGN_RIGHT)
        self.nspot_txt = wx.TextCtrl(
            self, wx.NewId(), value='0', 
            style=wx.RAISED_BORDER)
        

        # Run button
        
        self.run  = wx.Button(self, wx.NewId(), 'Add Spots')
        self.clear_but  = wx.Button(self, wx.NewId(), 'Clear Spots')

        # Spots for Indexing info

        self.amat_lab = wx.StaticText(
            self, wx.NewId(), 
            'Active Material', style=wx.ALIGN_RIGHT)
        self.amat_cho = wx.Choice(self, wx.NewId(), choices=['mat list'])
	self.Bind(wx.EVT_CHOICE, self.OnMatChoice, self.aread_cho)

        
        self.hkls_lab = wx.StaticText(
            self, wx.NewId(), 
            '', style=wx.ALIGN_RIGHT)
        self.hkls_but = wx.Button(self, wx.NewId(), 'HKLs')

        self.nspotind_lab = wx.StaticText(
            self, wx.NewId(), 
            'Number of Spots', style=wx.ALIGN_RIGHT)
        self.nspotind_txt = wx.TextCtrl(
            self, wx.NewId(), value='0', 
            style=wx.RAISED_BORDER)
        
        # Run button for indexing spots
        
        self.run_ind  = wx.Button(self, wx.NewId(), 'Process Spots\nfor Indexing')

        return

    def __makeBindings(self):
        """Bind interactors"""
        self.Bind(wx.EVT_TEXT_ENTER, self.OnThreshold,  self.thresh_txt)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnMinPX,      self.minpx_txt)

        self.Bind(wx.EVT_BUTTON, self.OnRun, self.run)
        self.Bind(wx.EVT_BUTTON, self.OnRunInd, self.run_ind)
        self.Bind(wx.EVT_BUTTON, self.OnRunHKLs, self.hkls_but)
        self.Bind(wx.EVT_BUTTON, self.OnClearBut, self.clear_but)

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
        self.fgSizer.Add(self.thresh_lab,  0, wx.ALIGN_RIGHT)
        self.fgSizer.Add(self.thresh_txt,  0, wx.ALIGN_RIGHT)
        #  min PX
        self.fgSizer.Add(self.minpx_lab,   0, wx.ALIGN_RIGHT)
        self.fgSizer.Add(self.minpx_txt,   0, wx.ALIGN_RIGHT)
        self.fgSizer.Add(self.run,         0, wx.ALIGN_RIGHT)
        self.fgSizer.Add(self.clear_but,   0, wx.ALIGN_RIGHT)

        # ========== Raw Info Sizer

        nrow = 0; ncol = 2; padx = 5; pady = 5
        self.rawinfo_sizer = wx.FlexGridSizer(nrow, ncol, padx, pady) 
        self.rawinfo_sizer.AddGrowableCol(1, 1)
        self.rawinfo_sizer.Add(self.aread_lab, 0, wx.ALIGN_RIGHT)
        self.rawinfo_sizer.Add(self.aread_cho, 0, wx.ALIGN_RIGHT)
        self.rawinfo_sizer.Add(self.rdr_lab, 0, wx.ALIGN_RIGHT)
        self.rawinfo_sizer.Add(self.rdr_lbx, 0, wx.ALIGN_RIGHT)
        self.rawinfo_sizer.Add(self.nspot_lab, 0, wx.ALIGN_RIGHT)
        self.rawinfo_sizer.Add(self.nspot_txt, 0, wx.ALIGN_RIGHT)

        # ========== Raw Options Sizer

	self.rawopt_sizer = wx.BoxSizer(wx.VERTICAL)
        self.rawopt_sizer.Add(self.cbsizer,   0, wx.ALIGN_LEFT|wx.TOP, padtop)
        self.rawopt_sizer.Add(self.fgSizer,   0, wx.ALIGN_LEFT|wx.TOP, 2*padtop)

        
        # ========== Raw Sizer

	self.raw_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.raw_sizer.Add(self.rawinfo_sizer, 0, wx.ALIGN_LEFT)
        self.raw_sizer.Add(EmptyWindow(self), 1, wx.EXPAND)
        self.raw_sizer.Add(self.rawopt_sizer, 0)

        # ========== Indexing Info Sizer

        nrow = 0; ncol = 2; padx = 5; pady = 5
        self.indinfo_sizer = wx.FlexGridSizer(nrow, ncol, padx, pady) 
        self.indinfo_sizer.AddGrowableCol(1, 1)
        self.indinfo_sizer.Add(self.amat_lab, 0, wx.ALIGN_RIGHT)
        self.indinfo_sizer.Add(self.amat_cho, 0, wx.ALIGN_RIGHT)
        self.indinfo_sizer.Add(self.hkls_lab, 0, wx.ALIGN_RIGHT)
        self.indinfo_sizer.Add(self.hkls_but, 0, wx.ALIGN_RIGHT)
        self.indinfo_sizer.Add(self.nspotind_lab, 0, wx.ALIGN_RIGHT)
        self.indinfo_sizer.Add(self.nspotind_txt, 0, wx.ALIGN_RIGHT)
        self.indinfo_sizer.Add(EmptyWindow(self), 0, wx.TOP, padtop)
        self.indinfo_sizer.Add(self.run_ind, 0, wx.ALIGN_RIGHT)

        # ========== Main Sizer

	self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.sizer.Add(self.tbarSizer, 0, wx.EXPAND|wx.ALIGN_CENTER)
	self.sizer.Add(self.tbar_raw,  0, wx.EXPAND|wx.ALIGN_CENTER)
        self.sizer.Add(self.raw_sizer, 0, wx.EXPAND|wx.TOP, padtop)
	self.sizer.Add(self.tbar_ind,  0, wx.EXPAND|wx.ALIGN_CENTER)
        self.sizer.Add(self.indinfo_sizer, 0, wx.TOP, padtop)

	return

    def _run_spots(self, log=None):
        """run spot finder for use in log window"""
        exp = wx.GetApp().ws
        log.write('running spot finder ...')
        exp.find_raw_spots()
        log.write('DONE')
        
        return
    #
    # ============================== API
    #
    #                     ========== *** Access Methods
    #
    def updateFromExp(self):
        """Update all subwindows"""
        exp  = wx.GetApp().ws
        
        ResetChoice(self.aread_cho, exp.readerNames, exp.activeReader.name)
        ResetChoice(self.amat_cho, exp.matNames, exp.activeMaterial.name)
        self.nspot_txt.ChangeValue(str(len(exp.raw_spots)))

        spots_ind = exp.spots_for_indexing
        if hasattr(spots_ind, 'nTThAssoc'):
            nspot_assoc = len(spots_ind) - numpy.sum(spots_ind.nTThAssoc == 0)
        else:
            nspot_assoc = 0

        self.nspotind_txt.ChangeValue(str(nspot_assoc))
        
        self.rdr_lbx.Set(exp.spot_readers)
        
        # Reset sizers

        self.Layout()
        
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

    def OnMatChoice(self, e):
        """Select new material"""
        exp = wx.GetApp().ws

        sel = self.amat_cho.GetSelection()
        if sel >= 0:
            exp.activeMaterial = sel
            self.updateFromExp()
            pass
    
        return

    def OnClearBut(self, evt):
        """Clear spots array"""
        exp = wx.GetApp().ws
        exp.clear_spots()
        self.updateFromExp()
        
        return
    
    def OnRun(self, evt):
        """Callback for run"""
        #
        # Fill in spot options from the form
        #
        exp = wx.GetApp().ws
        opts = exp.spotOpts
        #
        opts.thresh = int(self.thresh_txt.GetValue())
        opts.minpx = int(self.minpx_txt.GetValue())
        #
        action = {
            'exec'  : self._run_spots,
            'args'  : (),
            'kwargs': dict()
            }
        logwin = logWindow(self, wx.NewId(), action, 'Finding Spots')
        logwin.ShowModal()

        self.updateFromExp()
        return

    def OnRunInd(self, evt):
        """Run processing for indexing"""
        exp = wx.GetApp().ws
        exp.get_spots_ind()
        
        self.updateFromExp()
        return

    def OnRunHKLs(self, evt):
        """Select HKLs to use for indexing"""
        exp = wx.GetApp().ws
        hkls_dlg = HklsDlg(self, wx.NewId(), exp.activeMaterial)

        if hkls_dlg.ShowModal() == wx.ID_OK:
            exp.activeMaterial.planeData.exclusions = hkls_dlg.getExclusions()
            pass
        
        return

    pass # end class
#
# -----------------------------------------------END CLASS:  spotsPanel

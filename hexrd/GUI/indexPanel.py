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
    #
    # Class data
    #
    INDEX_CHOICES = ['Fiber Search', 'GrainSpotter']
    INDEX_CHOICE_IDS = [IND_FIBER, IND_GSPOT] = range(2)
    
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

        self.sz_titlebar = makeTitleBar(self, 'Indexing')
        self.method_cho = wx.Choice(self, wx.NewId(),
                                    choices=self.INDEX_CHOICES)
        self.method_cho.SetSelection(self.IND_FIBER)
        self.run_but  = wx.Button(self, wx.NewId(), 'Run Indexer')

        self.fiber_pan = FiberSearchPanel(self, wx.NewId())
        self.gspot_pan = GrainSpotterPanel(self, wx.NewId())

        return

    def __makeBindings(self):
        """Bind interactors"""
	self.Bind(wx.EVT_CHOICE, self.ChooseMethod, self.method_cho)
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
    def ChooseMethod(self, e):
        """Choose method button has been pressed"""
        if self.method_cho.GetSelection() == self.IND_FIBER:
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

        self.tbarSizer = makeTitleBar(self, 'Fiber Search Options', 
                                      color=WP.BG_COLOR_PANEL1_TITLEBAR)

        # checkboxes
        
        self.friedel_cbox = wx.CheckBox(self, wx.NewId(), 'Friedel Only')
        self.claims_cbox = wx.CheckBox(self, wx.NewId(), 'Preserve Claiims')
        self.refine_cbox  = wx.CheckBox(self, wx.NewId(), 'Do Refinement')
        self.multi_cbox  = wx.CheckBox(self, wx.NewId(), 'Use Multiprocessing')

        # value boxes
        
        self.etol_lab = wx.StaticText(self, wx.NewId(), 'Eta Tolerance',
                                      style=wx.ALIGN_CENTER)
        self.etol_txt = wx.TextCtrl(self, wx.NewId(), value='',
                                    style=wx.RAISED_BORDER)
        
        self.otol_lab = wx.StaticText(self, wx.NewId(), 'Omega Tolerance',
                                      style=wx.ALIGN_CENTER)
        self.otol_txt = wx.TextCtrl(self, wx.NewId(), value='',
                                    style=wx.RAISED_BORDER)
        
        self.steps_lab = wx.StaticText(self, wx.NewId(), 'Number of Steps',
                                      style=wx.ALIGN_CENTER)
        self.steps_spn = wx.SpinCtrl(self, wx.NewId(),
                                     min=10, max=360000, initial=360)

        label = 'Minimum Completeness'
        self.comp_lab = wx.StaticText(self, wx.NewId(), label,
                                      style=wx.ALIGN_CENTER)
        self.comp_txt = wx.TextCtrl(self, wx.NewId(), value='0.5',
                                    style=wx.RAISED_BORDER)

        label = 'Minimum Percent Claimed'
        self.claim_lab = wx.StaticText(self, wx.NewId(), label,
                                      style=wx.ALIGN_CENTER)
        self.claim_txt = wx.TextCtrl(self, wx.NewId(), value='50',
                                    style=wx.RAISED_BORDER)

        return

    def __makeBindings(self):
        """Bind interactors"""
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

    #
    #                     ========== *** Event Callbacks
    #
    
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

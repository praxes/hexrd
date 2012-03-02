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
#  Copyright-Info-Goes-Here
#
"""Panel for materials input
"""
#
import wx

from hexrd.XRD.Material import Material

from hexrd.GUI.guiConfig    import WindowParameters as WP
from hexrd.GUI.guiUtilities import ResetChoice, AddSpacer
#
# ---------------------------------------------------CLASS:  matPanel
#
class matPanel(wx.Panel):
    """matPanel """
    def __init__(self, parent, id, **kwargs):
	"""Constructor for matPanel."""
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
        exp = wx.GetApp().ws
        mat = exp.activeMaterial
        
        self.__makeTitleBar('Materials')
        #
        #  ========== Header
        #
        #  Material List
        #
        self.curr_lab = wx.StaticText(self, wx.NewId(), 
                                        'Active Material', style=wx.ALIGN_CENTER)
        self.mats_cho = wx.Choice(self, wx.NewId(), 
                                  choices=[m.name for m in exp.matList])
        self.new_but  = wx.Button(self, wx.NewId(), 'New Material')
        #
        #  Material Name
        #
        self.name_lab = wx.StaticText(self, wx.NewId(), 
                                        'MATERIAL NAME', style=wx.ALIGN_CENTER)
        self.name_txt = wx.TextCtrl(self, wx.NewId(), value=Material.DFLT_NAME, 
                                      style=wx.RAISED_BORDER|wx.TE_PROCESS_ENTER)
        #
        #  Categories
        #
        #  ========== Lattice Params
        #
        self.lp_a_lab = wx.StaticText(self, wx.NewId(), 'a', style=wx.ALIGN_CENTER)
        self.lp_a_txt = wx.TextCtrl(self, wx.NewId(), value='0', 
                                    style=wx.RAISED_BORDER|wx.TE_PROCESS_ENTER)

        self.lp_b_lab = wx.StaticText(self, wx.NewId(), 'b', style=wx.ALIGN_CENTER)
        self.lp_b_txt = wx.TextCtrl(self, wx.NewId(), value='0', 
                                    style=wx.RAISED_BORDER|wx.TE_PROCESS_ENTER)

        self.lp_c_lab = wx.StaticText(self, wx.NewId(), 'c', style=wx.ALIGN_CENTER)
        self.lp_c_txt = wx.TextCtrl(self, wx.NewId(), value='0', 
                                    style=wx.RAISED_BORDER|wx.TE_PROCESS_ENTER)

        self.alpha_lab = wx.StaticText(self, wx.NewId(), 'alpha', style=wx.ALIGN_CENTER)
        self.alpha_txt = wx.TextCtrl(self, wx.NewId(), value='90', 
                                     style=wx.RAISED_BORDER|wx.TE_PROCESS_ENTER)

        self.beta_lab = wx.StaticText(self, wx.NewId(), 'beta', style=wx.ALIGN_CENTER)
        self.beta_txt = wx.TextCtrl(self, wx.NewId(), value='90', 
                                    style=wx.RAISED_BORDER|wx.TE_PROCESS_ENTER)
        
        self.gamma_lab = wx.StaticText(self, wx.NewId(), 'gamma', style=wx.ALIGN_CENTER)
        self.gamma_txt = wx.TextCtrl(self, wx.NewId(), value='90', 
                                     style=wx.RAISED_BORDER|wx.TE_PROCESS_ENTER)
        
        self.units_lab  = wx.StaticText(self, wx.NewId(), 'UNITS', style=wx.ALIGN_CENTER)
        self.dunits_cho = wx.Choice(self, wx.NewId(), choices=['angstroms'])
        self.dunits_cho.SetSelection(0)
        self.aunits_cho = wx.Choice(self, wx.NewId(), choices=['degrees'])
        self.aunits_cho.SetSelection(0)
        #
        #  Save list of lattice parameter windows.
        #
        self.lpWins = [self.lp_a_txt,  self.lp_b_txt, self.lp_c_txt,
                       self.alpha_txt, self.beta_txt, self.gamma_txt]
        self.wDict = {
            0: self.lp_a_txt,
            1: self.lp_b_txt,
            2: self.lp_c_txt,
            3: self.alpha_txt,
            4: self.beta_txt,
            5: self.gamma_txt
            }

        #
        #  ========== Space group info
        #
        self.sg_lab = wx.StaticText(self, wx.NewId(), 'Space Group', 
                                    style=wx.ALIGN_CENTER)
        self.sg_spn = wx.SpinCtrl(self, wx.NewId(), min=1, max=230, initial=mat.spaceGroup.sgnum)

        self.hall_lab = wx.StaticText(self, wx.NewId(), 'Hall Symbol', 
                                      style=wx.ALIGN_CENTER)
        self.hall_txt = wx.TextCtrl(self, wx.NewId(), value='', 
                                    style=wx.RAISED_BORDER|wx.TE_READONLY)

        self.herm_lab = wx.StaticText(self, wx.NewId(), 'Hermann-Mauguin', 
                                      style=wx.ALIGN_CENTER)
        self.herm_txt = wx.TextCtrl(self, wx.NewId(), value='', 
                                    style=wx.RAISED_BORDER|wx.TE_READONLY)

        self.laue_lab = wx.StaticText(self, wx.NewId(), 'Laue Group', 
                                      style=wx.ALIGN_CENTER)
        self.laue_txt = wx.TextCtrl(self, wx.NewId(), value='', 
                                    style=wx.RAISED_BORDER|wx.TE_READONLY)

        self.ltype_lab = wx.StaticText(self, wx.NewId(), 'Lattice Type', 
                                       style=wx.ALIGN_CENTER)
        self.ltype_txt = wx.TextCtrl(self, wx.NewId(), value='', 
                                     style=wx.RAISED_BORDER|wx.TE_READONLY)

        self.hkls_lab = wx.StaticText(self, wx.NewId(), 'HKLs Max (sum of squares)', 
                                      style=wx.ALIGN_CENTER)
        self.hkls_txt = wx.TextCtrl(self, wx.NewId(), value='10', 
                                    style=wx.RAISED_BORDER|wx.TE_PROCESS_ENTER)
        
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
        self.Bind(wx.EVT_CHOICE,     self.OnMatChoice,   self.mats_cho)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnNameChange,  self.name_txt)
        self.Bind(wx.EVT_BUTTON,     self.OnNewMaterial, self.new_but)
        self.Bind(wx.EVT_SPINCTRL,   self.OnSpaceGroup,  self.sg_spn)

        # HKL sum of squares
        self.Bind(wx.EVT_TEXT_ENTER, self.OnHklMax, self.hkls_txt)
        
        # Lattice Parameters
        for w in self.lpWins:
            self.Bind(wx.EVT_TEXT_ENTER, self.OnLPChange,  w)
            pass

        return

    def __makeSizers(self):
	"""Lay out the interactors"""
	#  Header
        ncol = 3; nrow = 3; padx = pady = 5;
        headsizer = wx.FlexGridSizer(nrow, ncol, padx, pady) # m x n, paddings
	#self.fgsizer.AddGrowableRow(num, proportion)
	headsizer.AddGrowableCol(1, 1)
	headsizer.AddGrowableCol(2, 1)
        # material list
        headsizer.Add(self.curr_lab, 0, wx.EXPAND|wx.ALIGN_CENTER)
        headsizer.Add(self.mats_cho, 0, wx.EXPAND|wx.ALIGN_CENTER)
        headsizer.Add(self.new_but,  0, wx.ALIGN_RIGHT)
        # material name
        headsizer.Add(self.name_lab, 0, wx.EXPAND|wx.ALIGN_CENTER)
        headsizer.Add(wx.Window(self, -1), 0, wx.EXPAND|wx.ALIGN_CENTER)
        headsizer.Add(self.name_txt, 0, wx.EXPAND|wx.ALIGN_CENTER)
        #
        # lattice params
        #
        nrow = 4; ncol = 4; padx = pady = 5;
        lpsizer = wx.FlexGridSizer(nrow, ncol, padx, pady) # m x n, paddings
	lpsizer.AddGrowableCol(1, 1)
	lpsizer.AddGrowableCol(3, 1)
        # units
        lpsizer.Add(self.units_lab,  0,wx.ALIGN_CENTER)
        lpsizer.Add(self.dunits_cho,  1, wx.ALIGN_CENTER)
        lpsizer.AddStretchSpacer()
        lpsizer.Add(self.aunits_cho,  1, wx.ALIGN_CENTER)
        # a / alpha
        lpsizer.Add(self.lp_a_lab,  0,wx.ALIGN_CENTER)
        lpsizer.Add(self.lp_a_txt,  1,wx.EXPAND|wx.ALIGN_CENTER)
        lpsizer.Add(self.alpha_lab, 0,wx.ALIGN_CENTER)
        lpsizer.Add(self.alpha_txt, 1,wx.EXPAND|wx.ALIGN_CENTER)
        # b / beta
        lpsizer.Add(self.lp_b_lab, 0,wx.ALIGN_CENTER)
        lpsizer.Add(self.lp_b_txt, 1,wx.EXPAND|wx.ALIGN_CENTER)
        lpsizer.Add(self.beta_lab, 0,wx.ALIGN_CENTER)
        lpsizer.Add(self.beta_txt, 1,wx.EXPAND|wx.ALIGN_CENTER)
        # c / gamma
        lpsizer.Add(self.lp_c_lab, 0,wx.ALIGN_CENTER)
        lpsizer.Add(self.lp_c_txt, 1,wx.EXPAND|wx.ALIGN_CENTER)
        lpsizer.Add(self.gamma_lab, 0,wx.ALIGN_CENTER)
        lpsizer.Add(self.gamma_txt, 1,wx.EXPAND|wx.ALIGN_CENTER)
        #
        #  Space group
        #
        nrow = 5; ncol = 3; padx = pady = 5;
        sgsizer = wx.FlexGridSizer(nrow, ncol, padx, pady) # m x n, paddings
	sgsizer.AddGrowableCol(1, 1)
	sgsizer.AddGrowableCol(2, 1)
        #  labels
        sgsizer.Add(self.sg_lab,   0, wx.ALIGN_CENTER)
        sgsizer.Add(self.hall_lab, 0, wx.ALIGN_CENTER)
        sgsizer.Add(self.herm_lab, 0, wx.ALIGN_CENTER)
        #  info
        sgsizer.Add(self.sg_spn,   0, wx.ALIGN_RIGHT)
        sgsizer.Add(self.hall_txt, 0, wx.ALIGN_CENTER)
        sgsizer.Add(self.herm_txt, 0, wx.ALIGN_CENTER)
        #  labels
        sgsizer.Add(wx.Window(self, -1), 0, wx.ALIGN_RIGHT)
        sgsizer.Add(self.laue_lab,  0, wx.ALIGN_CENTER)
        sgsizer.Add(self.ltype_lab, 0, wx.ALIGN_CENTER)
        #  info
        sgsizer.Add(wx.Window(self, -1), 0, wx.ALIGN_RIGHT)
        sgsizer.Add(self.laue_txt,  0, wx.ALIGN_CENTER)
        sgsizer.Add(self.ltype_txt, 0, wx.ALIGN_CENTER)
        #  HKLs lab
        sgsizer.Add(self.hkls_lab,  0, wx.ALIGN_CENTER)
        sgsizer.Add(wx.Window(self, -1), 0, wx.ALIGN_RIGHT)
        sgsizer.Add(wx.Window(self, -1), 0, wx.ALIGN_RIGHT)
        #  HKLs
        sgsizer.Add(self.hkls_txt,  0, wx.ALIGN_CENTER)
        sgsizer.Add(wx.Window(self, -1), 0, wx.ALIGN_RIGHT)
        sgsizer.Add(wx.Window(self, -1), 0, wx.ALIGN_RIGHT)
        #
        #  ==================== MAIN SIZER
        #
	self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.sizer.Add(self.titlebar, 0, wx.EXPAND|wx.ALIGN_CENTER)
	self.sizer.Add(headsizer,     1, wx.EXPAND|wx.ALIGN_CENTER)
        #
        #  ==================== MAIN SIZER
        #
	self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.sizer.Add(self.titlebar, 0,
                       wx.EXPAND|wx.ALIGN_CENTER|wx.BOTTOM, 10)
	self.sizer.Add(headsizer,     1, wx.EXPAND|wx.ALIGN_CENTER)

        AddSpacer(self, self.sizer, WP.BG_COLOR_TITLEBAR_PANEL1)
        
        self.sizer.Add(sgsizer, 1, wx.EXPAND|wx.ALIGN_CENTER)
        
        AddSpacer(self, self.sizer, WP.BG_COLOR_TITLEBAR_PANEL1)

	self.sizer.Add(lpsizer, 1, wx.EXPAND|wx.ALIGN_CENTER)

        #self.sizer.Show(self.titlebar, False)

	return
    #
    # ============================== API
    #
    #
    #                     ========== *** Access Methods
    #
    def updateFromExp(self):
        """Update interactors"""
        app = wx.GetApp()
        exp = app.ws
        mat = exp.activeMaterial

        # Materials line:  reset choice interactor
        ResetChoice(self.mats_cho, exp.matNames, mat.name)

        # Name line
        self.name_txt.ChangeValue(mat.name)

        # Space group info
        sg = mat.spaceGroup
        self.sg_spn.SetValue(mat.sgnum)
        self.hall_txt.ChangeValue(sg.HallSymbol)
        self.herm_txt.ChangeValue(sg.hermannMauguin)
        self.laue_txt.ChangeValue(sg.laueGroup)
        self.ltype_txt.ChangeValue(sg.latticeType)

        self.hkls_txt.ChangeValue(str(mat.hklMax))

        # Lattice parameters
        A  = 'angstrom'
        D  = 'degrees'
        BG = {True: (200,255, 200), False: (200, 200, 200)}

        reqP = mat.spaceGroup.reqParams
        lprm = mat.latticeParameters

        enable = lambda i: (self.wDict[i].Enable(i in reqP), 
                            self.wDict[i].SetBackgroundColour(BG[i in reqP]) )

        self.lp_a_txt.ChangeValue(str(lprm[0].getVal(A)))
        self.lp_b_txt.ChangeValue(str(lprm[1].getVal(A)))
        self.lp_c_txt.ChangeValue(str(lprm[2].getVal(A)))

        self.alpha_txt.ChangeValue(str(lprm[3].getVal(D)))
        self.beta_txt.ChangeValue (str(lprm[4].getVal(D)))
        self.gamma_txt.ChangeValue(str(lprm[5].getVal(D)))
        
        for i in range(6): enable(i)

        app.getCanvas().update()


        return
    #
    #                     ========== *** Event Callbacks
    #
    def OnHklMax(self, e):
        """Max sum of squares for HKLs has changed"""
        mat = wx.GetApp().ws.activeMaterial

        mat.hklMax = int(self.hkls_txt.GetValue())
        print mat.hklMax
        self.updateFromExp()

        return

    def OnLPChange(self, e):
        """Lattice parameter has changed"""
        mat = wx.GetApp().ws.activeMaterial

        lp = []
        for w in self.lpWins:
            if w.IsEnabled():
                lp.append(float(w.GetValue()))
                pass
            pass

        print 'new lattice params:  ', lp
        mat.latticeParameters = lp
        print 'new plane data params: ', mat.latticeParameters

        self.updateFromExp()

        return

    def OnSpaceGroup(self, e):
        """Space group number has been updated"""
        mat = wx.GetApp().ws.activeMaterial
        mat.sgnum = self.sg_spn.GetValue()
        self.updateFromExp()

        return

    def OnNewMaterial(self, e):
        """New material button has been pressed"""
        exp = wx.GetApp().ws
        exp.newMaterial()
        self.updateFromExp()

        return

    def OnNameChange(self, e):
        """Name entry has changed"""
        exp = wx.GetApp().ws
        n = str(self.name_txt.GetValue())
        names = exp.matNames
        if names.count(n) > 0:
            wx.MessageBox('name already in use')
        else:
            exp.activeMaterial.name = n
            pass

        self.updateFromExp()

        return

    def OnMatChoice(self, e):
        """Select new material"""
        exp = wx.GetApp().ws

        sel = self.mats_cho.GetSelection()
        if sel >= 0:
            exp.activeMaterial = sel
            self.updateFromExp()
            pass

        return

    pass # end class
#
# -----------------------------------------------END CLASS:  matPanel

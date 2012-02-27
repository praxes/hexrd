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
#  $Id: planeDataEditor.py 819 2011-04-25 16:30:55Z boyce6 $
#
"""Frame for editing plane data/calibrants

OBSOLETE:  due to be removed, but currently saving for reference 
           until Materials panel is complete
"""
import wx
import numpy

from guiConfig   import WindowParameters as WP
#
#  Module data
#
UNIT_DEG = 'degrees'
UNIT_RAD = 'radians'
AngleUnits = [UNIT_DEG, UNIT_RAD]
#
# ---------------------------------------------------CLASS:  PlaneDataPanel
#
class PlaneDataPanel(wx.Panel):
    #
    def __init__(self, parent, id, mat, title='Plane Data Panel'):
        """Instantiate a PlaneDataEditor

        parent -- parent window
        id     -- wx ID
        mat    -- material to be edited
"""        
        #
        wx.Panel.__init__(self, parent, id)
	#
        #  Data
        #
        self.mat   = mat
        self.pData = mat.planeData

        nHKL           = self.pData.getNhklRef()
        self.exclude   = numpy.zeros(nHKL, dtype=bool)
        self.exclude[5:] = True
        self.pData.set_exclusions(self.exclude)
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
        self.Show(True)

	return

    pass # class
    #
    # ============================== Internal Methods
    #
    def __makeObjects(self):
        """Add interactors"""

        self.__makeTitleBar('Plane Data Panel')
	#
        #  Text control for name
        #
        self.name_lab = wx.StaticText(self, wx.NewId(), 'Name', style=wx.ALIGN_CENTER)
        self.name_txt = wx.TextCtrl(self, wx.NewId(), value='', 
                                    style=wx.RAISED_BORDER|wx.TE_PROCESS_ENTER)
        self.name_txt.ChangeValue(self.mat.name)
	#
        #  Two-Theta and wavelength selectors, with units.
        #
        self.tthmin_lab = wx.StaticText(self, wx.NewId(), 'Two Theta Min', 
                                            style=wx.ALIGN_CENTER)
        self.tthmin_txt = wx.TextCtrl(self, wx.NewId(), value='0.0', 
                                          style=wx.RAISED_BORDER|wx.TE_PROCESS_ENTER)
        self.tthmin_uni = wx.Choice(self, wx.NewId(), choices=AngleUnits)
	#self.Bind(wx.EVT_CHOICE, self.OnChoice, self.choice)


        self.tthmax_lab = wx.StaticText(self, wx.NewId(), 'Two Theta Max', 
                                        style=wx.ALIGN_CENTER)
        self.tthmax_txt = wx.TextCtrl(self, wx.NewId(), value='20.0', 
                                      style=wx.RAISED_BORDER|wx.TE_PROCESS_ENTER)
        self.tthmax_uni = wx.Choice(self, wx.NewId(), choices=AngleUnits)

        self.wave_lab = wx.StaticText(self, wx.NewId(), 'Wavelength', 
                                      style=wx.ALIGN_CENTER)
        self.wave_txt = wx.TextCtrl(self, wx.NewId(), value='', 
                                    style=wx.RAISED_BORDER|wx.TE_PROCESS_ENTER)
        self.wave_txt.ChangeValue(str(self.pData.wavelength))
        self.wave_uni = wx.Choice(self, wx.NewId(), choices=AngleUnits)
        #
        #  Group selectors
        #
        self.laue_lab = wx.StaticText(self, wx.NewId(), 
                                        'Select the Laue group', 
                                        style=wx.ALIGN_RIGHT)
        self.laue_cho = wx.Choice(self, wx.NewId(), choices=['Laue Groups'])

        self.space_lab = wx.StaticText(self, wx.NewId(), 
                                        'Select the Space group', 
                                        style=wx.ALIGN_RIGHT)
        self.space_cho = wx.Choice(self, wx.NewId(), choices=['Space Groups'])
        #
        #  Add HKL list
        #
        self.hkls_clb =  wx.CheckListBox(self, wx.NewId(), choices = self.__getHKLs())
        [self.hkls_clb.Check(i, not self.exclude[i]) for i in range(len(self.exclude))]
        #
        #  Lattice Parameters
        #
        self.lp_a_lab = wx.StaticText(self, wx.NewId(), 'a', style=wx.ALIGN_CENTER)
        self.lp_a_txt = wx.TextCtrl(self, wx.NewId(), value='', 
                                        style=wx.RAISED_BORDER|wx.TE_PROCESS_ENTER)
        self.lp_b_lab = wx.StaticText(self, wx.NewId(), 'b', style=wx.ALIGN_CENTER)
        self.lp_b_txt = wx.TextCtrl(self, wx.NewId(), value='', 
                                        style=wx.RAISED_BORDER|wx.TE_PROCESS_ENTER)
        self.lp_c_lab = wx.StaticText(self, wx.NewId(), 'c', style=wx.ALIGN_CENTER)
        self.lp_c_txt = wx.TextCtrl(self, wx.NewId(), value='', 
                                        style=wx.RAISED_BORDER|wx.TE_PROCESS_ENTER)
        self.lp_alpha_lab = wx.StaticText(self, wx.NewId(), 'alpha', style=wx.ALIGN_CENTER)
        self.lp_alpha_txt = wx.TextCtrl(self, wx.NewId(), value='', 
                                        style=wx.RAISED_BORDER|wx.TE_PROCESS_ENTER)
        self.lp_beta_lab = wx.StaticText(self, wx.NewId(), 'beta', style=wx.ALIGN_CENTER)
        self.lp_beta_txt = wx.TextCtrl(self, wx.NewId(), value='', 
                                        style=wx.RAISED_BORDER|wx.TE_PROCESS_ENTER)
        self.lp_gamma_lab = wx.StaticText(self, wx.NewId(), 'gamma', style=wx.ALIGN_CENTER)
        self.lp_gamma_txt = wx.TextCtrl(self, wx.NewId(), value='', 
                                        style=wx.RAISED_BORDER|wx.TE_PROCESS_ENTER)

        return

    def __makeTitleBar(self, t):
        """Add titlebar"""
	self.titlebar = wx.StaticText(self, -1, t, 
					 style=wx.ALIGN_CENTER|wx.SIMPLE_BORDER) 
	self.titlebar.SetBackgroundColour(WP.TITLEBAR_BG_COLOR_FRAME)
        myToolTip = r"""
FRAME FOR editing the list of calibrants
"""
        self.titlebar.SetToolTipString(myToolTip)

        return

    def __makeBindings(self):
        """Bind interactors"""
        self.Bind(wx.EVT_TEXT_ENTER, self.OnRename, self.name_txt)

        return

    def __makeSizers(self):
	"""Lay out the interactors"""
        # for text label and control
        nrow = 2; ncol = 3; padx = pady = 5;
	self.paramSizer = wx.FlexGridSizer(nrow, ncol, padx, pady) 
	self.paramSizer.AddGrowableCol(1, 1)
	self.paramSizer.AddGrowableCol(2, 1)
        self.paramSizer.Add(self.name_lab, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.paramSizer.Add(self.name_txt, 1, wx.EXPAND|wx.ALIGN_CENTER)
        self.paramSizer.Add(wx.Window(self, -1), 1, wx.EXPAND|wx.ALIGN_CENTER)

        self.paramSizer.Add(self.wave_lab, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.paramSizer.Add(self.wave_txt, 1, wx.EXPAND|wx.ALIGN_CENTER)
        self.paramSizer.Add(self.wave_uni, 1, wx.EXPAND|wx.ALIGN_CENTER)

        self.paramSizer.Add(self.tthmin_lab, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.paramSizer.Add(self.tthmin_txt, 1, wx.EXPAND|wx.ALIGN_CENTER)
        self.paramSizer.Add(self.tthmin_uni, 1, wx.EXPAND|wx.ALIGN_CENTER)

        self.paramSizer.Add(self.tthmax_lab, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.paramSizer.Add(self.tthmax_txt, 1, wx.EXPAND|wx.ALIGN_CENTER)
        self.paramSizer.Add(self.tthmax_uni, 1, wx.EXPAND|wx.ALIGN_CENTER)
        #
        #  Group sizer
        #
        nrow = 2; ncol = 2; padx = pady = 5;
	self.gSizer = wx.FlexGridSizer(nrow, ncol, padx, pady) 
	self.gSizer.AddGrowableCol(1, 1)
        self.gSizer.Add(self.laue_lab, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.gSizer.Add(self.laue_cho, 1, wx.EXPAND|wx.ALIGN_CENTER)
        self.gSizer.Add(self.space_lab, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.gSizer.Add(self.space_cho, 1, wx.EXPAND|wx.ALIGN_CENTER)
        #
        #  Lattice Parameter sizer
        #
        nrow = 2; ncol = 2; padx = pady = 5;
	self.lpSizer = wx.FlexGridSizer(nrow, ncol, padx, pady) 
	self.lpSizer.AddGrowableCol(1, 1)        
        self.lpSizer.Add(self.lp_a_lab, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.lpSizer.Add(self.lp_a_txt, 1, wx.EXPAND|wx.ALIGN_CENTER)
        self.lpSizer.Add(self.lp_b_lab, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.lpSizer.Add(self.lp_b_txt, 1, wx.EXPAND|wx.ALIGN_CENTER)
        self.lpSizer.Add(self.lp_c_lab, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.lpSizer.Add(self.lp_c_txt, 1, wx.EXPAND|wx.ALIGN_CENTER)
        self.lpSizer.Add(self.lp_alpha_lab, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.lpSizer.Add(self.lp_alpha_txt, 1, wx.EXPAND|wx.ALIGN_CENTER)
        self.lpSizer.Add(self.lp_beta_lab, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.lpSizer.Add(self.lp_beta_txt, 1, wx.EXPAND|wx.ALIGN_CENTER)
        self.lpSizer.Add(self.lp_gamma_lab, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.lpSizer.Add(self.lp_gamma_txt, 1, wx.EXPAND|wx.ALIGN_CENTER)
        #
        #  Main fg sizer
        #
        nrow = 2; ncol = 2; padx = pady = 5;
	self.fgSizer = wx.FlexGridSizer(nrow, ncol, padx, pady) 
	self.fgSizer.AddGrowableRow(1, 1)
	self.fgSizer.Add(self.paramSizer,   1, wx.EXPAND|wx.ALIGN_CENTER)
	self.fgSizer.Add(self.gSizer,    1, wx.EXPAND|wx.ALIGN_CENTER)
	self.fgSizer.Add(self.hkls_clb,  1, wx.EXPAND | wx.ALIGN_CENTER | wx.TOP, 5)
	self.fgSizer.Add(self.lpSizer,  1, wx.EXPAND | wx.ALIGN_CENTER | wx.TOP, 5)
        #
        #  Main Sizer
        #
	self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.sizer.Add(self.titlebar,  0, wx.EXPAND|wx.ALIGN_CENTER)
	self.sizer.Add(self.fgSizer,   1, wx.EXPAND|wx.ALIGN_CENTER)

	return

    def __getHKLs(self):
        """Return HKLs as a string list"""
        hkls = self.pData.getHKLs(asStr=True, allHKLs=True)

        return hkls
    #
    # ============================== API
    #
    #
    #                     ========== *** Event Callbacks
    #
    def OnRename(self, e):
        """Name control has been modified"""
        self.mat.setName(self.name_txt.GetValue())

        return

#
# -----------------------------------------------END CLASS:  Plane Data Panel
#
#
# ---------------------------------------------------CLASS:  PlaneDataDialog
#
class PlaneDataDialog(wx.Dialog):
    """PlaneDataDialog """
    def __init__(self, parent, id, mat, **kwargs):
	"""Constructor for PlaneDataDialog"""
	#
	wx.Dialog.__init__(self, parent, id, 
                           style=wx.RESIZE_BORDER|wx.CLOSE_BOX, **kwargs)
	#
	#  Data Objects.
	#
	self.mat = mat
	#
	#  Windows.
	#
	self.titlebar = wx.StaticText(self, -1, 'PlaneDataDialog', 
				      style=wx.ALIGN_CENTER|wx.SIMPLE_BORDER)
        self.pdPanel = PlaneDataPanel(self, wx.NewId(), self.mat)
        self.quitBut  = wx.Button(self, wx.NewId(), 'QUIT')
	#
	#  Bindings.
	#
	self.makeBindings()
	#
	#  Sizing.
	#
	self.makeSizers()
	#
	self.SetAutoLayout(True)
        self.SetSizerAndFit(self.sizer)
	#
	#  Events.
	#
	
	#
	return

    def makeBindings(self):
	"""Bind interactors to functions"""
        self.Bind(wx.EVT_BUTTON, self.quit, self.quitBut)
        
	return

    def makeSizers(self):
	"""Lay out windows"""
	self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.sizer.Add(self.quitBut, 0, wx.EXPAND|wx.ALIGN_CENTER)
	self.sizer.Add(self.titlebar, 0, wx.EXPAND|wx.ALIGN_CENTER)
	self.sizer.Add(self.pdPanel,  1, wx.EXPAND|wx.ALIGN_CENTER)

	return

    def quit(self,e):
        """Quit the dialogue"""
        self.Destroy()
        return
    
    pass # end class
#
# -----------------------------------------------END CLASS:  PlaneDataDialog
#

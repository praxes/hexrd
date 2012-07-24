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
"""List Editor

Basic panel to  be used for managing generic
lists, including reordering and deletion of items.
"""
import copy

import wx

from hexrd.wx.guiconfig import WindowParameters as WP

from hexrd.xrd.experiment import newName  # need better place for newName function
#
# ---------------------------------------------------CLASS:  ListEditor
#
class ListEditor(wx.Panel):
    """ListEditor """
    def __init__(self, parent, id, mylist, **kwargs):
	"""Constructor for ListEditor."""
	#
	wx.Panel.__init__(self, parent, id, **kwargs)
	#
        #  Data
        #
        self.mylist = mylist
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

        #self.__makeTitleBar('List Editor')
        self.main_lbx =  wx.ListBox(self, wx.NewId(),
                                    style = wx.LB_SINGLE,
                                    choices = [item.name for item in self.mylist])
        self.up_but   = wx.Button(self, wx.NewId(), 'up')
        self.down_but = wx.Button(self, wx.NewId(), 'down')
        self.del_but  = wx.Button(self, wx.NewId(), 'del')
        self.copy_but = wx.Button(self, wx.NewId(), 'copy')

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
        self.Bind(wx.EVT_BUTTON, self.OnUpButton,   self.up_but)
        self.Bind(wx.EVT_BUTTON, self.OnDownButton, self.down_but)
        self.Bind(wx.EVT_BUTTON, self.OnDelButton,  self.del_but)
        self.Bind(wx.EVT_BUTTON, self.OnCopyButton, self.copy_but)

        return

    def __makeSizers(self):
	"""Lay out the interactors"""
        #
        #  Button sizer
        #
        butAlign = wx.ALIGN_LEFT
	self.butSizer = wx.BoxSizer(wx.VERTICAL)

	self.butSizer.Add(self.up_but,   0, butAlign)
	self.butSizer.Add(self.down_but, 0, butAlign)
	self.butSizer.Add(self.del_but,  0, butAlign)
	self.butSizer.Add(self.copy_but, 0, butAlign)
	#
        #  Main sizer
        #
	nrow = 1; ncol = 2; padx = 5; pady = 5
	self.sizer = wx.FlexGridSizer(nrow, ncol, padx, pady)
	self.sizer.AddGrowableCol(0, 1)
	self.sizer.AddGrowableRow(0, 1)

        self.sizer.Add(self.main_lbx,  1, wx.EXPAND | wx.ALIGN_CENTER | wx.RIGHT, 5)
        self.sizer.Add(self.butSizer,  0)

	return

    def _resetLbox(self):
        """Reset entries in listbox to match list"""
        self.main_lbx.Set([i.name for i in self.mylist])

        return

    def _listSwap(self, i, j):
        """Swap two elements of mylist"""
        tmp = self.mylist[i]
        self.mylist[i] = self.mylist[j]
        self.mylist[j] = tmp

        return
    #
    # ============================== API
    #
    def OnCopyButton(self, e):
        """Copy button pressed"""
        sel = self.main_lbx.GetSelection()

        if sel < 0: return  # none selected

        selItem = self.mylist[sel]
        newItem = copy.deepcopy(selItem)
        newItem.name = newName(newItem.name, [i.name for i in self.mylist])
        newSel = sel + 1
        self.mylist.insert(newSel, newItem)
        self._resetLbox()
        self.main_lbx.SetSelection(newSel)

        return

    def OnUpButton(self, e):
        """Up button pressed"""
        sel = self.main_lbx.GetSelection()
        if sel < 0: return

        if sel > 0:
            self._listSwap(sel, sel - 1)
            self._resetLbox()
            self.main_lbx.SetSelection(sel - 1)
            pass

        return

    def OnDownButton(self, e):
        """Up button pressed"""
        sel = self.main_lbx.GetSelection()
        if sel < 0: return

        if sel < len(self.mylist) - 1:
            self._listSwap(sel, sel + 1)
            self._resetLbox()
            self.main_lbx.SetSelection(sel + 1)
            pass

        return

    def OnDelButton(self, e):
        """Up button pressed"""
        sel = self.main_lbx.GetSelection()
        if sel < 0: return # no selection

        del self.mylist[sel]

        self._resetLbox()

        if sel > 0:
            newsel = sel - 1
        else:
            newsel = 0
            pass

        self.main_lbx.SetSelection(newsel)

        return

    pass # end class
#
# -----------------------------------------------END CLASS:  ListEditor
# ---------------------------------------------------CLASS:  ListEditDlg
#
class ListEditDlg(wx.Dialog):
    """ListEditDlg """
    def __init__(self, parent, id, mylist, **kwargs):
	"""Constructor for ListEditDlg"""
	#
	wx.Dialog.__init__(self, parent, id,
                           style=wx.RESIZE_BORDER|wx.DEFAULT_DIALOG_STYLE)
	#
	#  Data Objects.
	#
	#
	#  Windows.
	#
	self.titlebar = wx.StaticText(self, -1, 'List Editor',
				      style=wx.ALIGN_CENTER|wx.SIMPLE_BORDER)
        #
	self.list_ed = ListEditor(self, wx.NewId(), mylist)
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
	return

    def makeSizers(self):
	"""Lay out windows"""
	self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.sizer.Add(self.titlebar, 0, wx.EXPAND|wx.ALIGN_CENTER)
	self.sizer.Add(self.list_ed,  1, wx.EXPAND|wx.ALIGN_CENTER)

	return

    pass # end class
#
# -----------------------------------------------END CLASS:  ListEditDlg
#

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
"""wx Utilities
"""
import wx.html
#
from hexrd.wx.guiconfig import WindowParameters as WP, onLinux
#
def ResetChoice(cho, names, sname):
    """Reset the choice interactor.

    Surprisingly, there is no wx function for this.

    INPUTS
    . cho   -- the choice control
    . names -- the list of strings to populate it with
    . sname  --(string) the selected name

    OUTPUTS (NONE)

    DESCRIPTION
    . Resets the choice list and sets the selection.
"""
    # Reset strings
    cho.Clear()
    sel = 0
    for i in range(len(names)):
        n = names[i]
        cho.Append(names[i])
        if sname == n:
            sel = i
            pass
        pass

    cho.SetSelection(sel)

    return


def AddSpacer(p, sizer, color):
    """Add a spacer to a sizer

    p     - parent window
    sizer - where to put the spacer
    color - color to use
"""
    #  *** convert to wx.StaticLine()
    pad   = 10
    vsize = 5
    #  Create a sizer for this empty window
    sep = wx.BoxSizer(wx.HORIZONTAL)
    win = wx.Window(p, -1, size=(0,vsize))
    win.SetBackgroundColour(color)
    sep.Add(win, 1, wx.EXPAND|wx.TOP)
    #
    sizer.Add(sep, 0, wx.EXPAND|wx.ALL|wx.ALIGN_CENTER, pad)
    #
    return

def EmptyWindow(p): return wx.Window(p, -1)

def makeTitleBar(p, t, **kwargs):
    """Add titlebar

    INPUTS
    p - parent window
    t - title

    --------- Keyword Args
    tooltip - (str) tool tip string
    color   - (color) background color for titlebar

    OUTPUTS
    tsizer - titlebar sizer

    NOTES

    Static text alignment and color doesn't seem to work in linux.
    We use a workaround by creating a sizer with colored boxes
    on either side.
"""
    titlebar = wx.StaticText(p, wx.NewIdRef(), t,
                             style=wx.ALIGN_CENTER|wx.SIMPLE_BORDER)
    #
    #  Keyword args
    #
    tt = 'tooltip'
    if tt in kwargs:  titlebar.SetToolTip(kwargs[tt])

    cl = 'color'
    if cl in kwargs:
        bgColor = kwargs[cl]
    else:
        bgColor = WP.BG_COLOR_PANEL_TITLEBAR
        pass
    #
    #  Create sizer
    #
    tsizer = wx.BoxSizer(wx.HORIZONTAL)
    if onLinux:
        w1 = wx.StaticLine(p, -1, size=(5,5))
        w1.SetBackgroundColour(bgColor)
        w2 = wx.StaticLine(p, -1, size=(5,5))
        w2.SetBackgroundColour(bgColor)

        tsizer.Add(w1, 1, wx.EXPAND)
        tsizer.Add(titlebar, 0, wx.ALIGN_CENTER)
        #                  10 pixel padding on bottom
        tsizer.Add(w2, 1, wx.EXPAND)
    else:
        titlebar.SetBackgroundColour(bgColor)
        tsizer.Add(titlebar, 1, wx.EXPAND|wx.ALIGN_CENTER)
        pass

    return tsizer

def callJoel(p):
    """Return message to display on empty pages"""
    hpage = wx.html.HtmlWindow(p, wx.NewIdRef())

    msg = r"""<html>
<body>
<h1> NOT YET IMPLEMENTED
</body>
</html>
"""
    hpage.SetPage(msg)

    return hpage

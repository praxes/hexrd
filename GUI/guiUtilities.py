#! /usr/bin/env python
#
"""GUI Utilities
"""
import wx.html
#
from guiConfig import WindowParameters as WP, onLinux
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
    titlebar = wx.StaticText(p, wx.NewId(), t, 
                             style=wx.ALIGN_CENTER|wx.SIMPLE_BORDER) 
    #
    #  Keyword args
    #
    tt = 'tooltip'
    if tt in kwargs:  titlebar.SetToolTipString(kwargs[tt])

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
    hpage = wx.html.HtmlWindow(p, wx.NewId())

    msg = r"""<html>
<body>
<h1> NOT YET IMPLEMENTED
</body>
</html>
"""
    hpage.SetPage(msg)
    
    return hpage

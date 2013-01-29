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
"""Main visualization window for hydra
"""
import wx
#
#  Matplotlib stuff
#
import matplotlib;  matplotlib.use('WXAgg')

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg
from matplotlib.figure                 import Figure
from matplotlib                        import pyplot
from matplotlib                        import cm
#
#  XRD package
#
from hexrd.wx.guiconfig  import WindowParameters as WP
from hexrd.wx.guiutil import ResetChoice

from hexrd.wx.floatcontrol import *
#
# ---------------------------------------------------CLASS:  hydraCanvasPanel
#
class hydraCanvasPanel(wx.Panel):
    """hydraCanvasPanel """
    def __init__(self, parent, id, **kwargs):
        """Constructor for hydraCanvasPanel."""
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

        self.__makeTitleBar('hydraCanvasPanel')
        self.SetBackgroundColour(WP.CANVAS_BG_COLOR)

        #
        #  ===== FIGURE CANVAS
        #
        #  Create figure, making room for 4 axes.
        #
        self.figure = Figure()
        #
        #  * rectangles are (l, b, w, h)
        #
        self.axes = 4*[None]
        rectWidth = 0.4; rectHeight = 0.4
        self.axes[0] = self.figure.add_axes((0.05, 0.55, rectWidth, rectHeight))
        self.axes[1] = self.figure.add_axes((0.55, 0.55, rectWidth, rectHeight))
        self.axes[2] = self.figure.add_axes((0.05, 0.05, rectWidth, rectHeight))
        self.axes[3] = self.figure.add_axes((0.55, 0.05, rectWidth, rectHeight))
        #self.axes.set_aspect('equal')
        #
        self.canvas = FigureCanvas(self, wx.NewId(), self.figure)
        self.__add_toolbar()  # comment this out for no toolbar

        def on_press(event):
            myFrame = self.GetParent()
            if hasattr(event, 'xdata') and event.xdata:
                myFrame.SetStatusText('Hi from %g, %g' % (event.xdata, event.ydata))
                pass
            pass

        #cid = self.canvas.mpl_connect('button_press_event', on_press)
        cid = self.canvas.mpl_connect('motion_notify_event', on_press)

        return

    def __add_toolbar(self):
        """add toolbar"""
        self.toolbar = NavigationToolbar2WxAgg(self.canvas)
        self.toolbar.Realize()
        self.toolbar.update()

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
        return

    def __makeSizers(self):
        """Lay out the interactors"""

        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.titlebar, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.sizer.Add(self.canvas,   1, wx.LEFT | wx.TOP | wx.GROW)
        self.sizer.Add(self.toolbar,  0, wx.LEFT | wx.EXPAND)

        return
    #
    # ============================== API
    #
    #                     ========== *** Access Methods
    #
    def loadImages(self):
        """Load images into axes"""
        exp = wx.GetApp().ws
        h   = exp.hydra
        for i in range(4):
            self.axes[i].imshow(h.images[i],
                                origin='upper',
                                interpolation='nearest',
                                cmap=getattr(cm, 'bone'))
            self.axes[i].set_autoscale_on(False)
            self.axes[i].format_coord =  lambda x, y: str(x) + str(y)
            pass

        return
    #
    #                     ========== *** Event Callbacks
    #

    pass # end class
#
# -----------------------------------------------END CLASS:  hydraCanvasPanel
# ---------------------------------------------------CLASS:  hydraCanvasFrame
#
class hydraCanvasFrame(wx.Frame):
    #
    def __init__(self, parent, id, title='Hydra Canvas'):
        #
        #  Pass option dictionary or string.
        #
        wx.Frame.__init__(self, parent, id, title)
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
        self.Show(True)

        return

    pass # class
    #
    # ============================== Internal Methods
    #
    def __makeObjects(self):
        """Add interactors"""

        self.__makeTitleBar('Hydra Canvas')
        #
        # A Statusbar in the bottom of the window
        #
        self.CreateStatusBar()
        #
        # Create the panel.
        #
        self.hcanvas = hydraCanvasPanel(self, wx.NewId())
        #
        # Creating the menubar.
        #
        menuBar = wx.MenuBar()
        # self.CreateFileMenu()
        # menuBar.Append(self.filemenu,  "&File")
        # self.CreateTableMenu()
        # menuBar.Append(self.tablemenu, "&Table")
        #
        # self.SetMenuBar(menuBar)  # Adding the MenuBar to the Frame content.

        return

    def __makeTitleBar(self, t):
        """Add titlebar"""
        self.titlebar = wx.StaticText(self, -1, t,
                                         style=wx.ALIGN_CENTER|wx.SIMPLE_BORDER)
        self.titlebar.SetBackgroundColour(WP.TITLEBAR_BG_COLOR_FRAME)
        myToolTip = r"""
FRAME FOR ...
"""
        self.titlebar.SetToolTipString(myToolTip)

        return

    def __makeBindings(self):
        """Bind interactors"""
        return

    def __makeSizers(self):
        """Lay out the interactors"""

        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.titlebar, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.sizer.Add(self.hcanvas, 1, wx.EXPAND|wx.ALIGN_CENTER)

        return
    #
    # ============================== API
    #
    #                     ========== *** Access Methods
    #

    #
    #                     ========== *** Event Callbacks
    #
    pass
#
# -----------------------------------------------END CLASS:  hydraCanvasFrame
#

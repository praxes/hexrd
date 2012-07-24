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
"""Canvas Panel for matplotlib graphics

We will have two subpanels: an options panel at the
top and a graphics canvas below it.
"""
import wx

import numpy
#
#  XRD package
#
from hexrd.wx.guiconfig import WindowParameters as WP
from hexrd.wx.guiutil import ResetChoice,makeTitleBar
from hexrd.wx.canvasutil import *

from hexrd.wx.floatcontrol import *
#
#  Data
#
r2d = 180.0/numpy.pi
#
# ---------------------------------------------------CLASS:  CanvasPanel
#
class CanvasPanel(wx.Panel):
    """CanvasPanel """
    def __init__(self, parent, id, **kwargs):
	"""Constructor for CanvasPanel."""
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
        self.tbarSizer = makeTitleBar(self, 'Graphics Canvas')
        self.SetBackgroundColour(WP.CANVAS_BG_COLOR)
        #
        #  Make options sizer
        #
        nrow = 0; ncol = 2; padx = 5; pady = 5
	self.optSizer = wx.FlexGridSizer(nrow, ncol, padx, pady)
	#self.optSizer.AddGrowableCol(num, proportion)
	#self.optSizer.AddGrowableRow(num, proportion)
        #
        #  ===== OPTIONS
        #
        #  * show image
        #
        self.showImage_box = wx.CheckBox(self, wx.NewId(),
                                         'Show Image')
        self.showImage_box.SetValue(True)
        self.optSizer.Add(self.showImage_box, 0, wx.LEFT | wx.EXPAND)
        self.optSizer.AddSpacer(1)
        #
        #  * show rings
        #
        self.showCalRings_box = wx.CheckBox(self, wx.NewId(),
                                            'Show Rings')
        self.showCalRings_box.SetValue(False) # default
        self.optSizer.Add(self.showCalRings_box, 0, wx.LEFT | wx.EXPAND)
        self.optSizer.AddSpacer(1)
        #
        #  * show ranges
        #
        self.showCalRanges_box = wx.CheckBox(self, wx.NewId(),
                                            'Show Ranges')
        self.showCalRanges_box.SetValue(False) # default
        self.optSizer.Add(self.showCalRanges_box, 0, wx.LEFT | wx.EXPAND)
        self.optSizer.AddSpacer(1)
        #
        #  Add colormap panel
        #
        self.cmPanel = cmapPanel(self, wx.NewId())
        #
        #  ===== FIGURE CANVAS
        #
        self.figure = Figure()
        self.axes   = self.figure.gca()
        self.axes.set_aspect('equal')
        self.canvas = FigureCanvas(self, wx.NewId(), self.figure)
        self.__add_toolbar()  # comment this out for no toolbar


        def on_press(event):
            exp = wx.GetApp().ws
            det = exp.detector
            pd  = exp.activeMaterial.planeData
            img = exp.activeImage
            if img is None:  return

            mainFrame = wx.GetApp().GetTopWindow()
            if hasattr(event, 'xdata') and event.xdata:
                x = event.xdata; xadj = x + 0.5; xint = numpy.floor(xadj)
                y = event.ydata; yadj = y + 0.5; yint = numpy.floor(yadj)
                tth, eta = numpy.array(det.xyoToAng_V(y, x))
                cartx, carty = det.cartesianCoordsOfPixelIndices(y, x)
                cx = (cartx - det.xc)/det.pixelPitch
                cy = (carty - det.yc)/det.pixelPitch
                rho = numpy.sqrt(cx*cx + cy*cy)
                dsp = 0.5*pd.wavelength/numpy.sin(0.5*tth)
                intens = img[yadj, xadj]
                hkls = str(pd.getHKLs(asStr=True, allHKLs=True, thisTTh=tth))
                statText = "px=%g, py=%g, x=%g, y=%g, rho=%g, d=%g, tth=%g, eta=%g, int=%g, HKLs=%s" %\
                           (x, y, cartx, carty, rho, dsp, r2d*tth, r2d*eta, intens, hkls)

                mainFrame.SetStatusText(statText)
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

    def __makeBindings(self):
        """Bind interactors"""

        self.Bind(wx.EVT_CHECKBOX, self.OnCheckImage,     self.showImage_box)
        self.Bind(wx.EVT_CHECKBOX, self.OnCheckCalRings,  self.showCalRings_box)
        self.Bind(wx.EVT_CHECKBOX, self.OnCheckCalRanges, self.showCalRanges_box)
        return

    def __makeSizers(self):
	"""Lay out the interactors"""

	self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.sizer.Add(self.tbarSizer, 0,
                       wx.EXPAND| wx.BOTTOM, 10)
        self.sizer.Add(self.optSizer, 0, wx.LEFT | wx.TOP | wx.GROW)
        self.sizer.Add(self.cmPanel,  0, wx.LEFT | wx.TOP | wx.GROW)
        self.sizer.Add(self.canvas,   1, wx.LEFT | wx.TOP | wx.GROW)
        self.sizer.Add(self.toolbar,  0, wx.LEFT | wx.EXPAND)

	return
    #
    #  Helper functions
    #
    def __drawRings(self, mat=None):
        """Draw rings for a particular material"""
        exp = wx.GetApp().ws
        #
        #  Material
        #
        if mat is None:
            pdat = exp.activeMaterial.planeData
            #pdat.tthMax = exp.detector.getTThMax()
        else:
            pdat = mat.planeData
            pass

        xyRings = exp.detector.getRings(pdat)
        opts = dict(linewidth=2, color='g', linestyle='-')
        self.addXYplot(xyRings, opts=opts)

        return

    def __drawRanges(self, mat=None):
        """Draw rings for a particular material"""
        exp = wx.GetApp().ws
        #
        #  Material
        #
        if mat is None:
            pdat = exp.activeMaterial.planeData
            #pdat.tthMax = exp.detector.getTThMax()
        else:
            pdat = mat.planeData
            pass

        xyRings = exp.detector.getRings(pdat, ranges=True)
        opts = dict(linewidth=2, color='y', linestyle='--')
        self.addXYplot(xyRings, opts=opts)

        return

    def __clearAxesLines(self):
        """Remove lines from axes object"""
        for l in self.axes.get_lines():
            l.remove()
            pass

        return
    #
    # ============================== API
    #
    #                     ========== *** Utility Methods
    def update(self, **kwargs):
        """Update the image according to the options in the option panel

        KEYWORD ARGS
        newImage -- if True then display image for first time on these axes
"""
        #
        #  Show image if box is checked.
        #
        app = wx.GetApp()
        exp = app.ws
        img = exp.activeImage

        if img is None:
            #wx.MessageBox('no image'); return
            pass
        else:
            if 'newImage' in kwargs:
                ni = kwargs['newImage']
            else:
                ni = False
                pass

            si = self.showImage_box.IsChecked()

            if ni:
                # delete old image first
                # print 'deleting old images'
                self.axes.images = []
                # show new image
                self.axes.imshow(img,
                                 origin='upper',
                                 interpolation='nearest',
                                 cmap=self.cmPanel.cmap,
                                 vmin=self.cmPanel.cmin_val,
                                 vmax=self.cmPanel.cmax_val,
                                 visible=si)
                self.axes.set_autoscale_on(False)
                self.axes.format_coord =  lambda x, y: str(x) + str(y)
                # print 'number of images in current axes:  %d' % len(self.axes.get_images())
            else:
                # set visibility of axes image
                img0 = self.axes.get_images()[0]
                img0.set_visible(si)
                pass
            pass
        #
        #  Show calibrant rings/ranges if box is checked.
        #
        self.__clearAxesLines()

        if (self.showCalRings_box.IsChecked()):  self.__drawRings()
        if (self.showCalRanges_box.IsChecked()):  self.__drawRanges()
        #
        #  Update ring list
        #
        #rcho = self.rings_cho
        #ResetChoice(rcho, exp.matNames, rcho.GetStringSelection)
        #
        self.draw()

        return

    def draw(self): self.canvas.draw()

    def clear(self):
        """Clear axes (including imag"""

        self.axes.cla()
        self.draw()

        return

    def clearLines(self):
        """Remove line objects from current axes"""
        self.__clearAxesLines()
        return

    def addXYplot(self, rings, opts=None):
        """Add axes to figure

        INPUTS
        rings - is a list of 2-tuples, with first being the x, second the y
        opts  - a dictionary of plot options, which gets passed to axes.plot()
              keys/default values are:
              . linewidth:  3
              . linestyle:  ':'
              .     color:  'g'
"""
        #  NOTE:  image coordinates are y, x instead of x, y
        colorS = 'g'
        if opts is None:
            optDict = dict(linewidth=3, linestyle=':')
        else:
            optDict = opts
            if 'color' in opts:
                colorS = opts['color']
                pass

            pass

        for xy in rings:
            self.axes.plot(xy[1], xy[0], colorS, **optDict)
            pass

        self.draw()

        return
    #                     ========== *** Event Callbacks

    def OnCheckCalRings(self, evt):
        """Callback for calRings_box"""
        self.update()

        return

    def OnCheckCalRanges(self, evt):
        """Callback for calRanges_box"""
        self.update()

        return

    def OnCheckImage(self, evt):
        """Callback for showImage_box"""
        self.update()

        return

    def showImage(self, img):
        """Show image on current axes"""
        self.axes.imshow(img)
        self.draw()

        return

    pass # end class
#
# -----------------------------------------------END CLASS:  CanvasPanel

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
from hexrd.wx.guiutil import ResetChoice, makeTitleBar, EmptyWindow
from hexrd.wx.listeditor import NamedItem, ListEditDlg
from hexrd.wx.canvasutil import FigureCanvas, NavigationToolbar2WxAgg,\
     Figure, cmapPanel
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
        #  Add image list management
        #
        self.ail_lab = wx.StaticText(self, wx.NewId(), 'Load Image', style=wx.ALIGN_CENTER)
        self.ail_cho = wx.Choice(self, wx.NewId(), choices=[])
        self.nam_lab = wx.StaticText(self, wx.NewId(), 'Name Image', style=wx.ALIGN_CENTER)
        self.nam_txt = wx.TextCtrl(self, wx.NewId(), value='<none>',
                                    style=wx.RAISED_BORDER|wx.TE_PROCESS_ENTER)
        self.eil_but  = wx.Button(self, wx.NewId(), 'Edit List')
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
            img = exp.active_img
            if img is None:  return

            mainFrame = wx.GetApp().GetTopWindow()
            if hasattr(event, 'xdata') and event.xdata:
                x = event.xdata; xadj = x + 0.5; xint = numpy.floor(xadj)
                y = event.ydata; yadj = y + 0.5; yint = numpy.floor(yadj)
                tth, eta = numpy.array(det.xyoToAng(y, x))
                cartx, carty = det.cartesianCoordsOfPixelIndices(y, x)
                cx = (cartx - det.xc)/det.pixelPitch
                cy = (carty - det.yc)/det.pixelPitch
                rho = numpy.sqrt(cx*cx + cy*cy)
                dsp = 0.5*pd.wavelength/numpy.sin(0.5*tth)
                intens = img[yint, xint]
                hkls = str(pd.getHKLs(asStr=True, allHKLs=True, thisTTh=tth))
                statText = "px=%g, py=%g, x=%g, y=%g, rho=%g, d=%g, "\
                           "tth=%g, eta=%g, int=%g, HKLs=%s" %\
                           (x, y, cartx, carty, rho, dsp,
                            r2d*tth, r2d*eta, intens, hkls)

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

        self.Bind(wx.EVT_TEXT_ENTER, self.OnNameImg, self.nam_txt)
        self.Bind(wx.EVT_CHOICE, self.OnLoadImg, self.ail_cho)
        self.Bind(wx.EVT_BUTTON, self.OnEditImg, self.eil_but)


        return

    def __makeSizers(self):
        """Lay out the interactors"""
        nrow = 0; ncol = 2; padx = pady = 5
        self.ilsizer = wx.FlexGridSizer(nrow, ncol, padx, pady)
        self.ilsizer.AddGrowableCol(1, 1)

        self.ilsizer.Add(self.ail_lab,  0, wx.EXPAND|wx.ALIGN_CENTER|wx.TOP, 5)
        self.ilsizer.Add(self.ail_cho,  1, wx.EXPAND|wx.ALIGN_CENTER|wx.TOP, 5)
        self.ilsizer.Add(self.nam_lab,  0, wx.EXPAND|wx.ALIGN_CENTER|wx.TOP, 5)
        self.ilsizer.Add(self.nam_txt,  1, wx.EXPAND|wx.ALIGN_CENTER|wx.TOP, 5)
        self.ilsizer.Add(EmptyWindow(self),  0, wx.ALIGN_CENTER|wx.TOP, 5)
        self.ilsizer.Add(self.eil_but,  1, wx.ALIGN_CENTER|wx.TOP, 5)

        self.osizer = wx.BoxSizer(wx.HORIZONTAL)
        self.osizer.Add(self.optSizer, 0)
        self.osizer.Add(self.ilsizer, 1)

        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.tbarSizer, 0,
                       wx.EXPAND| wx.BOTTOM, 10)
        self.sizer.Add(self.osizer,   0, wx.LEFT | wx.TOP | wx.GROW)
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
        newImage -- if True then display image afresh (False by default)
        loadImage -- if True, then image was loaded from saved list (False by default)
"""
        kwargs.setdefault('newImage', False)
        kwargs.setdefault('loadImage', False)
        kwargs.setdefault('updateImage', False)
        kwargs.setdefault('onInit', False)

        ni = kwargs['newImage']
        li = kwargs['loadImage']
        ui = kwargs['updateImage']
        oninit = kwargs['onInit']
        print 'li: ', li
        #
        #  Show image if box is checked.
        #
        app = wx.GetApp()
        exp = app.ws
        img = exp.active_img
        

        if img is None:
            # no active image, but possibly one on the axes
            images = self.axes.get_images()
            if images:
                img0 = images[0]
                img0.set_visible(False)
        else:
            si = self.showImage_box.IsChecked()
            
            if ni or ui:
                # not using axes image list
                if ni: self.axes.set_autoscale_on(True)
        
                self.axes.images = []

                self.axes.imshow(img,
                                 origin='upper',
                                 interpolation='nearest',
                                 cmap=self.cmPanel.cmap,
                                 vmin=self.cmPanel.cmin_val,
                                 vmax=self.cmPanel.cmax_val,
                                 visible=si)
                self.axes.set_autoscale_on(False)
                self.axes.format_coord =  lambda x, y: str(x) + str(y)
                #
            else:
                # set visibility of axes image
                images = self.axes.get_images()
                if images:
                    img0 = images[0]
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
        if li: 
            print 'loading image: working with axes'
            self.axes.set_autoscale_on(True)
      
        self.draw()
        #
        # Update image list
        #
        acho = self.ail_cho
        if li:
            # set text box to name of interactors
            name = acho.GetStringSelection()
            ResetChoice(acho, exp.img_names, name)
            self.nam_txt.ChangeValue(name)
        else:
            if ni or oninit:
                if exp.img_names: # to handle init case on load exp
                    ResetChoice(acho, exp.img_names, exp.img_names[0])
                    acho.SetSelection(wx.NOT_FOUND)
                    self.nam_txt.ChangeValue('<unnamed image>')

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
    def OnEditImg(self, evt):
        """Edit image list"""
        exp = wx.GetApp().ws

        #
        # Since images do not have a name attribute, use NamedItem class
        # to make list with names
        #
        ilist = exp.img_list
        iname = exp.img_names
        nilist = [NamedItem(iname[i], ilist[i]) for i in range(len(iname))]

        ssel = self.ail_cho.GetStringSelection()

        dlg = ListEditDlg(self, wx.NewId(), nilist)
        dlg.ShowModal()
        dlg.Destroy()

        exp.img_list[:] = [item.data for item in nilist]
        exp.img_names[:] = [item.name for item in nilist]

        #
        # Now reset choice by hand in case current image was dropped
        #
        if (not ssel):
            # new names, but no selection
            ResetChoice(self.ail_cho, exp.img_names, ssel)
            self.ail_cho.SetSelection(wx.NOT_FOUND)
            li = False
        elif (ssel in exp.img_names):
            # new names, keep old selection
            exp.active_img = ssel
            ResetChoice(self.ail_cho, exp.img_names, ssel)
            li = True
        else:
            # new names, old selection gone
            exp.active_img = None
            ResetChoice(self.ail_cho, exp.img_names, '')
            self.ail_cho.SetSelection(wx.NOT_FOUND)
            li = False

        ni = not li
        self.update(newImage=ni, loadImage=li)

        return

    def OnLoadImg(self, evt):
        """Load image from list"""
        exp = wx.GetApp().ws
        exp.active_img = evt.GetSelection()

        self.update(loadImage=True, newImage=True)

        return

    def OnNameImg(self, evt):
        """Name the curent Image"""
        exp = wx.GetApp().ws

        name = evt.GetString()
        acho = self.ail_cho
        achosel = acho.GetSelection()
        if achosel == wx.NOT_FOUND:
            exp.add_to_img_list(name)
            ResetChoice(acho, exp.img_names, name)
            self.update(loadImage=True, newImage=True)
        else:
            exp.img_names[achosel] = name
            ResetChoice(acho, exp.img_names, name)

        return

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

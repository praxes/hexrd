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
"""Display for caking output.
"""
import wx

import numpy
from matplotlib.axes import Axes
from matplotlib.patches import Rectangle, Circle, Polygon
from matplotlib.collections import PatchCollection

from hexrd.matrixutil import columnNorm, unitVector

from hexrd.xrd.experiment  import PolarRebinOpts as prOpts
from hexrd.xrd.xrdutil import makeMeasuredScatteringVectors as makeMSV

from hexrd.wx.guiconfig import WindowParameters as WP
from hexrd.wx.guiutil import makeTitleBar
from hexrd.wx.canvasutil import *
#
# ---------------------------------------------------CLASS:  cakeCanvas
#
class cakeCanvas(wx.Panel):
    """cakeCanvas """
    def __init__(self, parent, id, cakeType, data, **kwargs):
	"""Constructor for cakeCanvas panel

        **INPUTS**
        *cakeType* - flag for type of analysis: full image caking, ring-based or spherical (omega-eta map)
        *data* - the object with data to be plotted caking, depending on *cakeType*
                 .. for rings case, a MultiRingBinned instance
                 .. for full image case, img_info (output of polarRebin)

        
        **OUTPUTS**
        *self* - this panel

        
        **DESCRIPTION**
        Shows data on output canvas.
"""
	#
	wx.Panel.__init__(self, parent, id, **kwargs)
	#
        #  Data
        #
        self.cakeType = cakeType
        self.data     = data
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
        self.opt_pan.update() # make initial plot
 	#
	return
    #
    # ============================== Internal Methods
    #
    def __makeObjects(self):
        """Add interactors"""

        self.tbarSizer = makeTitleBar(self, 'cakeCanvas')
        #
        if   self.cakeType == prOpts.CAKE_IMG:
            self.opt_pan = imgOpts(self, wx.NewId())
            # full image caking panel
            pass
        elif self.cakeType == prOpts.CAKE_RNG:
            # multiring panel
            self.opt_pan = rngOpts(self, wx.NewId())
            pass
        elif self.cakeType == prOpts.CAKE_SPH:
            # omega-eta panel
            self.opt_pan = sphOpts(self, wx.NewId())
            pass
        
        self._makeFigureCanvas()
        #

        return

    def _makeFigureCanvas(self):
        """Build figure canvas"""
        self.figure = Figure()
        self.canvas = FigureCanvas(self, wx.NewId(), self.figure)

        self.axes = self.figure.gca()
        self.axes.set_aspect('equal')

        self.toolbar = NavigationToolbar2WxAgg(self.canvas)
        self.toolbar.Realize()
        self.toolbar.update()

        return

    def __makeBindings(self):
        """Bind interactors"""
        return

    def __makeSizers(self):
	"""Lay out the interactors"""

	self.sizer = wx.BoxSizer(wx.VERTICAL)
        #
	self.sizer.Add(self.tbarSizer, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.sizer.Show(self.tbarSizer, False)
	self.sizer.Add(self.opt_pan,   0, wx.EXPAND|wx.ALIGN_CENTER)
	self.sizer.Add(self.canvas,    1, wx.EXPAND|wx.ALIGN_CENTER)
	self.sizer.Add(self.toolbar,   0, wx.EXPAND|wx.ALIGN_CENTER)

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
# -----------------------------------------------END CLASS:  cakeCanvas
# ---------------------------------------------------CLASS:  cakeDisplay
#
class cakeDisplay(wx.Frame):
    #
    def __init__(self, parent, id, cakeType, data, title=''):
        """
        INPUTS
        cakeType -- full caking, ring-based or spherical (omega-eta map)
        data     -- caking output, depending on cakeType

        OUTPUTS

        DESCRIPTION
        Passes args to cakeCanvas
"""
        #
        #  Pass option dictionary or string.
        #
        wx.Frame.__init__(self, parent, id, title)
	#
        #  Data
        #
        self.cakeType = cakeType
        self.data     = data
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
    #
    # ============================== Internal Methods
    #
    def __makeObjects(self):
        """Add interactors"""

        self.tbarSizer = makeTitleBar(self, 'Rebin Canvas')
        #
        #  Add canvas panel
        #
        self.cpan = cakeCanvas(self, wx.NewId(),  self.cakeType, self.data)
	#
        # A Statusbar in the bottom of the window
        #
        self.CreateStatusBar()
	#
        # Creating the menubar.
        #
        # menuBar = wx.MenuBar()
        # self.CreateFileMenu()
        # menuBar.Append(self.filemenu,  "&File")
        # self.CreateTableMenu()
        # menuBar.Append(self.tablemenu, "&Table")
        #
        # self.SetMenuBar(menuBar)  # Adding the MenuBar to the Frame content.

        return

    def __makeBindings(self):
        """Bind interactors"""
        return

    def __makeSizers(self):
	"""Lay out the interactors"""

	self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.sizer.Add(self.tbarSizer, 0, wx.EXPAND|wx.ALIGN_CENTER)
	self.sizer.Add(self.cpan,      1, wx.EXPAND|wx.ALIGN_CENTER)

	return
    #
    # ============================== API
    #
    #                     ========== *** Access Methods
    #

    #
    #                     ========== *** Event Callbacks
    #
    pass # class
#
# -----------------------------------------------END CLASS:
# ---------------------------------------------------CLASS:  rngOpts
#
class rngOpts(wx.Panel):
    """rngOpts """
    unitList = ['d-spacing', 'radians', 'strain', 'degrees']

    def __init__(self, parent, id, **kwargs):
	"""Constructor for rngOpts."""
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

	return
    #
    # ============================== Internal Methods
    #
    def __makeObjects(self):
        """Add interactors"""

        self.tbarSizer = makeTitleBar(self, 'Multiring Rebinning Results')
        #
        self.unit_cho = wx.Choice(self, wx.NewId(), choices=rngOpts.unitList)
        self.exp_but  = wx.Button(self, wx.NewId(), 'Export')
	#self.Bind(wx.EVT_CHOICE, self.OnChoice, self.choice)

        return

    def __makeBindings(self):
        """Bind interactors"""
        self.Bind(wx.EVT_CHOICE, self.OnUpdate, self.unit_cho)
        self.Bind(wx.EVT_BUTTON, self.OnExport, self.exp_but)
        return

    def __makeSizers(self):
	"""Lay out the interactors"""

	self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.sizer.Add(self.tbarSizer, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.sizer.Show(self.tbarSizer, False)
        self.sizer.Add(self.unit_cho, 0, wx.ALIGN_LEFT)
        self.sizer.Add(self.exp_but,  0, wx.ALIGN_LEFT|wx.TOP, 5)

	return
    #
    # ============================== API
    #
    #                     ========== *** Access Methods
    #
    def update(self):
        """Update canvas"""
        p = self.GetParent()
        u = self.unit_cho.GetStringSelection()

        self.errs = p.data.getTThErrors(units=u)

        p.figure.delaxes(p.axes)
        p.axes   = p.figure.gca()
        p.axes.plot(self.errs)
        p.axes.set_title('Errors')
        p.axes.set_ylabel(u)

        p.canvas.draw()

        return
    #
    #                     ========== *** Event Callbacks
    #
    def OnUpdate(self, e):
        """Update canvas"""
        self.update()
        return

    def OnExport(self, e):
        """Export results to a text file"""
        # export self.errs to a file
        dlg = wx.FileDialog(self, 'Export Binning Errors', style=wx.FD_SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            f = dlg.GetPath()
            try:
                # numpy.savetxt(f, self.errs)
                p = self.GetParent()
                p.data.getTThErrors(outputFile=f)
            except:
                wx.MessageBox('failed to save file:  %s' % f)
                pass
            pass


        dlg.Destroy()

        return

    pass # end class
#
# -----------------------------------------------END CLASS:  rngOpts
# ---------------------------------------------------CLASS:  imgOpts
#
class imgOpts(wx.Panel):
    """Options panel for IMG (standard) caking canvas"""

    def __init__(self, parent, id, **kwargs):
	"""Constructor for imgOpts."""
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

	return
    #
    # ============================== Internal Methods
    #
    def __makeObjects(self):
        """Add interactors"""
        self.tbarSizer = makeTitleBar(self, 'Full Image Rebinning Results')
        self.cmPanel = cmapPanel(self, wx.NewId())
        #
        return

    def __makeBindings(self):
        """Bind interactors"""
        return

    def __makeSizers(self):
	"""Lay out the interactors"""
	self.sizer = wx.BoxSizer(wx.HORIZONTAL)
	self.sizer.Add(self.tbarSizer, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.sizer.Show(self.tbarSizer, False)
	self.sizer.Add(self.cmPanel,   1, wx.EXPAND|wx.ALIGN_CENTER)        

	return
    #
    # ============================== API
    #
    #                     ========== *** Access Methods
    #
    def update(self, **kwargs):
        """Update canvas"""
        p   = self.GetParent()
        # tlp = wx.GetTopLevelParent(self)
        
        intensity = p.data['intensity']

        p.axes = p.figure.gca()
        p.axes.set_aspect('equal')


        p.axes.images = []
        # show new image
        p.axes.imshow(intensity, origin='upper',
                      interpolation='nearest',
                      aspect='auto',
                      cmap=self.cmPanel.cmap,
                      vmin=self.cmPanel.cmin_val,
                      vmax=self.cmPanel.cmax_val)
        p.axes.set_autoscale_on(False)

        p.canvas.draw()

        return
    #
    #                     ========== *** Event Callbacks
    #
    def OnUpdate(self, e):
        """Update canvas"""
        self.update()
        return

    def OnExport(self, e):
        """Export results to a text file"""
        # export self.errs to a file
        dlg = wx.FileDialog(self, 'Export Binning Data', style=wx.FD_SAVE)
        if dlg.ShowModal() == wx.ID_OK:
            f = dlg.GetPath()
            try:
                wx.MessageBox('would save to file:  %s' % f)
                #exp.saveDetector(f)
            except:
                wx.MessageBox('failed to save file:  %s' % f)
                pass
            pass

        dlg.Destroy()

        return

    pass # end class
#
# -----------------------------------------------END CLASS:  imgOpts
# ---------------------------------------------------CLASS:  sphOpts
#
class sphOpts(wx.Panel):
    """Options panel for SPH (omega-eta) caking canvas"""
    [DISP_RAW, DISP_QUICK, DISP_FULL] = range(3)
    DISP_METHODS = ['Raw', 'Quick Render', 'Full Render']
    
    def __init__(self, parent, id, **kwargs):
	"""Constructor for sphOpts."""
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

	return
    #
    # ============================== Internal Methods
    #
    def __makeObjects(self):
        """Add interactors"""
        exp = wx.GetApp().ws
        self.cmPanel = cmapPanel(self, wx.NewId())
        self.tbarSizer = makeTitleBar(self, 'Omega-Eta Plots',
                                      color=WP.BG_COLOR_TITLEBAR_PANEL1)

        # choice interactor for HKL
        hkls = exp.activeMaterial.planeData.getHKLs(asStr=True)
        self.hkl_cho = wx.Choice(self, wx.NewId(), choices=hkls)
        self.hkl_cho.SetSelection(0)

        self.disp_cho = wx.Choice(self, wx.NewId(), choices=self.DISP_METHODS)
        self.disp_cho.SetSelection(0)

        self.idata = 0
        self.dispm = self.DISP_RAW
        
        return

    def __makeBindings(self):
        """Bind interactors"""
	self.Bind(wx.EVT_CHOICE, self.OnHKLChoice, self.hkl_cho)
	self.Bind(wx.EVT_CHOICE, self.OnDispChoice, self.disp_cho)
        return

    def __makeSizers(self):
	"""Lay out the interactors"""
	self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.sizer.Add(self.tbarSizer, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.sizer.Show(self.tbarSizer, True)
        self.osizer = wx.BoxSizer(wx.VERTICAL)
        self.osizer.Add(self.hkl_cho,  1, wx.ALIGN_LEFT|wx.TOP, 5)
        self.osizer.Add(self.disp_cho, 1, wx.ALIGN_LEFT|wx.TOP, 5)
        self.csizer =wx.BoxSizer(wx.HORIZONTAL)
        self.csizer.Add(self.osizer, 1, wx.ALIGN_RIGHT|wx.TOP, 5)
        self.csizer.Add(self.cmPanel, 1, wx.ALIGN_LEFT|wx.TOP, 5)
        self.sizer.Add(self.csizer, 1, wx.ALIGN_CENTER|wx.EXPAND)
        
	return
    #
    # ============================== API
    #
    #                     ========== *** Access Methods
    #
    def update(self, **kwargs):
        """Update canvas"""
        p = self.GetParent()
        exp = wx.GetApp().ws
        
        ome_eta = p.data
        hkldata = ome_eta.getData(self.idata)

        if self.dispm == self.DISP_RAW:
            
            p.figure.delaxes(p.axes)
            p.axes = p.figure.gca()
            p.axes.set_aspect('equal')
            p.axes.images = []
            # show new image
            p.axes.imshow(hkldata, origin='upper',
                          interpolation='nearest',
                          aspect='auto', 
                          cmap=self.cmPanel.cmap,
                          vmin=self.cmPanel.cmin_val,
                          vmax=self.cmPanel.cmax_val)
            p.axes.set_autoscale_on(False)
            
            p.canvas.draw()

        elif self.dispm == self.DISP_QUICK:
            res = 100
            kwargs = dict(
                pfig=res,
                iData=self.idata
                )
            pfig = ome_eta.display(**kwargs)

        elif self.dispm == self.DISP_FULL:

            self.oe_pfig()

            pass

        return

    def oe_pfig(self):
        """Make an omega-eta polefigure"""
        # some constants

        deg2rad = numpy.pi /180.0
        radius = numpy.sqrt(2)
        offset = 2.5*radius
        nocolor = 'none'
        pdir_cho = 'X'  # to come from choice interactor
        
        # parent window and data

        exp = wx.GetApp().ws
        p = self.GetParent()
        ome_eta = p.data
        hkldata = ome_eta.getData(self.idata)

        # axes/figure

        p.figure.delaxes(p.axes)
        p.axes = p.figure.gca()
        p.axes.set_autoscale_on(True)
        p.axes.set_aspect('equal')
        p.axes.axis([-2.0, offset+2.0, -1.5, 1.5])

        # outlines for pole figure
        
        C1 = Circle((0,0), radius)
        C2 = Circle((offset,0), radius)
        outline_circles = [C1, C2]
        pc_outl = PatchCollection(outline_circles, facecolors=nocolor)
        p.axes.add_collection(pc_outl)

        # build the rectangles
        
        pf_rects = []

        tTh = exp.activeMaterial.planeData.getTTh()[self.idata]
            
        etas  = ome_eta.etaEdges 
        netas = len(etas) - 1
        deta = abs(etas[1] - etas[0])

        omes  = ome_eta.omeEdges 
        nomes = len(omes) - 1
        dome = abs(omes[1] - omes[0])
        
        if pdir_cho == 'X':
            pdir = numpy.c_[1, 0, 0].T  # X
        elif pdir_cho == 'Y':
            pdir = numpy.c_[0, 1, 0].T  # X
        elif pdir_cho == 'Z':
            pdir = numpy.c_[0, 0, 1].T  # X
            pass

        ii = 0
        for i in range(nomes):
            for j in range(netas):
                qc = makeMSV(tTh, etas[j] + 0.5*deta, omes[i] + 0.5*dome)

                qll = makeMSV(tTh,        etas[j],        omes[i])                
                qlr = makeMSV(tTh, etas[j] + deta,        omes[i])
                qur = makeMSV(tTh, etas[j] + deta, omes[i] + dome)
                qul = makeMSV(tTh,        etas[j], omes[i] + dome)

                pdot_p = numpy.dot(qll.T, pdir) >= 0 \
                    and numpy.dot(qlr.T, pdir) >= 0 \
                    and numpy.dot(qur.T, pdir) >= 0 \
                    and numpy.dot(qul.T, pdir) >= 0

                pdot_m = numpy.dot(qll.T, pdir) < 0 \
                    and numpy.dot(qlr.T, pdir) < 0 \
                    and numpy.dot(qur.T, pdir) < 0 \
                    and numpy.dot(qul.T, pdir) < 0

                if pdot_p:
                    sgn = 1.0
                    ii += 1
                elif pdot_m:
                    sgn = -1.0
                    ii += 1
                elif not pdot_p and not pdot_m:
                    continue

                # the vertex chords
                qll = makeMSV(tTh,        etas[j],        omes[i]) - sgn * pdir
                qlr = makeMSV(tTh, etas[j] + deta,        omes[i]) - sgn * pdir
                qur = makeMSV(tTh, etas[j] + deta, omes[i] + dome) - sgn * pdir
                qul = makeMSV(tTh,        etas[j], omes[i] + dome) - sgn * pdir

                nll = columnNorm(qll)
                nlr = columnNorm(qlr)
                nur = columnNorm(qur)
                nul = columnNorm(qul)

                if pdir_cho == 'X':
                    pqll = nll*unitVector(qll[[1, 2]].reshape(2, 1))
                    pqlr = nlr*unitVector(qlr[[1, 2]].reshape(2, 1))
                    pqur = nur*unitVector(qur[[1, 2]].reshape(2, 1))
                    pqul = nul*unitVector(qul[[1, 2]].reshape(2, 1))
                elif pdir_cho == 'Y':
                    pqll = nll*unitVector(qll[[0, 2]].reshape(2, 1))
                    pqlr = nlr*unitVector(qlr[[0, 2]].reshape(2, 1))
                    pqur = nur*unitVector(qur[[0, 2]].reshape(2, 1))
                    pqul = nul*unitVector(qul[[0, 2]].reshape(2, 1))
                elif pdir_cho == 'Z':
                    pqll = nll*unitVector(qll[[0, 1]].reshape(2, 1))
                    pqlr = nlr*unitVector(qlr[[0, 1]].reshape(2, 1))
                    pqur = nur*unitVector(qur[[0, 1]].reshape(2, 1))
                    pqul = nul*unitVector(qul[[0, 1]].reshape(2, 1))

                xy = numpy.hstack([pqll, pqlr, pqur, pqul]).T
                if sgn == -1:
                    xy[:, 0] = xy[:, 0] + offset
                pf_rects.append(Polygon(xy, aa=False))
                pass
            pass
        
        cmap=matplotlib.cm.jet
        pf_coll = PatchCollection(pf_rects, cmap=cmap, edgecolors='None')
        pf_coll.set_array(numpy.array(hkldata.T.flatten()))
        p.axes.add_collection(pf_coll)

        p.canvas.draw()
        p.axes.axis('off')
        return
    #
    #                     ========== *** Event Callbacks
    #
    def OnHKLChoice(self, e):
        """HKL selection"""
        self.idata = self.hkl_cho.GetSelection()
        self.update()
        
        return
    
    def OnDispChoice(self, e):
        """HKL selection"""
        self.dispm = self.disp_cho.GetSelection()
        self.update()
        
        return
    def OnUpdate(self, e):
        """Update canvas"""
        self.update()
        return

    pass # end class
#
# -----------------------------------------------END CLASS:  sphOpts

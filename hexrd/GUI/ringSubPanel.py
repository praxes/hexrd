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
"""
"""
import wx

from hexrd.valUnits import valWUnit
from hexrd.XRD.Material import Material

from hexrd.GUI.guiConfig    import WindowParameters as WP
from hexrd.GUI.guiUtilities import makeTitleBar
from hexrd.GUI.selectHKLs   import selectHKLsDialog as hklsDlg
#
#  Data
#
from hexrd.XRD.crystallography import dUnit as WAVELENGTH_UNIT
from hexrd.XRD.crystallography import processWavelength

# AngstromTimesKev = 12.39854 # from APS site (lose digits accuracy this way)
AngstromTimesKev = processWavelength(1.0)
#
widChoices = ['two theta', 'strain']
widStrain = 1; widTwoTheta = 0
#
#  Tooltip
#
dfltToolTip = r"""Press this button to make the
current value (wavelength/strain/2-theta) the default
value for new materials.  Applies only for the session.
"""
#
# ---------------------------------------------------CLASS:  ringPanel
#
class ringPanel(wx.Panel):
    """ringPanel """
    def __init__(self, parent, id, **kwargs):
	"""Constructor for ringPanel."""
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

        self.tbarSizer = makeTitleBar(self, 'Rings',
                                      color=WP.TITLEBAR_BG_COLOR_PANEL1)
        #
        #  b.  Wavelength
        #
        self.dfwv_but  = wx.Button(self, wx.NewId(), 'Make Default')
        self.dfwv_but.SetToolTipString(dfltToolTip)

        self.wave_lab = wx.StaticText(self, wx.NewId(), 
                                        'Wavelength:', 
                                        style=wx.ALIGN_RIGHT)
        self.waveAng_txt = wx.TextCtrl(self, wx.NewId(), 
                                       value='0', 
                                       style=wx.RAISED_BORDER | wx.TE_PROCESS_ENTER)
        self.waveAng_lab = wx.StaticText(self, wx.NewId(), 
                                         WAVELENGTH_UNIT, 
                                         style=wx.ALIGN_RIGHT)
        self.waveKEV_txt = wx.TextCtrl(self, wx.NewId(), 
                                       value='0', 
                                       style=wx.RAISED_BORDER | wx.TE_PROCESS_ENTER)
        self.waveKEV_lab = wx.StaticText(self, wx.NewId(), 
                                         'keV', 
                                         style=wx.ALIGN_RIGHT)
        #
        #  c.  Edit HKLs
        #
        self.hkl_but  = wx.Button(self, wx.NewId(), 'Edit HKLs')
        #
        #  d.  Ring widths
        #
        self.dfwd_but  = wx.Button(self, wx.NewId(), 'Make Default')
        self.dfwd_but.SetToolTipString(dfltToolTip)

        self.width_lab = wx.StaticText(self, wx.NewId(), 
                                        'Ring Width:', 
                                        style=wx.ALIGN_RIGHT)
        self.width_cho = wx.Choice(self, wx.NewId(), choices=widChoices)

        self.width_txt = wx.TextCtrl(self, wx.NewId(), 
                                     value='1.0e-3', 
                                     style=wx.RAISED_BORDER | wx.TE_PROCESS_ENTER)
        #
        self.updateFromExp()

        return

    def __makeBindings(self):
        """Bind interactors"""
        self.Bind(wx.EVT_TEXT_ENTER, self.OnAngstromsTxt,    self.waveAng_txt)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnKEVTxt,          self.waveKEV_txt)
        
	self.Bind(wx.EVT_CHOICE,     self.OnWidthChoice,     self.width_cho)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnWidthChoice,     self.width_txt)

        self.Bind(wx.EVT_BUTTON,     self.OnEditHKLs,        self.hkl_but)
        #
        #  Defaults
        #
        self.Bind(wx.EVT_BUTTON, self.OnDefaultWavelength, self.dfwv_but)
        self.Bind(wx.EVT_BUTTON, self.OnDefaultWidth,      self.dfwd_but)

        return

    def __makeSizers(self):
	"""Lay out the interactors"""
        #
        #  Options sizer
        #
        nrow = 3; ncol = 6; padx = pady = 5
	optSizer = wx.FlexGridSizer(nrow, ncol, padx, pady)
	optSizer.AddGrowableCol(0, 1)
        #  * Row 1:  wavelength
        optSizer.Add(self.dfwv_but,    0, wx.ALIGN_LEFT)
        optSizer.Add(self.wave_lab,    1, wx.ALIGN_RIGHT)
        optSizer.Add(self.waveAng_txt, 0, wx.ALIGN_RIGHT|wx.LEFT,   5)
        optSizer.Add(self.waveAng_lab, 0, wx.ALIGN_LEFT |wx.RIGHT,  5)
        optSizer.Add(self.waveKEV_txt, 0, wx.ALIGN_RIGHT|wx.LEFT,   5)
        optSizer.Add(self.waveKEV_lab, 0, wx.ALIGN_LEFT)
        #  * Row 3:  ring widths
        optSizer.Add(self.dfwd_but,    0, wx.ALIGN_LEFT)
        optSizer.Add(self.width_lab,   1, wx.ALIGN_RIGHT)
        optSizer.Add(self.width_cho,   0,  wx.ALIGN_RIGHT|wx.LEFT,  5)
        optSizer.AddSpacer(1)
        optSizer.Add(self.width_txt,   0, wx.ALIGN_RIGHT|wx.LEFT,  5)
        optSizer.AddSpacer(1)
        #  * Row 4:  hkls button
        optSizer.AddSpacer(1)
        optSizer.AddSpacer(1)
        optSizer.Add(self.hkl_but, 0, wx.ALIGN_RIGHT   | wx.TOP, 5)
        optSizer.AddSpacer(1)
        optSizer.AddSpacer(1)
        optSizer.AddSpacer(1)
        #
        #  Main sizer
        #
	self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.sizer.Add(self.tbarSizer, 0, wx.EXPAND|wx.ALIGN_CENTER|wx.BOTTOM, 5)
        self.sizer.Add(optSizer,       0, wx.EXPAND|wx.ALIGN_CENTER)

	return
    
    def __showWavelength(self):
        """Set wavelength values in boxes

        NOTE:  wavelength is returned in angstroms, but set in keV
"""
        exp = wx.GetApp().ws
        ang = exp.activeMaterial.planeData.wavelength
        
        kev = AngstromTimesKev/ang

        str_kev = '%.6g' % kev
        str_ang = '%.6g' % ang
        self.waveKEV_txt.ChangeValue(str_kev)
        self.waveAng_txt.ChangeValue(str_ang)

        return

    def __showRingWidth(self):
        """Show the ring width"""
        app = wx.GetApp()
        exp = app.ws
        mpd = exp.activeMaterial.planeData
        
        if mpd.tThWidth is None:
            self.width_cho.SetSelection(widStrain)
            self.width_txt.ChangeValue(str(mpd.strainMag))
        else:
            self.width_cho.SetSelection(widTwoTheta)
            self.width_txt.ChangeValue(str(mpd.tThWidth))
            pass

        return
    #
    # ============================== API
    #
    #                     ========== *** Access Methods
    #
    def updateFromExp(self):
        """Set interactor values from experiment"""
        self.__showWavelength()
        self.__showRingWidth()
        return
    #
    #                     ========== *** Event Callbacks
    #
    def OnEditHKLs(self, evt):
        """Callback for hkl_but"""

        app = wx.GetApp()
        exp = app.ws
        mat = exp.activeMaterial

        dlg = hklsDlg(self, wx.NewId(), mat)

        if dlg.ShowModal() == wx.ID_OK:
            mat.planeData.exclusions = dlg.getExclusions()
            pass

        app.getCanvas().update()

        evt.Skip()

        return

    def OnKEVTxt(self, evt):
        """Callback for waveAng_txt choice"""
        app = wx.GetApp()
        exp = app.ws
        mat = exp.activeMaterial

        try:
            kev = float(self.waveKEV_txt.GetValue())
            print 'setting kev: ', kev
            mat.planeData.wavelength = kev
        except Exception as e:
            msg = 'Failed to set wavelength:\n%s' % str(e)
            wx.MessageBox(msg)
            pass
        
        self.__showWavelength()

        app.getCanvas().update()

        evt.Skip()

        return

    def OnAngstromsTxt(self, evt):
        """Callback for waveAng_txt choice"""

        app = wx.GetApp()
        exp = app.ws
        mat = exp.activeMaterial

        try:
            ang = float(self.waveAng_txt.GetValue())
            kev = AngstromTimesKev/ang
            mat.planeData.wavelength = float(kev)
        except Exception as e:
            msg = 'Failed to set wavelength:\n%s' % str(e)
            wx.MessageBox(msg)
            pass
            
        self.__showWavelength()

        app.getCanvas().update()

        evt.Skip()

        return
    
    def OnWidthChoice(self, evt):
        """Two theta magnitude set"""
        app = wx.GetApp()
        exp = app.ws
        mpd = exp.activeMaterial.planeData

        cho = self.width_cho.GetSelection()
        if cho == widStrain:
            try:
                val = float(self.width_txt.GetValue())
                mpd.tThWidth  = None
                mpd.strainMag = val
            except:
                wx.MessageBox('failed to set ring widths')
                self.width_txt.ChangeValue(str(mpd.strainMag))
                pass
        else:
            try:
                val = float(self.width_txt.GetValue())
                # mat.strainMag = no change
                mpd.tThWidth  = val
            except:
                wx.MessageBox('failed to set ring widths')
                self.width_txt.ChangeValue(str(mpd.tthWidth))
                pass
            pass

        # update anything?
        app.getCanvas().update()

        evt.Skip()

        return

    def OnDefaultWavelength(self, e):
        """Make the current wavelength the default for the session"""
        kev = float(self.waveKEV_txt.GetValue())
        Material.DFLT_KEV = valWUnit('wavelength', 'energy', kev, 'keV')
        print 'new default for kev', Material.DFLT_KEV
        
        return
    
    def OnDefaultWidth(self, e):
        """Make the current width (stran/2-theta) the default for the session"""
        cho = self.width_cho.GetSelection()
        if cho == widStrain:
            try:
                # set default 2-theta to None and default strain to value
                val = float(self.width_txt.GetValue())
                Material.DFLT_STR = val
                Material.DFLT_TTH = None
            except:
                wx.MessageBox('failed to set default strain width')
                pass
        else:
            try:
                # set default 2-theta to value
                val = float(self.width_txt.GetValue())
                # mat.strainMag = no change
                Material.DFLT_TTH = val
            except:
                wx.MessageBox('failed to set default 2-theta width')
                pass
            pass

        msg  = 'new default:  (tth) %s, (strain) %s'
        vals = (str(Material.DFLT_TTH), str(Material.DFLT_STR))
        print msg % vals
        
        return
    
    pass # end class
#
# -----------------------------------------------END CLASS:  ringPanel

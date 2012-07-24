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
"""Panel for reading frames from GE detector
"""
import os, sys
import wx
import wx.lib.mixins.listctrl  as  listMixins

from hexrd.xrd import detector
from hexrd.xrd.Experiment import *

from hexrd.wx.guiconfig    import WindowParameters as WP
from hexrd.wx.guiUtilities import ResetChoice, makeTitleBar
from hexrd.wx.canvaspanel  import CanvasPanel
#
#  DATA
#
#  * Image Mode Choices
#
IMGMODE_SF  = 'single frame'
IMGMODE_MF  = 'multiframe'
MODE_CHOICES = [IMGMODE_SF, IMGMODE_MF]
IMG_MODES = [ImageModes.SINGLE_FRAME, ImageModes.MULTI_FRAME]
#
IMAGE_MODE_DICT = {
    IMGMODE_SF: IMG_MODES[0],
    IMGMODE_MF: IMG_MODES[1]
    }
IMAGE_MODE_DICT_SEL = dict(zip(IMG_MODES, range(len(MODE_CHOICES))))
#
#  * Dark file choices
#
DARK_CHO_NONE  = 'no dark image'
DARK_CHO_FILE  = 'dark file'
DARK_CHO_EMPTY = 'empty frames'
DARK_CHOICES = [DARK_CHO_NONE, DARK_CHO_FILE, DARK_CHO_EMPTY]
DARK_MODES   = [ReaderInput.DARK_MODE_NONE, ReaderInput.DARK_MODE_FILE, ReaderInput.DARK_MODE_EMPTY]
DARK_MODE_DICT = dict(zip(DARK_CHOICES, DARK_MODES))
DARK_MODE_DICT_INV = dict(zip(DARK_MODES, DARK_CHOICES))
#
#  * Aggregation choices
#
AGG_CHO_NONE = 'SINGLE FRAMES'
AGG_CHO_SUM  = 'average over all frames'
AGG_CHO_MAX  = 'max over all frames'
AGG_CHO_MIN  = 'min over all frames'
AGG_CHOICES  = [AGG_CHO_NONE, AGG_CHO_SUM, AGG_CHO_MAX, AGG_CHO_MIN]
AGG_MODE_DICT     = dict(zip(AGG_CHOICES, ReaderInput.AGG_MODES))
AGG_MODE_DICT_INV = dict(zip(ReaderInput.AGG_MODES, AGG_CHOICES))
#
#  * FLIP choices
#
FLIP_CHO_NONE  = 'no flip'
FLIP_CHO_V     = 'vertical'
FLIP_CHO_H     = 'horizontal'
FLIP_CHO_180   = '180 degrees'
FLIP_CHO_M90   = '-90 degrees'
FLIP_CHO_P90   = '+90 degrees'
FLIP_CHOICES   = [FLIP_CHO_NONE, FLIP_CHO_V, FLIP_CHO_H, FLIP_CHO_180, FLIP_CHO_M90, FLIP_CHO_P90]
FLIP_MODE_DICT = dict(zip(FLIP_CHOICES, ReaderInput.FLIP_MODES))
FLIP_MODE_DICT_INV = dict(zip(ReaderInput.FLIP_MODES, FLIP_CHOICES))
#
#  Utility vFunctions
#
def getValStr(r, i):
    """Return a string of values for display

    r -- a list of tuples
    i -- index

    RETURNS

    s -- a string showing tuple-element i from each member of r
"""
    s = ''
    for ri in r:
        s += ' %s' % str(ri[i])
        pass

    return s

def updVals(r, ind, s):  # Is this used?
    """Return a set of values from a string

    r   -- a list of tuples, to be modified
    s   -- a string showing tuple-element i from each member of r
    ind -- index

    RETURNS

    rn -- the new, updated r
"""
    rn = []
    vals = [float(si) for si in s.split()]
    for i in range(len(r)):
        li = list(r[i]) # from tuple
        li[ind] = vals[i]
        rn.append(tuple(li))
        pass

    return rn
#
# ---------------------------------------------------CLASS:  geReaderPanel
#
class geReaderPanel(wx.Panel):
    """geReaderPanel """
    def __init__(self, parent, id, **kwargs):
	"""Constructor for geReaderPanel

        The approach here is to maintain an active reader
        in the underlying Experiment module.  Each interactor
        then modifies attributes of that reader.
"""
	#
	wx.Panel.__init__(self, parent, id, **kwargs)
	#
        #  Data
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
        self.update()
	#
	return
    #
    # ============================== Internal Methods
    #
    def __makeObjects(self):
        """Add interactors"""
        exp = wx.GetApp().ws

        self.tbarSizer = makeTitleBar(self, 'GE Reader Panel',
                                      color=WP.TITLEBAR_BG_COLOR_PANEL1)
        #
        #  Reader List
        #
        self.curr_lab = wx.StaticText(self, wx.NewId(),
                                        'Current Reader', style=wx.ALIGN_CENTER)
        self.rdrs_cho = wx.Choice(self, wx.NewId(),
                                  choices=[r.name for r in exp.savedReaders])
        #self.save_but = wx.Button(self, wx.NewId(), 'Save Reader')
        self.new_but  = wx.Button(self, wx.NewId(), 'New Reader')
        #
        #  Reader Name
        #
        self.name_lab = wx.StaticText(self, wx.NewId(),
                                        'READER NAME', style=wx.ALIGN_CENTER)
        self.name_txt = wx.TextCtrl(self, wx.NewId(), value=ReaderInput.DFLT_NAME,
                                      style=wx.RAISED_BORDER|wx.TE_PROCESS_ENTER)
        #
        #  Mode interactors
        #
        self.mode_lab = wx.StaticText(self, wx.NewId(), 'Image Mode',
                                      style=wx.ALIGN_RIGHT)
        self.mode_cho = wx.Choice(self, wx.NewId(), choices=MODE_CHOICES)
        #
        #
        #  Image and dark file names
        #
        self.img_but    = wx.Button(self, wx.NewId(), 'Select Image Files')
        self.dir_but    = wx.Button(self, wx.NewId(), 'Change Image Folder')

        self.drk_lab = wx.StaticText(self, wx.NewId(), 'Dark Mode',
                                      style=wx.ALIGN_RIGHT)
        self.drk_cho = wx.Choice(self, wx.NewId(), choices=DARK_CHOICES)
        self.drk_but = wx.Button(self, wx.NewId(), 'Select Dark File')
        #
        #  Aggregation
        #
        self.agg_lab = wx.StaticText(self, wx.NewId(), 'Frame Aggregation',
                                      style=wx.ALIGN_RIGHT)
        self.agg_cho = wx.Choice(self, wx.NewId(), choices=AGG_CHOICES)
        #
        #  Action buttons
        #
        self.files_lab = wx.StaticText(self, wx.NewId(), 'Image Files',
                                      style=wx.ALIGN_RIGHT)
        self.read_lab = wx.StaticText(self, wx.NewId(), 'Read',
                                      style=wx.ALIGN_RIGHT)

        self.read_but   = wx.Button(self, wx.NewId(), 'Load')
        self.browse_lab = wx.StaticText(self, wx.NewId(), 'Browse Frames',
                                      style=wx.ALIGN_RIGHT)
        self.browse_spn = wx.SpinCtrl(self, wx.NewId(), min=0, initial=0)
        self.browse_inf = wx.TextCtrl(self, wx.NewId(), value='',
                                      style=wx.RAISED_BORDER|wx.TE_READONLY)
	self.sizer = wx.BoxSizer(wx.HORIZONTAL|wx.VERTICAL)
        #
        #  Orientation
        #
        self.flip_lab = wx.StaticText(self, wx.NewId(), 'Image Orientation',
                                      style=wx.ALIGN_RIGHT)
        self.flip_cho = wx.Choice(self, wx.NewId(), choices=FLIP_CHOICES)
        #
        #  Subpanels
        #
        self.sp_single = SF_Subpanel(self, wx.NewId())
        self.sp_multi  = MF_Subpanel(self, wx.NewId())
        self.sp_info   = infoPanel(self, wx.NewId())

        return

    def __makeBindings(self):
        """Bind interactors"""
        self.Bind(wx.EVT_CHOICE, self.OnModeChoice, self.mode_cho)

        self.Bind(wx.EVT_TEXT_ENTER, self.OnNameChange, self.name_txt)

        self.Bind(wx.EVT_BUTTON, self.OnDarkBut,    self.drk_but)
        self.Bind(wx.EVT_BUTTON, self.OnImgBut,     self.img_but)
        self.Bind(wx.EVT_BUTTON, self.OnImgDirBut,  self.dir_but)
        self.Bind(wx.EVT_BUTTON, self.OnReadBut,    self.read_but)

        self.Bind(wx.EVT_CHOICE,     self.OnDarkChoice, self.drk_cho)
        self.Bind(wx.EVT_CHOICE,     self.OnAggChoice,  self.agg_cho)
        self.Bind(wx.EVT_CHOICE,     self.OnFlipChoice,  self.flip_cho)

        self.Bind(wx.EVT_SPINCTRL,   self.OnBrowseSpin,  self.browse_spn)

        self.Bind(wx.EVT_CHOICE, self.OnReaderChoice, self.rdrs_cho)
        #self.Bind(wx.EVT_BUTTON, self.OnReaderSave,   self.save_but)
        self.Bind(wx.EVT_BUTTON, self.OnReaderNew,    self.new_but)

        return

    def __makeSizers(self):
	"""Lay out the interactors"""

        nrow = 9; ncol = 4; padx = 5; pady = 5
	self.fgsizer = wx.FlexGridSizer(nrow, ncol, padx, pady) # m x n, paddings
	self.fgsizer.AddGrowableCol(2, 1)
	self.fgsizer.AddGrowableCol(3, 1)
	#self.fgsizer.AddGrowableRow(num, proportion)
        self.fgsizer.Add(self.curr_lab, 0, wx.ALIGN_RIGHT)
        self.fgsizer.Add(wx.Window(self, -1), 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.fgsizer.Add(self.rdrs_cho, 0, wx.ALIGN_RIGHT)
        #self.fgsizer.Add(self.save_but, 0, wx.ALIGN_RIGHT)
        self.fgsizer.Add(self.new_but,  0, wx.ALIGN_RIGHT)

        self.fgsizer.Add(self.name_lab,       0, wx.ALIGN_RIGHT)
        self.fgsizer.Add(wx.Window(self, -1), 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.fgsizer.Add(wx.Window(self, -1), 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.fgsizer.Add(self.name_txt,       0, wx.EXPAND)

        self.fgsizer.Add(self.mode_lab,       0, wx.ALIGN_RIGHT)
        self.fgsizer.Add(wx.Window(self, -1), 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.fgsizer.Add(wx.Window(self, -1), 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.fgsizer.Add(self.mode_cho,       0, wx.ALIGN_RIGHT)

        self.fgsizer.Add(self.agg_lab,    0, wx.ALIGN_RIGHT)
        self.fgsizer.Add(wx.Window(self, -1), 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.fgsizer.Add(wx.Window(self, -1), 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.fgsizer.Add(self.agg_cho,    0, wx.EXPAND|wx.ALIGN_CENTER)

        self.fgsizer.Add(self.flip_lab,    0, wx.ALIGN_RIGHT)
        self.fgsizer.Add(wx.Window(self, -1), 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.fgsizer.Add(wx.Window(self, -1), 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.fgsizer.Add(self.flip_cho,    0, wx.EXPAND|wx.ALIGN_CENTER)

        self.fgsizer.Add(self.drk_lab,    0, wx.ALIGN_RIGHT)
        self.fgsizer.Add(wx.Window(self, -1), 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.fgsizer.Add(self.drk_cho,    0, wx.EXPAND|wx.ALIGN_CENTER)
        self.fgsizer.Add(self.drk_but,    0, wx.ALIGN_RIGHT)

        self.fgsizer.Add(self.files_lab,      0, wx.ALIGN_RIGHT)
        self.fgsizer.AddSpacer(1)
        self.fgsizer.Add(self.dir_but,        0, wx.ALIGN_RIGHT)
        self.fgsizer.Add(self.img_but,        0, wx.ALIGN_RIGHT)

        self.fgsizer.Add(wx.Window(self, -1), 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.fgsizer.Add(self.read_lab,       0, wx.ALIGN_RIGHT)
        self.fgsizer.Add(wx.Window(self, -1), 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.fgsizer.Add(self.read_but,       0, wx.ALIGN_RIGHT)

        self.fgsizer.Add(wx.Window(self, -1), 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.fgsizer.Add(self.browse_lab,       0, wx.ALIGN_RIGHT)
        self.fgsizer.Add(self.browse_spn,       0, wx.ALIGN_RIGHT)
        self.fgsizer.Add(self.browse_inf,       0, wx.ALIGN_CENTER | wx.EXPAND)

	self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.sizer.Add(self.tbarSizer, 0,
                       wx.EXPAND|wx.ALIGN_CENTER|wx.BOTTOM, 10)
	self.sizer.Add(self.fgsizer,    0, wx.EXPAND|wx.ALIGN_RIGHT)

	self.sizer.Add(self.sp_single,   0, wx.EXPAND|wx.ALIGN_RIGHT)
	self.sizer.Add(self.sp_multi,    0, wx.EXPAND|wx.ALIGN_RIGHT)

	self.sizer.Add(self.sp_info,    0, wx.EXPAND|wx.ALIGN_RIGHT)

	return
    #
    # ============================== API
    #
    #                     ========== *** Access Methods
    #
    def update(self):
        """Update panel and interactors from the current exp"""
        exp = wx.GetApp().ws
        rdr = exp.activeReader

        # Reader line:  reset choice interactor
        rnames = [r.name for r in exp.savedReaders]
        ResetChoice(self.rdrs_cho, rnames, rdr.name)

        # Name line
        self.name_txt.ChangeValue(rdr.name)
        #
        mode = rdr.imageMode
        self.mode_cho.SetSelection(IMAGE_MODE_DICT_SEL[mode])
        # Agg choice
        self.agg_cho.SetStringSelection(AGG_MODE_DICT_INV[rdr.aggMode])
        # Image Orientation
        self.flip_cho.SetStringSelection(FLIP_MODE_DICT_INV[rdr.flipMode])
        # Dark mode
        self.drk_cho.SetStringSelection(DARK_MODE_DICT_INV[rdr.darkMode])

        # Mode Subpanel
        self.sizer.Show(self.sp_single, (mode == ImageModes.SINGLE_FRAME))
        self.sizer.Show(self.sp_multi,  (mode == ImageModes.MULTI_FRAME))
        # Image names
        if mode == ImageModes.SINGLE_FRAME:
            lctrl = self.sp_single.file_lctrl
        else:
            lctrl = self.sp_multi.file_lctrl
            pass

        lctrl.DeleteAllItems()
        for n in rdr.imageNames:
            index = lctrl.InsertStringItem(sys.maxint, n)
            if mode == ImageModes.MULTI_FRAME:
                vals = rdr.imageNameD[n]
                lctrl.SetStringItem(index, 1, str(vals[0]))
                lctrl.SetStringItem(index, 3, str(vals[1]))
                lctrl.SetStringItem(index, 4, str(vals[2]))
                lctrl.SetStringItem(index, 5, str(vals[3]))
                # add total number of frames available
                try:
                    d = rdr.imageDir
                    r = detector.ReadGE((os.path.join(d, n), 0))
                    nframe = r.getNFrames()
                    lctrl.SetStringItem(index, 2, str(nframe))
                except:
                    lctrl.SetStringItem(index, 2, '(error)')
                    pass
                pass
            pass

        # info panel
        self.sp_info.update()

        self.sizer.Layout()

        return

    updateFromExp = update
    #
    #                     ========== *** Event Callbacks
    #
    def NotYetImpl(self, e):
        """Not implemented"""

        print dir(e)
        wx.MessageBox('not yet implemented')

        return

    def OnImgDirBut(self, e):
        """Change image directory"""
        #
        #  Change the image directory
        #
        dlg = wx.DirDialog(self, 'Change Image Directory',
                           style=wx.DD_DIR_MUST_EXIST)
        if dlg.ShowModal() == wx.ID_OK:
            dir = str(dlg.GetPath())
            if (dir):
                #
                #  Update image directory in Experiment
                #
                exp = wx.GetApp().ws
                exp.activeReader.imageDir = dir
                self.sp_info.update()
                pass
            pass
        dlg.Destroy()

        return

    def OnReaderChoice(self, e):
        """Change the active reader"""
        exp = wx.GetApp().ws

        sel = self.rdrs_cho.GetSelection()
        if sel >= 0:
            exp.activeReader = sel
            self.update()

            self.browse_inf.Enable(False)
            self.browse_spn.Enable(False)

            pass

        return

    def OnReaderNew(self, e):
        """A new reader has been requested"""
        exp = wx.GetApp().ws
        exp.newReader()
        self.update()

        self.browse_inf.Enable(False)
        self.browse_spn.Enable(False)

        return

    def OnNameChange(self, e):
        """Name entry has changed"""
        exp = wx.GetApp().ws
        n = str(self.name_txt.GetValue())
        names = [r.name for r in exp.savedReaders]
        if names.count(n) > 0:
            wx.MessageBox('name already in use')
        else:
            exp.activeReader.name = n
            pass

        self.update()

        return

    def OnBrowseSpin(self, e):
        """Browse spinner has been pressed"""
        app = wx.GetApp()
        exp = app.ws
        exp.readImage(self.browse_spn.GetValue())
        app.getCanvas().update(newImage=True)

        return

    def OnDarkChoice(self, e):
        """Dark mode choice has been made"""
        val = e.GetString()
        mode = DARK_MODE_DICT[val]
        exp = wx.GetApp().ws
        exp.activeReader.darkMode = mode
        #
        #  Enable/disable other interactors
        #
        self.drk_but.Enable(mode == ReaderInput.DARK_MODE_FILE)

        #  Update info window

        self.sp_info.update()

        return

    def OnAggChoice(self, e):
        """Aggregation function selection"""
        val = e.GetString()
        mode =  AGG_MODE_DICT[val]
        exp = wx.GetApp().ws
        exp.activeReader.aggFun = mode
        #
        #  Enable/disable other interactors
        #
        self.browse_spn.Enable(mode == ReaderInput.AGG_FUN_NONE)

        #  Update info window

        self.sp_info.update()

        return

    def OnModeChoice(self, e):
        """Detector mode selected"""
        #
        mode = e.GetString()
        # Set mode in active reader
        exp = wx.GetApp().ws
        exp.activeReader.imageMode = IMAGE_MODE_DICT[mode]
        # Show panels according to mode
        self.update()

        return

    def OnFlipChoice(self, e):
        """Flip mode chosen"""
        print 'flip mode:  ', e.GetString()
        wx.GetApp().ws.activeReader.flipMode = FLIP_MODE_DICT[e.GetString()]
        return

    def OnDarkBut(self, e):
        """Load dark file names with file dialogue"""
        #
        #  !! Check that "subtract dark" is true
        #
        dlg = wx.FileDialog(self, 'Select Dark Image')
        if dlg.ShowModal() == wx.ID_OK:
            dir = str(dlg.GetDirectory())
            fil = str(dlg.GetFilename())
            if (fil):
                #
                #  Set dark file and display name in info box.
                #
                exp = wx.GetApp().ws
                exp.activeReader.darkDir = dir
                exp.activeReader.darkName = fil
                self.sp_info.update()
                pass
            pass
        dlg.Destroy()

        return

    def OnImgBut(self, e):
        """Load image file names with file dialogue

        NOTE:  converts filenames to str from unicode
"""
        dlg = wx.FileDialog(self, 'Select Images',
                            style=wx.FD_MULTIPLE)
        if dlg.ShowModal() == wx.ID_OK:
            d = str(dlg.GetDirectory())
            fnames = [str(p) for p in dlg.GetFilenames()]
            if (fnames):
                #
                #  Set image file list and display name in box.
                #
                fnames.sort()
                exp = wx.GetApp().ws
                print d, fnames
                exp.activeReader.imageDir = d
                exp.activeReader.imageNames = fnames

                pass

            self.update()
            pass
        dlg.Destroy()

        return

    def OnDrkSubtract(self, e):
        """Subtract dark checkbox

        * No other effects until read image button is pressed
"""
        wx.GetApp().ws.drkSubtract = self.drk_box.GetValue()

        return

    def OnReadBut(self, e):
        """Read the frames"""

        """Create a GE Reader"""

        mainFrame = wx.GetApp().GetTopWindow()
        mainFrame.SetStatusText('reading images ...')

        app = wx.GetApp()
        exp = app.ws

        try:
            exp.readImage()
        except Exception as ex:
            msg = 'Failed to read image:\n%s' % str(ex)

            wx.MessageBox(msg)
            #raise
            pass

        mainFrame.SetStatusText('done; updating canvas now ...')

        # update canvas

        app.getCanvas().update(newImage=True)

        # update browse spinner if more frames are available

        self.browse_inf.Enable()
        self.browse_spn.Enable()

        self.browse_inf.SetValue('total frames: %d' % exp.numFramesTotal)
        self.browse_spn.SetValue(exp.curFrameNumber)
        self.browse_spn.SetRange(1, exp.numFramesTotal)

        mainFrame.SetStatusText('done')

        return

    pass # end class
#
# -----------------------------------------------END CLASS:  geReaderPanel
# ---------------------------------------------------CLASS:  MF_Subpanel
#
class MF_Subpanel(wx.Panel):
    """MF_Subpanel """
    def __init__(self, parent, id, **kwargs):
	"""Constructor for MF_Subpanel."""
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

        self.tbarSizer = makeTitleBar(self, ' Multiframe Options ',
                                      color=WP.TITLEBAR_BG_COLOR_PANEL1)

        #
        self.file_lctrl =  self.__makeListCtrl()


        return

    def __makeListCtrl(self):
	"""Make the list control"""
        #
	LStyle = wx.LC_REPORT|wx.LC_SINGLE_SEL
	#
	listctrl = myListCtrl(self, wx.NewId(), style=LStyle)
	listctrl.InsertColumn(0, 'Image File')
	listctrl.InsertColumn(1, 'Empty Frames')
	listctrl.InsertColumn(2, 'Total Frames')
	listctrl.InsertColumn(3, 'Omega Min')
	listctrl.InsertColumn(4, 'Omega Max')
	listctrl.InsertColumn(5, 'Delta')
        listctrl.SetColumnWidth(1, wx.LIST_AUTOSIZE_USEHEADER)
        listctrl.SetColumnWidth(2, wx.LIST_AUTOSIZE_USEHEADER)
        listctrl.SetColumnWidth(3, wx.LIST_AUTOSIZE_USEHEADER)
        listctrl.SetColumnWidth(4, wx.LIST_AUTOSIZE_USEHEADER)
        listctrl.SetColumnWidth(5, wx.LIST_AUTOSIZE_USEHEADER)
        listctrl.SetColumnWidth(0, 200)

	return listctrl

    def __makeBindings(self):
        """Bind interactors"""
        #self.Bind(wx.EVT_LIST_ITEM_DESELECTED, self.OnEditListCtrl, self.file_lctrl)
        # event posted by textedit mixin
        #self.Bind(wx.wxEVT_COMMAND_LIST_END_LABEL_EDIT, self.OnEditListCtrl, self.file_lctrl)

        return

    def __makeSizers(self):
	"""Lay out the interactors"""

	self.sizer = wx.BoxSizer(wx.VERTICAL)

	self.sizer.Add(self.tbarSizer,  0, wx.EXPAND|wx.ALIGN_CENTER)
        self.sizer.Add(self.file_lctrl, 1, wx.EXPAND | wx.ALIGN_CENTER | wx.TOP, 5)

	return
    #
    # ============================== API
    #
    #                     ========== *** Access Methods
    #

    #
    #                     ========== *** Event Callbacks
    #
    def OnEditListCtrl(self, e):
        """List control item Activated

        Update information for that item.
"""
        print 'item selected'
        print e.GetIndex(), e.GetColumn(), e.Column
        #print 'text: ', self.file_lctrl.GetItem(e.GetIndex(), e.GetColumn()).GetText()
        #print 'text: ', self.file_lctrl.GetItem(e.GetIndex(),1).GetText()
        i = e.GetIndex()
        fname = self.file_lctrl.GetItem(i, 0).GetText()
        numEm = self.file_lctrl.GetItem(i, 1).GetText()
        # field 2 is informational/derived

        exp = wx.GetApp().ws
        exp.activeReader.imageNameD[fname] = (int(numEm), omin, omax, odel)

        return

    pass # end class
#
# -----------------------------------------------END CLASS:  MF_Subpanel
# ---------------------------------------------------CLASS:  SF_Subpanel
#
class SF_Subpanel(wx.Panel):
    """SF_Subpanel """
    def __init__(self, parent, id, **kwargs):
	"""Constructor for SF_Subpanel."""
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

        self.tbarSizer = makeTitleBar(self, ' Single Frame Options ',
                                      color=WP.TITLEBAR_BG_COLOR_PANEL1)


        self.file_lctrl = self.__makeListCtrl()

        return

    def __makeListCtrl(self):
	"""Make the list control"""
        #
	LStyle = wx.LC_REPORT|wx.LC_SINGLE_SEL
	#
	listctrl = wx.ListCtrl(self, wx.NewId(), style=LStyle)
	listctrl.InsertColumn(0, 'Image File')
	listctrl.InsertColumn(1, 'Omega')
        listctrl.SetColumnWidth(1, wx.LIST_AUTOSIZE_USEHEADER)
        listctrl.SetColumnWidth(0, 200)

	return listctrl

    def __makeBindings(self):
        """Bind interactors"""
        return

    def __makeSizers(self):
	"""Lay out the interactors"""

	self.sizer = wx.BoxSizer(wx.VERTICAL)

	self.sizer.Add(self.tbarSizer, 0, wx.EXPAND|wx.ALIGN_CENTER)
	self.sizer.Add(self.file_lctrl, 0, wx.EXPAND|wx.ALIGN_CENTER)

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
# -----------------------------------------------END CLASS:  SF_Subpanel
# ---------------------------------------------------CLASS:  infoPanel
#
class infoPanel(wx.Panel):
    """infoPanel """
    def __init__(self, parent, id, **kwargs):
	"""Constructor for infoPanel."""
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

        self.tbarSizer = makeTitleBar(self, ' Information ',
                                      color=WP.TITLEBAR_BG_COLOR_PANEL1)

        #
        #  File lists for display.
        #

        self.drk_txt_lab = wx.StaticText(self, wx.NewId(), 'Dark File',
                                         style=wx.ALIGN_CENTER)
        self.drk_txt = wx.TextCtrl(self, wx.NewId(), value='<no dark file>',
                                   style=wx.RAISED_BORDER)

        self.img_txt_lab = wx.StaticText(self, wx.NewId(), 'Image Directory',
                                         style=wx.ALIGN_CENTER)
        self.img_txt = wx.TextCtrl(self, wx.NewId(), value='<no images loaded>',
                                   style=wx.RAISED_BORDER)

        return

    def __makeBindings(self):
        """Bind interactors"""
        return

    def __makeSizers(self):
	"""Lay out the interactors"""
        nrow = 1; ncol = 2; padx = 5; pady = 5
	self.fgsizer = wx.FlexGridSizer(nrow, ncol, padx, pady) # m x n, paddings
	self.fgsizer.AddGrowableCol(1, 1)
	#self.fgsizer.AddGrowableRow(num, proportion)

	self.fgsizer.Add(self.img_txt_lab, 0, wx.ALIGN_RIGHT|wx.RIGHT, 5)
	self.fgsizer.Add(self.img_txt,     1, wx.EXPAND|wx.ALIGN_CENTER)

	self.fgsizer.Add(self.drk_txt_lab, 0, wx.ALIGN_RIGHT|wx.RIGHT, 5)
	self.fgsizer.Add(self.drk_txt,     1, wx.EXPAND|wx.ALIGN_CENTER)
        #
        #  Main sizer
        #
	self.sizer = wx.BoxSizer(wx.VERTICAL)

	self.sizer.Add(self.tbarSizer, 0, wx.EXPAND|wx.ALIGN_CENTER)
	self.sizer.Add(self.fgsizer,   0, wx.EXPAND|wx.ALIGN_CENTER)

	return
    #
    # ============================== API
    #
    #                     ========== *** Access Methods
    #
    def update(self):
        """Update information"""
        exp = wx.GetApp().ws
        self.drk_txt.SetValue(exp.activeReader.darkFile)
        self.img_txt.SetValue(exp.activeReader.imageDir)

        return
    #
    #                     ========== *** Event Callbacks
    #

    pass # end class
#
# -----------------------------------------------END CLASS:  infoPanel
# ---------------------------------------------------CLASS:  ListCtrl
#
class myListCtrl(wx.ListCtrl,
                 listMixins.ListCtrlAutoWidthMixin,
                 listMixins.TextEditMixin):
    """myListCtrl:  subclassing to include mixins"""
    def __init__(self, parent, ID, **kwargs):
	"""Constructor for ListCtrl"""
	#
	wx.ListCtrl.__init__(self, parent, ID, **kwargs)
        listMixins.ListCtrlAutoWidthMixin.__init__(self)
        listMixins.TextEditMixin.__init__(self)
	#

	return

    def CloseEditor(self, e=None):
        exp = wx.GetApp().ws

        listMixins.TextEditMixin.CloseEditor(self, e)
        i = self.curRow; j = self.curCol

        item = self.GetItem(i, j)

        print 'item (%d, %d) modified' % (i, j)
        print 'item text: ', item.GetText()

        fname = self.GetItem(i, 0).GetText()

        try:
            numEm = int(self.GetItem(i, 1).GetText())
        except:
            numEm = exp.activeReader.imageNameD[fname][0]
            self.SetStringItem(i, j, str(numEm))
            pass

        omin  = self.GetItem(i, 3).GetText()
        omax  = self.GetItem(i, 4).GetText()
        odel  = self.GetItem(i, 5).GetText()

        print 'omega info:  ', omin, omax, odel


        exp.activeReader.imageNameD[fname] = (numEm, omin, omax, odel)
        return


    pass # end class
#
# -----------------------------------------------END CLASS:  ListCtrl

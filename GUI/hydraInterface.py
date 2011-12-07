#! /usr/bin/env python
#
"""User controls for the hydra interface
"""
import wx

from guiConfig import WindowParameters as WP

from hydraCanvas import hydraCanvasFrame
#
# ---------------------------------------------------CLASS:  HydraControlPanel
#
class HydraControlPanel(wx.Panel):
    """HydraControlPanel """
    def __init__(self, parent, id, **kwargs):
	"""Constructor for HydraControlPanel."""
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

        self.__makeTitleBar('Hydra Control Panel')
        #
        #  Reader selectors
        #
        exp = wx.GetApp().ws
        #
        self.img1_lab = wx.StaticText(self, wx.NewId(), 'Image 1', 
                                      style=wx.ALIGN_CENTER)
        self.img1_cho = wx.Choice(self, wx.NewId(), choices=exp.readerNames)
        #
        self.img2_lab = wx.StaticText(self, wx.NewId(), 'Image 2', 
                                      style=wx.ALIGN_CENTER)
        self.img2_cho = wx.Choice(self, wx.NewId(), choices=exp.readerNames)
        #
        self.img3_lab = wx.StaticText(self, wx.NewId(), 'Image 3', 
                                      style=wx.ALIGN_CENTER)
        self.img3_cho = wx.Choice(self, wx.NewId(), choices=exp.readerNames)
        #
        self.img4_lab = wx.StaticText(self, wx.NewId(), 'Image 4', 
                                      style=wx.ALIGN_CENTER)
        self.img4_cho = wx.Choice(self, wx.NewId(), choices=exp.readerNames)
        #
        #  Load Button
        #
        self.load_but  = wx.Button(self, wx.NewId(), 'Load')

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
        self.Bind(wx.EVT_CHOICE, self.OnChooseImg1, self.img1_cho)
        self.Bind(wx.EVT_CHOICE, self.OnChooseImg2, self.img2_cho)
        self.Bind(wx.EVT_CHOICE, self.OnChooseImg3, self.img3_cho)
        self.Bind(wx.EVT_CHOICE, self.OnChooseImg4, self.img4_cho)
        #
        self.Bind(wx.EVT_BUTTON, self.OnLoad, self.load_but)

        return

    def __makeSizers(self):
	"""Lay out the interactors"""
        nrow = 4; ncol = 2; padx = pady = 5
	self.fgsizer = wx.FlexGridSizer(nrow, ncol, padx, pady)
        self.fgsizer.Add(self.img1_lab, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.fgsizer.Add(self.img1_cho, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.fgsizer.Add(self.img2_lab, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.fgsizer.Add(self.img2_cho, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.fgsizer.Add(self.img3_lab, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.fgsizer.Add(self.img3_cho, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.fgsizer.Add(self.img4_lab, 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.fgsizer.Add(self.img4_cho, 0, wx.EXPAND|wx.ALIGN_CENTER)
        #
        self.fgsizer.Add(wx.Window(self, -1), 0, wx.EXPAND|wx.ALIGN_CENTER)
        self.fgsizer.Add(self.load_but, 0, wx.EXPAND|wx.ALIGN_CENTER)
        
	#self.fgsizer.AddGrowableCol(num, proportion)
	#self.fgsizer.AddGrowableRow(num, proportion)
	
	self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.sizer.Add(self.titlebar, 0, wx.EXPAND|wx.ALIGN_CENTER)
	self.sizer.Add(self.fgsizer,  0, wx.EXPAND|wx.ALIGN_CENTER)

	return
    #
    # ============================== API
    #
    #                     ========== *** Access Methods
    #

    #
    #                     ========== *** Event Callbacks
    #
    def OnLoad(self, e):
        """Load button pressed"""
        h = wx.GetApp().ws.hydra
        h.loadImages()
        self.hydraCanvas = hydraCanvasFrame(self, wx.NewId())
        self.hydraCanvas.hcanvas.loadImages()

        return

    def OnChooseImg1(self, e):
        """Image 1 reader chosen"""
        exp = wx.GetApp().ws
        h = exp.hydra
        rname = self.img1_cho.GetStringSelection()
        h.readers[0] = exp.getSavedReader(rname)

        return

    def OnChooseImg2(self, e):
        """Image 2 reader chosen"""
        exp = wx.GetApp().ws
        h = exp.hydra
        rname = self.img2_cho.GetStringSelection()
        h.readers[1] = exp.getSavedReader(rname)

        return

    def OnChooseImg3(self, e):
        """Image 3 reader chosen"""
        exp = wx.GetApp().ws
        h = exp.hydra
        rname = self.img3_cho.GetStringSelection()
        h.readers[2] = exp.getSavedReader(rname)

        return

    def OnChooseImg4(self, e):
        """Image 4 reader chosen"""
        exp = wx.GetApp().ws
        h = exp.hydra
        rname = self.img4_cho.GetStringSelection()
        h.readers[3] = exp.getSavedReader(rname)

        return

    pass # end class
#
# -----------------------------------------------END CLASS:  HydraControlPanel
# ---------------------------------------------------CLASS:  HydraControlFrame
#
class HydraControlFrame(wx.Frame):
    #
    def __init__(self, parent, id, title='Hydra'):
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

        self.hydraPanel = HydraControlPanel(self, wx.NewId())
	#
        # A Statusbar in the bottom of the window
        #
        self.CreateStatusBar()
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
FRAME FOR: hydra user controls
"""
        self.titlebar.SetToolTipString(myToolTip)

        return

    def __makeBindings(self):
        """Bind interactors"""
        return

    def __makeSizers(self):
	"""Lay out the interactors"""
	
	self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.sizer.Add(self.hydraPanel, 0, wx.EXPAND|wx.ALIGN_CENTER)

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
# -----------------------------------------------END CLASS:  HydraControlFrame
#

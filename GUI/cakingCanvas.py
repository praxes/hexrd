#! /usr/bin/env python
#
"""Display for caking output.
"""
import wx

from XRD.Experiment  import PolarRebinOpts as prOpts

from guiConfig    import WindowParameters as WP
from guiUtilities import makeTitleBar

from canvasUtilities import *
#
# ---------------------------------------------------CLASS:  cakeCanvas
#
class cakeCanvas(wx.Panel):
    """cakeCanvas """
    def __init__(self, parent, id, cakeType, data, **kwargs):
	"""Constructor for cakeCanvas panel

        INPUTS
        cakeType -- full caking, ring-based or spherical (omega-eta map)
        data     -- caking output, depending on cakeType
                 .. for rings case, a MultiRingBinned instance
        
        OUTPUTS
        
        DESCRIPTION
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
        #  see if we need cmPanel
        #
        #self.cmPanel = cmapPanel(self, wx.NewId())
        #self.cmPanel.Disable()
        #
        if   self.cakeType == prOpts.CAKE_IMG:
            # full image caking panel
            pass
        elif self.cakeType == prOpts.CAKE_RNG:
            # multiring panel
            self.opt_pan = rngOpts(self, wx.NewId())
            pass
        elif self.cakeType == prOpts.CAKE_SPH:
            # omega-eta panel
            pass
        
        self.__makeFigureCanvas()
        #

        return

    def __makeFigureCanvas(self):
        """Build figure canvas"""
        self.figure = Figure()
        self.canvas = FigureCanvas(self, wx.NewId(), self.figure)

        self.axes   = self.figure.gca()
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
        self.sizer.Add(self. unit_cho, wx.EXPAND|wx.ALIGN_CENTER)
        self.sizer.Add(self.exp_but,  0, wx.EXPAND|wx.ALIGN_CENTER|wx.TOP, 5)
        
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
# -----------------------------------------------END CLASS:  rngOpts

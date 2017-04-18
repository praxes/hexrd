"""panel for reader input"""
import wx

from guiconfig    import WindowParameters as WP
from guiutil import makeTitleBar
#
# ---------------------------------------------------CLASS:  ReaderInfoPanel
#
class ReaderInfoPanel(wx.Panel):
    """ReaderInfoPanel """
    def __init__(self, parent, id, **kwargs):

	wx.Panel.__init__(self, parent, id, **kwargs)
	#
        #  Data
        #
        self.image_dir = ''
        self.image_fname = ''
        #
	#  Window Objects.
	#
        self.__make_objects()
	#
	#  Bindings.
	#
	self.__make_bindings()
	#
	#  Sizing.
	#
	self.__make_sizers()
	#
	self.SetAutoLayout(True)
        self.SetSizerAndFit(self.sizer)
	#
	return
    #
    # ============================== Internal Methods
    #
    def __make_objects(self):
        """Add interactors"""

        self.tbarSizer = makeTitleBar(self, 'Reader Info')
        self.file_but = wx.Button(self, wx.NewId(),
                                  'File'
                                  )
        self.file_txt = wx.TextCtrl(self, wx.NewId(),
                                    value="<none selected>",
                                    style=wx.RAISED_BORDER|wx.TE_READONLY
                                    )
        self.format_lab = wx.StaticText(self, wx.NewId(),
                                        'Format', style=wx.ALIGN_RIGHT
                                        )
        self.format_cho = wx.Choice(self, wx.NewId(),
                                    choices=['hdf5', 'frame-cache']
                                    )
        self.pixel_lab = wx.StaticText(self, wx.NewId(),
                                       'Pixel Pitch', style=wx.ALIGN_RIGHT
                                       )
        self.pixel_txt = wx.TextCtrl(self, wx.NewId(),
                                     value='0.2',
                                     style=wx.RAISED_BORDER
                                     )
        self.option_lab = wx.StaticText(self, wx.NewId(),
                                        'Option', style=wx.ALIGN_RIGHT
                                        )
        self.value_lab = wx.StaticText(self, wx.NewId(),
                                        'Value', style=wx.ALIGN_LEFT
                                        )
        self.option_cho = wx.Choice(self, wx.NewId(),
                                    choices=['path', 'pixel pitch']
                                    )
        self.value_txt = wx.TextCtrl(self, wx.NewId(),
                                     value="/imageseries",
                                     style=wx.RAISED_BORDER
                                     )

    def __make_bindings(self):
        """Bind interactors"""
        self.Bind(wx.EVT_BUTTON, self.OnFileBut, self.file_but)
        
    def __make_sizers(self):
	"""Lay out the interactors"""

	self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.sizer.Add(self.tbarSizer, 0, wx.EXPAND|wx.ALIGN_CENTER)

        nrow = 5; ncol = 2; padx = 5; pady = 5
        self.info_sz = wx.FlexGridSizer(nrow, ncol, padx, pady)
        self.info_sz.AddGrowableCol(0, 0)
        self.info_sz.AddGrowableCol(1, 1)
        self.info_sz.Add(self.file_but,   0, wx.ALIGN_RIGHT)
        self.info_sz.Add(self.file_txt,   0, wx.ALIGN_LEFT|wx.EXPAND)
        self.info_sz.Add(self.format_lab, 0, wx.ALIGN_RIGHT)
        self.info_sz.Add(self.format_cho, 0, wx.ALIGN_LEFT|wx.EXPAND)
        self.info_sz.Add(self.pixel_lab, 0, wx.ALIGN_RIGHT)
        self.info_sz.Add(self.pixel_txt, 0, wx.ALIGN_LEFT|wx.EXPAND)
        self.info_sz.Add(self.option_lab, 0, wx.ALIGN_RIGHT)
        self.info_sz.Add(self.value_lab,  0, wx.ALIGN_LEFT)
        self.info_sz.Add(self.option_cho, 0, wx.ALIGN_RIGHT)
        self.info_sz.Add(self.value_txt,  0, wx.ALIGN_LEFT|wx.EXPAND)


        self.sizer.Add(self.info_sz, 1, wx.EXPAND)
    #
    # ============================== API
    #
    #                     ========== *** Access Methods
    #

    #
    #                     ========== *** Event Callbacks
    #
    def OnFileBut(self, e):
        """Load image file name with file dialogue

        NOTE:  converts filenames to str from unicode
        """
        dlg = wx.FileDialog(self, 'Select Imageseries File',
                            style=wx.FD_FILE_MUST_EXIST)
        if dlg.ShowModal() == wx.ID_OK:
            self.image_dir = str(dlg.GetDirectory())
            self.image_fname = dlg.GetFilename()
            self.file_txt.SetValue(self.image_fname)
        dlg.Destroy()

    pass # end class
# -----------------------------------------------END CLASS:  ReaderInfoPanel
# ---------------------------------------------------CLASS:  ReaderInfoDlg
#
class ReaderInfoDialog(wx.Dialog):
    """Pop-Up for reader file, format and options """
    def __init__(self, parent, id, **kwargs):
	"""Constructor"""
	#
        myStyle = wx.RESIZE_BORDER|wx.DEFAULT_DIALOG_STYLE
	wx.Dialog.__init__(self, parent, id, style=myStyle)
	#
	#  Data Objects.
	#

	#
	#  Windows.
	#
        self.tbarSizer = makeTitleBar(self, 'Reader Info',
                                      color=WP.BG_COLOR_TITLEBAR_FRAME)
        self.infopanel = ReaderInfoPanel(self, -1)
	#
	#  Bindings.
	#
	self._makeBindings()
	#
	#  Sizing.
	#
	self._makeSizers()
	#
	self.SetAutoLayout(True)
        self.SetSizerAndFit(self.sizer)
	#
	#  Events.
	#

	#
	return
    #
    # ============================== Internal Methods
    #
    def _makeBindings(self):
	"""Bind interactors to functions"""

    def _makeSizers(self):
	"""Lay out windows"""
	self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.sizer.Add(self.tbarSizer,  0, wx.EXPAND|wx.ALIGN_CENTER)
        self.sizer.Add(self.CreateSeparatedButtonSizer(wx.OK | wx.CANCEL))
	self.sizer.Add(self.infopanel, 1, wx.EXPAND|wx.ALIGN_CENTER)
    #
    # ============================== API
    #
    #                     ========== *** Access Methods
    #
    def GetInfo(self):
        p = self.infopanel
        d = dict(
            directory = p.image_dir,
            file = p.image_fname,
            format = p.format_cho.GetStringSelection(),
            path = p.value_txt.GetValue(),
            pixel_size = p.pixel_txt.GetValue())
        return d

    pass # end class
#
# -----------------------------------------------END CLASS:  ReaderInfoDlg
#

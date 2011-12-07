#! /usr/bin/env python
#
#  Copyright-Info-Goes-Here
#
"""Log management and display
"""
import thread

import wx
from guiUtilities import makeTitleBar
#
# ---------------------------------------------------CLASS:  logWindow
#
class logWindow(wx.Dialog):
    #
    def __init__(self, parent, id, actionD, title='Log'):
        #
        #  Open slightly offset from parent.
        #
        newPos = parent.GetPosition()
        print 'new position:  ', newPos
        newPos[0] += 50; newPos[1] += 50;
        #
	wx.Dialog.__init__(self, parent, id, title,
                           pos=newPos,
                           style=wx.DEFAULT_DIALOG_STYLE|wx.RESIZE_BORDER)
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
        # self.Show(True)

        self.__run(actionD)
        
	return

    pass # class
    #
    # ============================== Internal Methods
    #
    def __makeObjects(self):
        """Add interactors"""

        self.tbarSizer = makeTitleBar(self, 'Log')
        #
        self.log_pan = logPanel(self, wx.NewId())
	#
        return

    def __makeBindings(self):
        """Bind interactors"""
        return

    def __makeSizers(self):
	"""Lay out the interactors"""
	
	self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.sizer.Add(self.tbarSizer, 0, wx.EXPAND|wx.ALIGN_CENTER)
	self.sizer.Add(self.log_pan,   1, wx.EXPAND|wx.ALIGN_RIGHT)

	return

    def __run(self, actionD):
        """run the action in a thread"""
        action = actionD['exec']
        args   = actionD['args']
        kw     = actionD['kwargs']
        kw.update({'log': self})

        thread.start_new_thread(action, args, kw)
        
        return
    #
    # ============================== API
    #
    #                     ========== *** Access Methods
    #
    def write(self, msg, **kwargs):
        """Write a message to this window"""

        self.log_pan.write(msg, **kwargs)

        return
    #
    #                     ========== *** Event Callbacks
    #
    pass
#
# -----------------------------------------------END CLASS:  logWindow
# ---------------------------------------------------CLASS:  logWindow
#
class logPanel(wx.Panel):
    #
    def __init__(self, parent, id):
        #
	wx.Panel.__init__(self, parent, id)
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

        self.log_txt = wx.TextCtrl(self, wx.NewId(), value='', size=(500,700),
                                   style=wx.RAISED_BORDER|wx.TE_MULTILINE|wx.TE_READONLY)
	#
        return

    def __makeBindings(self):
        """Bind interactors"""
        return

    def __makeSizers(self):
	"""Lay out the interactors"""
	
	self.sizer = wx.BoxSizer(wx.VERTICAL)
	self.sizer.Add(self.log_txt,   1, wx.EXPAND|wx.ALIGN_RIGHT)

	return
    #
    # ============================== API
    #
    #                     ========== *** Access Methods
    #
    def write(self, msg, done=False):
        """Write a message to this window"""
        wx.CallAfter(self.log_txt.AppendText, msg)
        # self.log_txt.AppendText(msg)
        self.Refresh()

        if done:
            doneMsg = '\nCLOSE THIS WINDOW TO CONTINUE\n'
            wx.CallAfter(self.log_txt.AppendText, doneMsg)
            pass
        
        return
    #
    #                     ========== *** Event Callbacks
    #
    pass
#
# -----------------------------------------------END CLASS:  logWindow
#

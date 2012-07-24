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
"""Control for Reading in floats

NOTE:  new syntax for custom events

NEW WAY:

SomeNewEvent, EVT_SOME_NEW_EVENT = \
     wx.lib.newevent.NewEvent()
SomeNewCommandEvent, EVT_SOME_NEW_COMMAND_EVENT = \
     wx.lib.newevent.NewCommandEvent()
#create the event
evt = SomeNewEvent(self.GetId(), attr1="hello", attr2=654)
[note that the id arg was left off of the PyWiki docs]

#post the event
wx.PostEvent(target, evt)

OLD WAY:  with subclassing PyEventBinder

EVT_FLOAT_CTRL_T = wx.NewEventType()
EVT_FLOAT_CTRL   = wx.PyEventBinder(EVT_FLOAT_CTRL_T, 1)
"""
import wx
import wx.lib.newevent

__all__ = ['FloatControl', 'EVT_FLOAT_CTRL']
FloatControlEvent, EVT_FLOAT_CTRL = wx.lib.newevent.NewCommandEvent()
#
# ---------------------------------------------------CLASS:  FloatControl
#
class FloatControl(wx.Panel):
    """FloatControl """
    __inParmDict = {
        'min':None,
        'max':None,
        'initial':None,
        'delta':1.0}
    def __init__(self, parent, id, **kwargs):
	"""Constructor for FloatControl."""
	#
	wx.Panel.__init__(self, parent, id) # , **kwargs)
	#
        #  Data
        #
        for parm, val in self.__inParmDict.iteritems():
            if kwargs.has_key(parm):
                val = kwargs.pop(parm)
            self.__setattr__(parm, val)
        if self.initial is not None:
            self.value = self.initial
        else:
            self.value = 0.0
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
        myToolTip = r"""FLOAT DATA ENTRY CONTROL

You can enter the number directly in the text
box or use the spinner to change the number.
The spinner increment is shown in the gray
box to the right of the spinner.
"""
        self.SetToolTipString(myToolTip)

	self.SetAutoLayout(True)
        self.SetSizerAndFit(self.sizer)
	#
	return
    #
    # ============================== Internal Methods
    #
    def __makeObjects(self):
        """Add interactors"""
        self.value_txt = wx.TextCtrl(self, wx.NewId(),
                                     value=str(self.value),
                                     style=wx.RAISED_BORDER| wx.TE_PROCESS_ENTER)
        self.delta_txt = wx.TextCtrl(self, wx.NewId(),
                                     value=str(self.delta),
                                     style=wx.RAISED_BORDER| wx.TE_PROCESS_ENTER)
        self.delta_txt.SetBackgroundColour( (230, 230, 230) )
        self.spin_but = wx.SpinButton(self, wx.NewId())
        self.spin_but.SetRange(-1,1)

        return

    def __makeBindings(self):
        """Bind interactors"""
        self.Bind(wx.EVT_TEXT_ENTER, self.OnValueTxt, self.value_txt)
        self.Bind(wx.EVT_TEXT_ENTER, self.OnDeltaTxt, self.delta_txt)
        self.Bind(wx.EVT_SPIN,       self.OnSpin,     self.spin_but)

        return

    def __makeSizers(self):
	"""Lay out the interactors"""

	self.sizer = wx.BoxSizer(wx.HORIZONTAL)
	self.sizer.Add(self.value_txt, 1, wx.EXPAND|wx.ALIGN_RIGHT)
	self.sizer.Add(self.spin_but,  0, wx.LEFT|wx.RIGHT, 5)
	self.sizer.Add(self.delta_txt, 1, wx.EXPAND|wx.ALIGN_LEFT)

	return

    def __showValues(self):
        """Show current values in the controls"""
        self.value_txt.ChangeValue(str(self.value))
        self.delta_txt.ChangeValue(str(self.delta))

        return
    #
    # ============================== API
    #
    #                     ========== *** Methods
    #
    def SetValue(self, val):
        """Set the value of the interactor"""
        self.value = val

        self.__showValues()

        return

    def SetDelta(self, delta):
        """Set the value of the interactor"""
        self.delta = delta

        self.__showValues()

        return
    #                     ========== *** Event Callbacks
    def OnSpin(self, evt):
        """Callback for spin button"""
        sval = self.spin_but.GetValue() # int

        try:
            tmpVal = self.value + sval * self.delta
            if self.min is not None:
                if tmpVal < min:
                    msg = 'Requested spinner value is less than minimum specified: %f' % (min)
                    wx.MessageBox(msg)
                    tmpVal = self.value
            if self.max is not None:
                if tmpVal > max:
                    msg = 'Requested spinner value is greater than maximum specified: %f' % (max)
                    wx.MessageBox(msg)
                    tmpVal = self.value
            self.value = tmpVal
            evt = FloatControlEvent(self.GetId(), floatValue=self.value)
            wx.PostEvent(self, evt)
            self.spin_but.SetValue(0)
        except Exception as e:
            msg = 'Failed to set float value from spinner\n%s' % (str(e))
            wx.MessageBox(msg)
            pass

        self.__showValues()

        return

    def OnValueTxt(self, evt):
        """Callback for value_txt choice"""
        try:
            tmpVal = float(self.value_txt.GetValue())
            if self.min is not None:
                if tmpVal < min:
                    msg = 'Requested spinner value is less than minimum specified: %f' % (min)
                    wx.MessageBox(msg)
                    tmpVal = self.value
            if self.max is not None:
                if tmpVal > max:
                    msg = 'Requested spinner value is greater than maximum specified: %f' % (max)
                    wx.MessageBox(msg)
                    tmpVal = self.value
            self.value = tmpVal
            evt = FloatControlEvent(self.GetId(), floatValue=self.value)
            wx.PostEvent(self, evt)
        except Exception as e:
            msg = 'Failed to set float value from text\n%s' % str(e)
            wx.MessageBox(msg)
            pass
        self.__showValues()

        return

    def OnDeltaTxt(self, evt):
        """Callback for delta_txt choice"""
        try:
            self.delta = float(self.delta_txt.GetValue())
        except Exception as e:
            msg = '\n%s' % str(e)
            wx.MessageBox(msg)
            pass

        self.__showValues()

        return
    pass # end class
#
# -----------------------------------------------END CLASS:  FloatControl

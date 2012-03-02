#! /usr/bin/env python
# ============================================================
# Copyright (c) 2012, Lawrence Livermore National Security, LLC. 
# Produced at the Lawrence Livermore National Laboratory. 
# Written by Joel Bernier <bernier2@llnl.gov> and others. 
# LLNL-CODE-529294. 
# All rights reserved.
# 
# This file is part of HEXRD. For details, see https://github.com/joelvbernier/hexrd.
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
"""Module for handling fit parameters
"""
# ---------------------------------------------------CLASS:  FitParams
#
class FitParams:
    """FitParams"""
    def __init__(self, names, values, rangeMin, rangeMax):
        """Constructor for FitParams

        INPUTS 
        
        names    -- a list of names for each parameter
        values   -- list of initial values
        rangeMin -- minimum values
        rangeMax -- maximum values
"""
        #
	self.numParam = len(names)
        #
        self.params = []
        for i in range(self.numParam):
            self.params.append(
                Param(names[i], values[i], rangeMin[i], rangeMax[i])
                )
            pass
        #
        self.paramDict = dict(zip([p.name for p in self.params], self.params))

        return
    
    def __iter__(self):
        return self.params.__iter__()
    #
    # ============================== API
    #
    def getNumParam(self):  return self.numParam

    def getParam(self, name):  return self.paramDict[name]

    def setProp(self, name, 
                value=None, min=None, max=None, active=None):
        
        """Set a property of the named fit parameter

        INPUTS 

        name -- the name of the parameter to alter
        value, min, max, active (boolean) -- one or more properties to be set
"""
        self.paramDict[name].setProp(value=value, min=min, max=max, active=active)

        return

    def getProp(self, name, propName):
        """Return value of requested property

        INPUTS 

        name     -- the name of the parameter to query
        propName -- the name of the property to get

        OUTPUTS

        1.  the value of parameter property
"""
        return self.paramDict[name].getProp(propName)
    #
    pass  # end class
#
# -----------------------------------------------END CLASS:  FitParams
# ---------------------------------------------------CLASS:  Param
#
class Param:
    """Param -- class for a single fit parameter.
"""
    def __init__(self, name, value, min, max):
        """Constructor for Param"""
        #
	self.name  = name
        self.value = value
        self.min   = min
        self.max   = max
        #
        return
    #
    # ============================== API
    #
    def setProp(self,
                value=None, min=None, max=None, active=None):
        
        """Set a property of the  fit parameter

        INPUTS 

        value, min, max, active (boolean) -- one or more properties to be set
"""
        if value  is not None:  self.value  = value
        if min    is not None:  self.min    = min
        if max    is not None:  self.max    = max
        if active is not None:  self.active = active

        return

    def getProp(self, propName):
        
        """Get a property by name

        INPUTS 

        propName -- name of the property
"""
        return getattr(self, propName)
    #
    pass  # end class
#
# -----------------------------------------------END CLASS:  Param

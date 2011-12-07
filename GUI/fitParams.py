#! /usr/bin/env python
#
#  $Id: fitParams.py 532 2010-04-27 16:40:03Z boyce6 $
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

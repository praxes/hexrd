#! /usr/bin/env python
# DO-NOT-DELETE revisionify.begin() 
#
#   Copyright (c) 2007-2009 Lawrence Livermore National Security,
#   LLC. Produced at the Lawrence Livermore National Laboratory (Nathan
#   Barton <barton22@llnl.gov>) CODE-OCEC-08-104.
#   
#   Please also read the file NOTICES.
#   
#   This file is part of the mdef package (version 0.2) and is
#   free software: you can redistribute it and/or modify it under the
#   terms of the GNU Lesser General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#   
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   Lesser General Public License for more details.
#   
#   A copy of the GNU Lesser General Public License may be found in the
#   file NOTICES. If this file is missing, see
#   <http://www.gnu.org/licenses/>.
#
# DO-NOT-DELETE revisionify.end() 
#
"""Module for associating units with scalar quantities

This module has been modified from its original form
by removing the call to the "units" executable and restricting
the units to only those used by the heXRD package.

"""

import math

__all__ = ['valWUnit', 'toFloat', 'valWithDflt']

# centralized unit types; chosen to have a match in the units command
energyUN = "ENERGY"
lengthUN = "LENGTH"
angleUN = "ANGLE"
#
# mapping to centralized unit types, for convenience
#
uTDict = {
    "length"     : lengthUN,
    "angle"      : angleUN,
    "energy"     : energyUN,
    }


class UNames(object):
    """Units used in this module"""
    degrees = 'degrees'
    radians = 'radians'

    m = 'm'
    mm = 'mm'
    meter = 'meter'
    angstrom = 'angstrom'

    J = 'J'
    
    pass

cv_dict = {
    (UNames.degrees, UNames.radians): math.pi/180.0,
    (UNames.radians, UNames.degrees): 180/math.pi,
    
    (UNames.m, UNames.mm):  1.0e3,
    (UNames.m, UNames.meter):  1.0,
    (UNames.m, UNames.angstrom):  1.0e10,
    
    (UNames.meter, UNames.mm):  1.0e3,
    (UNames.meter, UNames.m):  1.0,
    (UNames.meter, UNames.angstrom):  1.0e10,
    
    (UNames.mm, UNames.m):  1.0e-3,
    (UNames.mm, UNames.meter):  1.0e-3,
    (UNames.mm, UNames.angstrom):  1.0e7,
    
    (UNames.angstrom, UNames.m):  1.0e-10,
    (UNames.angstrom, UNames.meter):  1.0e-10,
    (UNames.angstrom, UNames.mm):  1.0e-7,
    
    }

class valWUnit(object):
    "Value with units"""
    def __init__(self, name, unitType, value, unit):
        """Initialization

        INPUTS
        name
           (str) name of the item
        unitType
           (str) class of units, e.g. 'length', 'angle'
        value
           (float) numerical value 
        unit
           (str) name of unit
"""
        self.name = name
        if uTDict.has_key(unitType):
            self.uT   = uTDict[unitType]
        else:
            # trust that unitType is correct -- may be a combined type
            self.uT   = unitType
            
        self.value = value
        self.unit = unit
        #
        #  Original checked if unit is of unitType
        #
        pass # end init
            
    def __str__(self):
        tmpl = """item named "%s" representing %g %s""" 
        return tmpl % (self.name, self.value, self.unit)
    
    def __repr__(self):
        tmpl = 'valWUnit("%s","%s",%s,"%s")' 
        return tmpl % (self.name, self.uT, self.value, self.unit)
    
    def __mul__(self, other):
        if isinstance(other, float):
            new = valWUnit(self.name, self.uT, self.value*other, self.unit)
            return new
        elif isinstance(other, valWUnit):
            new = valWUnit('%s_times_%s' % (self.name,other.name),
                           '%s %s' % (self.uT,other.uT),
                           self.value*other.value,
                           '(%s)*(%s)' % (self.unit,other.unit)
                           )
            # really need to put in here something to resolve new.uT
            return new
        else:
            raise RuntimeError("mul with unsupported operand")
        
    def __add__(self, other):
        if isinstance(other, float):
            new = valWUnit(self.name, self.uT, self.value + other, self.unit)
            return new
        elif isinstance(other, valWUnit):
            new = valWUnit(self.name, self.uT, self.value + other.getVal(self.unit), self.unit)
            return new
        else:
            raise RuntimeError("add with unsupported operand")
        
    def __sub__(self, other):
        if isinstance(other, float):
            new = valWUnit(self.name, self.uT, self.value - other, self.unit)
            return new
        elif isinstance(other, valWUnit):
            new = valWUnit(self.name, self.uT, self.value - other.getVal(self.unit), self.unit)
            return new
        else:
            raise RuntimeError("add with unsupported operand")

    def _convert(self, toUnit):
        """Convert unit"""
        if self.unit == toUnit:
            return self.value
        #
        #  Needs conversion
        #
        try:
            return cv_dict[(self.unit, toUnit)]*self.value
        except:
            raise RuntimeError("Unit conversion not recognized")
        
    def isLength(self):
        """Return true if quantity is a length"""
        retval = self.uT == uTDict['length']
        return retval
    def isAngle(self):
        """Return true if quantity is an angle"""
        retval = self.uT == uTDict['angle']
        return retval
    def isEnergy(self):
        """Return true if quantity  is an energy"""
        retval = self.uT == uTDict['energy']
        return retval
    
    def getVal(self, toUnit):
        """Return value in requested units

        INPUTS
        
        toUnit
           (str) requested unit for output
"""
        return self._convert(toUnit)
    
def toFloat(val, unitName):
    """Return the raw value of the object

    INPUTS
    
    val
       (float|valWUnit) object with value
    unitName
       (str) name of unit

    This function returns the raw value of the object, ignoring the
    unit, if it is numeric or converts it to the requested units and
    returns the magnitude if it is a valWUnit instance.

    For example:
    
    >>> print toFloat(1.1, 'radians')
    1.1
    >>> v = valWUnit('vee', 'angle', 1.1, 'radians')
    >>> print toFloat(v, 'degrees')
    63.0253574644
        
    """
    toFloatScalar = lambda v, u: v.getVal(u) if hasattr(v, 'getVal') else v
    
    if hasattr(val, '__len__'):
        retval = map(lambda x: toFloatScalar(x, unitName), val)
    else:
        retval = toFloatScalar(val, unitName)
    return retval

def valWithDflt(val, dflt, toUnit=None):
    """Return value or default value"""
    
    retval = val
    if retval is None:
        retval = dflt
    if not toUnit is None:
        retval = toFloat(retval, toUnit)
    return retval
    
if __name__ == '__main__':
    #
    #  doc testing
    #
    import doctest
    print "running doctest"
    doctest.testmod()
    #
    #  other tests
    #
    def testConversions():
        print '===== Testing unit conversions ...'
        print '..... angles:'
        v = valWUnit('180d', 'angle', 180.0, 'degrees')
        print v
        print '   in degrees:', v.getVal('degrees')
        print '   in radians: ', v.getVal('radians')

        print '..... lengths:'
        ulist = ['m', 'mm', 'meter', 'angstrom']
        v = valWUnit('one meter', 'length', 1.0, 'meter')
        print v
        for u in ulist:
            print '   in ', u, ': ', v.getVal(u)
        return

    testConversions()
    
    pass # end of tests

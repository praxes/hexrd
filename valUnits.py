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
# for quantities with units: valUnits.py
#
# see also
#	~/working/materials/mpsiForms.py
# 	~/working/materials/iron/shared_c/matVals.py
#	~/working/materials/iron/shared_c/make_cmat3MBTR_bdivk.py
#
# examples:
#	from valUnits import *
#
#	# common units, like pressure, have short names in uTDict:
#	c11 = valWUnit("c11", "p", 210., "GPa")
#	c11.getVal("Mbar")
#	c11.getVal(bdivk)
#
#	# oddball units are fine too:
#	ab_c = valWUnit("ab_c", stressUN+"/"+temperatureUN, 1., "Mbar/degK")
#	# but may not work properly with dictionary creation

# from scipy import *
import os, sys, re

if __name__ != '__main__':
    debug = 0

unitsCommand = 'units'
if os.path.exists('/sw/bin/units'):
    unitsCommand = '/sw/bin/units'
#
import fileUtil
haveGenerics = False
try:
    unitsVersion = map(int, fileUtil.getFromPipe('%s --version' % (unitsCommand), werr=True).split()[-1].split('.'))
    haveGenerics = unitsVersion[0] >= 1 and unitsVersion[1] >= 85
except:
    'do not bother with fall-back'

# centralized unit types; chosen to have a match in the units command
forceUN = "FORCE"
areaInvUN = "1/AREA"
energyUN = "ENERGY"
timeUN = "TIME"
massUN = "MASS"
lengthUN = "LENGTH"
angleUN = "ANGLE"
volumeUN = "VOLUME"
temperatureUN = "TEMPERATURE"
rateUN = "FREQUENCY" # strainrate
stressUN = "STRESS"
pressureUN = "STRESS"
densityUN = "DENSITY"
heatCapMUN = "ENERGY/(MASS TEMPERATURE)"
heatCapUN = "ENERGY/(VOLUME TEMPERATURE)"
activationUN = "TEMPERATURE/STRESS"
velocityUN = "LENGTH/TIME"
nullUN = "none"
#
# mapping to centralized unit types, for convenience
uTDict = {
    "areaInv"    : areaInvUN,
    "length"     : lengthUN,
    "angle"      : angleUN,
    "energy"     : energyUN,
    "volume"     : volumeUN,
    "stress"     : stressUN,
    "pressure"   : stressUN,
    "pres"       : stressUN,
    "p"          : stressUN,
    "strainrate" : rateUN,
    "time"       : timeUN,
    "timerate"   : rateUN,
    "frequency"  : rateUN,
    "temperature" : temperatureUN,
    "rate"       : rateUN,
    "heatCapM"   : heatCapMUN,
    "heatCap"    : heatCapUN,
    "density"    : densityUN,
    "velocity"   : velocityUN,
    "none"       : nullUN,
    "null"       : nullUN
    }

# for convenience
GPa = "GPa"
permicrosecond = "1/microsecond"

def fillUnitSystem(uS):
    def toPow(this,pow):
        val = "(%s)^(%g)" % (this,pow)
        return val
    def divBy(a,b):
        'a/b'
        val = "(%s)/(%s)" % (a,b)
        return val
    uS[nullUN] = nullUN
    if (not uS.has_key(volumeUN)):
        uS[volumeUN] = toPow(uS[lengthUN],3)
    if (not uS.has_key(densityUN)):
        uS[densityUN] = divBy(uS[massUN],uS[volumeUN])
    if (not uS.has_key(forceUN)):
        uS[forceUN] = divBy(uS[massUN]+" "+uS[lengthUN], toPow(uS[timeUN],2))
    if (not uS.has_key(energyUN)):
        uS[energyUN] = divBy(uS[massUN]+" "+toPow(uS[lengthUN],2), toPow(uS[timeUN],2))
    if (not uS.has_key(stressUN)):
        uS[stressUN] = divBy(uS[massUN], toPow(uS[timeUN],2)+" "+uS[lengthUN])
    if (not uS.has_key(rateUN)):
        uS[rateUN] = divBy("1", uS[timeUN])
    if (not uS.has_key(heatCapMUN)):
        uS[heatCapMUN] = divBy(uS[energyUN], uS[massUN]+" "+uS[temperatureUN])
    if (not uS.has_key(heatCapUN)):
        uS[heatCapUN] =  divBy(uS[energyUN], uS[volumeUN]+" "+uS[temperatureUN])
    if (not uS.has_key(activationUN)):
        uS[activationUN] = divBy(uS[temperatureUN], uS[stressUN])
    if (not uS.has_key(areaInvUN)):
        uS[areaInvUN] = toPow(uS[lengthUN],-2)
    if (not uS.has_key(velocityUN)):
        uS[velocityUN] = divBy(uS[lengthUN],uS[timeUN])
    return

# unit systems
bdivk = {
    timeUN         : "microsecond",
    massUN         : "gram",
    lengthUN       : "centimeter",
    temperatureUN  : "degK",
    stressUN       : "Mbar"
    }
fillUnitSystem(bdivk)
cm_g_mus = bdivk
#
mum_g_ms = {
    timeUN         : "millisecond",
    massUN         : "gram",
    lengthUN       : "micrometer",
    temperatureUN  : "degK"
    }
fillUnitSystem(mum_g_ms)
#
mum_kg_s = {
    timeUN         : "second",
    massUN         : "kilogram",
    lengthUN       : "micrometer",
    temperatureUN  : "degK"
    }
fillUnitSystem(mum_kg_s)
#
# for mm_g_ms:
#	stress ends up being in GPa
mm_g_ms = {
    timeUN         : "millisecond",
    massUN         : "gram",
    lengthUN       : "millimeter",
    temperatureUN  : "degK"
    }
fillUnitSystem(mm_g_ms)
#
mm_kg_ms = {
    timeUN         : "millisecond",
    massUN         : "kilogram",
    lengthUN       : "millimeter",
    temperatureUN  : "degK"
    }
fillUnitSystem(mm_kg_ms)
#
# Aaron Wemhoff's units:
#  Unit of time = 10^-10 s
#  Unit of distance = 10^-6 m
#  Unit of mass = 10^-12 g.
wemhoff_uS = { # mum_picog_1e-10s
    timeUN         : "0.1 nanosecond",
    massUN         : "picogram",
    lengthUN       : "micrometer",
    temperatureUN  : "degK"
    }
fillUnitSystem(wemhoff_uS)
# checks:
# valWUnit('pres', stressUN, 1.0, 'Mbar').getVal(wemhoff_uS)
#	-> 1.0
# valWUnit('dens', densityUN, 1.0, 'gram/(cm^3)').getVal(wemhoff_uS)
# 	-> 1.0
um_pg_ws = wemhoff_uS

def callUnits(uV, toUnit):
    negative = False
    val = uV.value
    if uV.value < 0:
        negative = True
        val = -val
    #
    # Command needs double quotes on windows
    #
    #command = "%s '%s %s' '%s'" % (unitsCommand, val, uV.unit, toUnit)
    command = '%s "%s %s" "%s"' % (unitsCommand, val, uV.unit, toUnit)
    unitsPipe = os.popen(command)
    tmp = unitsPipe.readline()
    unitRetVal = None
    while tmp:
        if debug > 1 :
            print "tmp:", tmp
        tmpb = re.split("\*",tmp)
        if len(tmpb) > 1:
            unitRetVal = float(tmpb[1])
        tmp = unitsPipe.readline()
    status = unitsPipe.close()
    if status != None:
        print >> sys.stderr, "error %s from units command %s" % (status,command)
        raise RuntimeError, "unrecoverable error"
    if unitRetVal is None:
        print >> sys.stderr, "unitRelVal not set from units command %s" % (command)
        raise RuntimeError, "unrecoverable error"
    if negative:
        unitRetVal = -unitRetVal
    return unitRetVal

class valWUnit:
    "value with units"
    def __init__(self, *args):
        if len(args) == 4:
            name     = args[0]
            unitType = args[1]
            value    = args[2]
            unit     = args[3]
        elif len(args) == 3:
            name     = args[0]
            unitType = args[1]
            value    = args[2]
            if unitType == nullUN:
                unit     = nullUN
            else:
                print >> sys.stderr, "unable to figure out arguments"
                raise RuntimeError, "unrecoverable error"
        else:
            print >> sys.stderr, "unable to figure out arguments"
            raise RuntimeError, "unrecoverable error"
        self.name = name
        if uTDict.has_key(unitType):
            self.uT   = uTDict[unitType]
        else:
            # trust that unitType is correct -- may be a combined type
            self.uT   = unitType
        self.value = value
        self.unit = unit
        if (not unitType == nullUN):
            # call callUnits to make sure that unit is of unitType
            if haveGenerics:
                dummy = callUnits(self, self.uT)
            else:
                dummy = callUnits(self, self.unit)
            # if this did not raise an exception, the unit is of the right type
    def __repr__(self):
        return 'valWUnit("%s","%s",%s,"%s")' % (self.name, self.uT, self.value, self.unit)
    def __str__(self):
        return "%s %s" % (self.value, self.unit)
    def isLength(self):
        retval = self.uT == uTDict['length']
        return retval
    def isAngle(self):
        retval = self.uT == uTDict['angle']
        return retval
    def isEnergy(self):
        retval = self.uT == uTDict['energy']
        return retval
    def getVal(self, asUnit):
        import string
        "asUnit can be a string or a unit system dictionary"
        if isinstance(asUnit, dict):
            if asUnit.has_key(self.uT):
                toUnit = asUnit[self.uT]
            else:
                "attempt to use dictionary to make toUnit"
                toUnit = self.uT
                for item in asUnit.iteritems():
                    toUnit = string.replace(toUnit, item[0], item[1])
                #print >> sys.stderr, "key not found in uTDict %s" % self.uT
                #raise RuntimeError, "unrecoverable error"
        else:
            toUnit = asUnit
        
        if toUnit == nullUN:
            return self.value
        
        unitRetVal = callUnits(self, toUnit)
        return unitRetVal
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
            raise RuntimeError, "mul with unsupported operand"
    def __add__(self, other):
        if isinstance(other, float):
            new = valWUnit(self.name, self.uT, self.value + other, self.unit)
            return new
        elif isinstance(other, valWUnit):
            new = valWUnit(self.name, self.uT, self.value + other.getVal(self.unit), self.unit)
            return new
        else:
            raise RuntimeError, "add with unsupported operand"
    def __sub__(self, other):
        if isinstance(other, float):
            new = valWUnit(self.name, self.uT, self.value - other, self.unit)
            return new
        elif isinstance(other, valWUnit):
            new = valWUnit(self.name, self.uT, self.value - other.getVal(self.unit), self.unit)
            return new
        else:
            raise RuntimeError, "add with unsupported operand"

def makeUnitValDict(vUList, unitSystem):
    "make a dictionary with values and units, given a unit system"
    dict = {}
    for vU in vUList:
        dict[vU.name] = vU.getVal(unitSystem)
        dict[vU.name + "_unit"] = unitSystem[vU.uT]
    return dict

def makeValDict(vUList, unitSystem, prefix=None):
    "make a dictionary with values, given a unit system"
    dict = {}
    for vU in vUList:
        if prefix:
            dict[prefix + vU.name] = vU.getVal(unitSystem)
        else:
            dict[vU.name] = vU.getVal(unitSystem)
    return dict

def getConversion(uS, unitType, unit):
    tempValWUnit = valWUnit('temp',
                            unitType,
                            1.0,
                            unit)
    conversion = tempValWUnit.getVal(uS)
    return conversion

def toFloatScalar(val, unitName):
    if hasattr(val, 'getVal'):
        retval = val.getVal(unitName)
    else:
        retval = val
    return retval
    
def toFloat(val, unitName):
    if hasattr(val, '__len__'):
        retval = map(lambda x: toFloatScalar(x, unitName), val)
    else:
        retval = toFloatScalar(val, unitName)
    return retval

def valWithDflt(val, dflt, toUnit=None):
    retval = val
    if retval is None:
        retval = dflt
    if not toUnit is None:
        retval = toFloat(retval, toUnit)
    return retval
    

'''
import units
bar = units.scaled_unit("br","Pa",1.0e5) # not happy with 'bar' for some reason
Pa = units.unit('Pa')
'''

# ============================================================
# Copyright (c) 2007-2012, Lawrence Livermore National Security, LLC.
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
haveSidl = True
try:
  import sidlPyArrays
  if sidlPyArrays.type == "numpy":
    import numpy as num
  else:
    if sidlPyArrays.type == "numeric":
      import numeric as num
except:
  haveSidl = False
  import numpy as num

if haveSidl:

   import warnings
   warnings.simplefilter('ignore',DeprecationWarning)

   def getMem(shape,asOnes=False,asZeros=False,typeInt=False):
     if typeInt:
       if asOnes:
         mem = sidlPyArrays.createInt(shape, ordering='column', value=1)
       elif asZeros:
         mem = sidlPyArrays.createInt(shape, ordering='column', value=0)
       else:
         mem = sidlPyArrays.createInt(shape, ordering='column')
     else:
       if asOnes:
         mem = sidlPyArrays.createDouble(shape, ordering='column', value=1.0e0)
       elif asZeros:
         mem = sidlPyArrays.createDouble(shape, ordering='column', value=0.0e0)
       else:
         mem = sidlPyArrays.createDouble(shape, ordering='column')
     return mem
else:
   def getMem(shape,asOnes=False,asZeros=False,typeInt=False):
     if typeInt:
       dtype = 'i'
     else:
       dtype = 'float64'
     if asOnes:
       mem = num.ones(shape, dtype=dtype)
     elif asZeros:
       mem = num.zeros(shape, dtype=dtype)
     else:
       mem = num.empty(shape, dtype=dtype)
     return mem

dtypeI = getMem([1],typeInt=True).dtype.name
dtypeF = getMem([1],typeInt=False).dtype.name
def toArray(a):
  numa = num.array(a) # wasteful?
  kind = numa.dtype.kind
  if kind == 'i':
    #mem = num.array(a,dtype=dtypeI)
    #numa.dtype = dtypeI
    mem = num.array(numa, order='F', dtype=dtypeI)
  elif kind == 'f':
    #mem = num.array(a,dtype=dtypeF)
    #numa.dtype = dtypeF
    mem = num.array(numa, order='F', dtype=dtypeF)
  else:
    raise RuntimeError, 'unknown data type'

  return mem

def writeArray(fid, *args, **dargs):
  '''
  Print to file, pasting arrays together;
  eventually replace with numpy.savetxt?
  '''
  f = open(fid,'w')
  if dargs.has_key('format'):
    format = dargs['format']
  else:
    format='%24.16e'
  if dargs.has_key('sep'):
    sep = dargs['sep']
  else:
    sep=' '
  if dargs.has_key('layout'):
    layout = dargs['layout']
  else:
    layout = 'C'
  if layout == 'C' or layout == 'c':
    for iLine in range(args[0].shape[1]):
      for iArg in range(len(args)):
        arg = args[iArg]
        for iArgCol in range(arg.shape[0]):
          print >> f, format % (arg[iArgCol,iLine]),
          print >> f, sep,
      print >> f, "" # to end line
  elif layout == 'F' or layout == 'f':
    for iLine in range(args[0].shape[0]):
      for iArg in range(len(args)):
        arg = args[iArg]
        if len(arg.shape) > 1:
          for iArgCol in range(arg.shape[1]):
            print >> f, format % (arg[iLine,iArgCol]),
            print >> f, sep,
        else:
            print >> f, format % (arg[iLine]),
            print >> f, sep,
      print >> f, "" # to end line
  else:
    raise RuntimeError, 'bad layout specification, try C or F'
  f.close()
  return

def arrayToString(a):
  import StringIO
  s = StringIO.StringIO()
  num.savetxt(s, a)
  return s.getvalue().replace('\n',' ')

def structuredSort(order, things):
  '''
  sort things by order, return sorted things
  '''
  assert len(order) == len(things), 'order and things must have same length'
  orderType = type(order[0])
  if hasattr(things,'dtype') and len(things.dtype) > 0:
    'this is realy ugly'
    dtype = [('order', orderType)] + things.dtype.descr
    a = num.array([ tuple( [order[i]] + list(things[i]) ) for i in range(len(order))], dtype=dtype)
    aSorted = num.sort(a, order='order')
    retval = aSorted
  elif hasattr(things[0],'lower'):
    thingType = object
    dtype = [('order', orderType), ('thing', thingType)]
    a = num.array([(order[i], things[i]) for i in range(len(order))], dtype=dtype)
    aSorted = num.sort(a, order='order')
    retval = aSorted['thing']
  else:
    thingType = type(things[0])
    dtype = [('order', orderType), ('thing', thingType)]
    a = num.array([(order[i], things[i]) for i in range(len(order))], dtype=dtype)
    aSorted = num.sort(a, order='order')
    retval = aSorted['thing']
  return retval

def histoFit(data, nBins, plot=False):
  '''
  Fit data using histogramming;
  useful for stuff pulled in using DataThief
  '''
  xVals = data[:,0]
  xMin = xVals.min()
  xMax = xVals.max()
  # nBins = 700
  delX = (xMax-xMin)/nBins
  xBins = num.arange(xMin, xMax+0.5*delX, delX)
  xCen  = 0.5*(xBins[1:]+xBins[:-1])
  dataHistO, edges = num.histogram(data[:,0],bins=xBins)
  dataHistW, edges = num.histogram(data[:,0],bins=xBins,weights=data[:,1])
  dataHist = dataHistW/dataHistO
  # take care of empty bins
  notNan = num.logical_not(num.isnan(dataHist))
  xCen = xCen[notNan]
  dataHist = dataHist[notNan]

  if plot:
    import hexrd.plotwrap as plotwrap
    pw = plotwrap.PlotWrap()
    pw(data[:,0],data[:,1],style='rx')
    pw(xCen, dataHist, style='k-')

  return xCen, dataHist

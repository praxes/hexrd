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

import numpy as num
import copy
import math

from hexrd.quadrature import q1db
#from hexrd.quadrature import q2db
#from scipy.optimize import fsolve
from scipy.optimize import leastsq




class Peak1DAtLoc:
  """
  base class for 1D peak shapes at fixed location;
  fixed that is unless newCenter is passed to the __call__ method
  """
  def __init__(self, centers, xVecDflt=None):
    """
    If __init__ is called with a list, then put one peak at each location
    """
    self.setCenters(centers)
    self.bRef = num.mean(self.centers) # does not change even if centers changed later
    self.setXVecDflt(xVecDflt)
    return
  def setXVecDflt(self, xVecDflt):
      if hasattr(xVecDflt, '__len__'):
          assert len(xVecDflt) == self.getNParams(),\
              'xVecDflt wrong length'
          self.xVecDflt = copy.deepcopy(xVecDflt)
      else:
          self.xVecDflt = xVecDflt
      return
  def getNParams(self):
    raise NotImplementedError
  def guessXVec(self, xs, vals, w=None):
      raise NotImplementedError
  def getNPeaks(self):
    return len(self.centers)
  def setCenters(self, centers):
    if hasattr(centers, '__len__'):
      self.centers = centers
    else:
      self.centers = [centers]
    return
  def __call__(self, xVec, p):
      # if newCenters is not None:
      #     self.setCenters(newCenters)
      if xVec is None or len(xVec) == 0:
          assert self.xVecDflt is not None,\
              'xVec is empty and xVecDflt is None'
          retval = self.eval(self.xVecDflt, p)
      else:
          retval = self.eval(xVec, p)
      return retval
  def eval(self, xVec, p):
    """
    xVec is parameters, p is positions
    """
    raise NotImplementedError
  def d_dx(self, xVec, p):
    'derivative of call with respect to xVec'
    raise NotImplementedError
  def d_dp(self, xVec, p):
    """derivative of call with respect to p
    assuming each eval depends only on its own point!
    """
    raise NotImplementedError
  def d_dCenters(self, xVec, p):
    'derivative of call with respect to centers'
    raise NotImplementedError
  def fitFloatingCenter(self, tThVals, intensityVals,
                        xVecGuess=None, centersGuess=None,
                        weights=4, tThWidth=None,
                        fitGoodnessTol=0.5):
      '''
      Note that centers are kept as they are -- if you want to
      actually change the centers of the function you need to call
      setCenters(cFit) after calling this function
      '''

      func = self

      if isinstance(weights, int):
          'interpret weights as number of quadrature points, need tThWidth'
          assert tThWidth is not None, \
              'need tThWidht for doing internal quadrature'
          assert len(tThVals.shape) == 1,\
              'tThVals wrong shape : '+str(tThVals.shape)
          quadr = weights
          xi1d, w1d = q1db.qLoc(quadr)
          'xi are in [0,1], so need to be centered'
          xi1d = xi1d - 0.5
          'and now make the right width for tTh bin sizes'
          xi1d = xi1d * tThWidth
          tThQP = num.tile(tThVals, (len(w1d),1)).T + num.tile(xi1d[:], (len(tThVals),1))
          w2d = num.tile(w1d, (len(tThVals),1))
      elif hasattr(weights, 'shape'):
          assert len(tThVals.shape) == 2,\
              'tThVals wrong shape : '+str(tThVals.shape)
          assert len(weights.shape) == 2,\
              'weights wrong shape'
          w2d = weights
          tThQP = tThVals
      else:
          raise RuntimeError, 'do not know what to do with weights of type '+str(type(weights))

      nPeaks = func.getNPeaks()

      if xVecGuess is None:
          #xVecGuess = func.xVecDflt
          xVecGuess = func.guessXVec(tThQP, intensityVals, w=w2d)
      if centersGuess is None:
          centersGuess = func.centers
      assert nPeaks == len(centersGuess), 'failed sanity check'
      xVec0 = num.hstack( (xVecGuess, centersGuess) )

      'keep centers so that can restore them later'
      centersRef = copy.deepcopy(func.centers)

      if len(xVec0) > len(intensityVals):
          raise RuntimeError, 'more DOF than equations'

      def _f_leastsq(x):
          # retval = num.empty(nBins)
          xFit = x[:-nPeaks]
          cFit = x[-nPeaks:] # copy.deepcopy(func.centers)
          func.setCenters(cFit)
          # evalIntensity = func(xFit, tThVals)
          funcQP = func(xFit, tThQP)
          evalIntensity = num.sum(funcQP * w2d, axis=1)
          retval = evalIntensity - intensityVals
          return retval
      def _df_leastsq(x):
          # retval = num.empty(nBins)
          xFit = x[:-nPeaks]
          cFit = x[-nPeaks:] # copy.deepcopy(func.centers)
          func.setCenters(cFit)
          retval = num.empty((len(intensityVals),len(x)))
          # funcQP = func(xFit, tThQP)
          d_evalQP_d_x       = func.d_dx(xFit, tThQP)
          d_evalQP_d_centers = func.d_dCenters(xFit, tThQP)
          for iX in range(0,len(xFit)):
              retval[:,iX]           = num.sum(d_evalQP_d_x[:,:,iX] * w2d, axis=1)
          for iC in range(0,len(cFit)):
              retval[:,len(xFit)+iC] = num.sum(d_evalQP_d_centers[:,:,iC] * w2d, axis=1)
          return retval

      x, ier = \
          leastsq(_f_leastsq, xVec0, Dfun=_df_leastsq)
      if not [1,2,3,4].count(ier):
          print >> sys.stderr, \
              'error code %d from leastsq' % (ier)
          raise RuntimeError, 'error from leastsq'

      xFit = x[:-nPeaks]
      cFit = x[-nPeaks:] # copy.deepcopy(func.centers)
      'return evalIntensity so caller does not have to monkey with changing centers back and forth'
      func.setCenters(cFit)
      residual = _f_leastsq(x)
      residNorm = num.linalg.norm(residual, ord=num.inf)
      fitTol = num.linalg.norm(intensityVals)*fitGoodnessTol
      if residNorm > fitTol :
          print 'fit might not be good enough : %g > %g' % (residNorm, fitTol)
          # raise RuntimeError, 'fit is not good enough : %g > %g' % (residNorm, fitTol)
      evalIntensity = residual + intensityVals # func(xFit, tThVals)
      func.setCenters(centersRef)
      return xFit, cFit, evalIntensity

class PeakPV1DAtLoc(Peak1DAtLoc):
    """
    the pseudo-Voigt:
    f = A*( n*fl + (1 - n)*fg )
    """
    def __init__(self, *args, **kwargs):
        Peak1DAtLoc.__init__(self, *args, **kwargs)
        self.C0 = 4.0*math.log(2.0)
        self.C1 = 4.0
        self.Al = num.sqrt(self.C1)/num.pi
        self.Ag = num.sqrt(self.C0/num.pi)
        return
    def getNParams(self):
        nParams = 2 + len(self.centers) * 3
        return nParams
    def guessXVec(self, xs, vals, w=None):
        xVec = num.empty(self.getNParams())
        if len(self.centers) == 1:
            xGauss = getGaussNDParams([xs], v=vals, w=w)
            xVec[0] = xGauss[3]
            xVec[1] = 0.0e0
            iPeak = 0
            H = xGauss[1] # FWHM
            xVec[2+iPeak*3]   = xGauss[2] * H / self.Ag
            xVec[2+iPeak*3+1] = H
            xVec[2+iPeak*3+2] = 0.0 # xn
        else:
            maxV = vals.max()
            minV = vals.min()
            width = (xs.max()-xs.min())/len(self.centers)
            xVec[0] = minV
            xVec[1] = 0.0e0
            for iPeak in range(len(self.centers)):
                xVec[2+iPeak*3]   = maxV * width # A
                xVec[2+iPeak*3+1] = 0.25 * width # H
                xVec[2+iPeak*3+2] = 0.0 # xn
        return xVec
    def eval(self, xVec, p):
        B  = xVec[0]
        dB = xVec[1]
        retval =  num.tile(B, num.shape(p))
        retval += dB * (p - self.bRef)
        for iPeak, center in enumerate(self.centers):
            A  = xVec[2+iPeak*3]
            H  = xVec[2+iPeak*3+1]
            xn = xVec[2+iPeak*3+2]
            n  = 0.5 + 0.5 * num.tanh(xn)
            delx = (p - center)

            fl = self.Al/H * 1 / (1 + self.C1*delx**2 / H**2)
            fg = self.Ag/H * num.exp(-self.C0*delx**2 / H**2)

            retval += A * (n*fl + (1.0 - n)*fg)
        return retval
    def d_dx(self, xVec, p):
        'allow p to be general shape'
        ps = p.shape
        p = p.reshape(p.size)
        retval = num.zeros((p.size, len(xVec)))

        # B  = xVec[0]
        # retval =  num.tile(B, num.shape(p))
        retval[:,0] = 1.0e0

        # dB = xVec[1]
        # retval += dB * (p - self.bRef)
        retval[:,1] = p - self.bRef

        for iPeak, center in enumerate(self.centers):
            A  = xVec[2+iPeak*3]
            H  = xVec[2+iPeak*3+1]
            xn = xVec[2+iPeak*3+2]
            n  = 0.5 + 0.5 * num.tanh(xn)
            delx = (p - center)

            fl = self.Al/H * 1 / (1 + self.C1*delx**2 / H**2)
            fg = self.Ag/H * num.exp(-self.C0*delx**2 / H**2)

            #####################
            # some intermediate partials
            dfl_dH = ( 2.0 * self.C1 * delx**2 * fl/(self.Al*H) - 1.0 ) * (fl / H)
            dfg_dH = ( 2.0 * self.C0 * delx**2 / H**2 - 1.0 ) * (fg / H)

            dfl_dx0 = 2.0 * self.C1 * delx * fl**2 / (self.Al * H)
            dfg_dx0 = 2.0 * self.C0 * delx * fg / H**2

            # mixing parm
            dn_dxn = 0.5 * (1.0/num.cosh(xn)**2)
            #####################

            # amplitude
            df_dA = (n*fl + (1.0 - n)*fg)

            # FWHM
            df_dH = A * ( n*dfl_dH + (1.0 - n)*dfg_dH )

            # mixing parm
            df_dxn = A * ( fl - fg ) * dn_dxn

            # assign jacobian vals
            retval[:,2+iPeak*3]   = df_dA
            retval[:,2+iPeak*3+1] = df_dH
            retval[:,2+iPeak*3+2] = df_dxn

        retval = retval.reshape(num.hstack((ps,len(xVec))))
        return retval
    def d_dp(self, xVec, p):
        retval = num.zeros(p.shape)

        # B  = xVec[0]
        dB = xVec[1]
        # retval =  num.tile(B, num.shape(p))
        # retval += dB * (p - self.bRef)
        retval = dB

        for iPeak, center in enumerate(self.centers):
            A  = xVec[2+iPeak*3]
            H  = xVec[2+iPeak*3+1]
            xn = xVec[2+iPeak*3+2]
            n  = 0.5 + 0.5 * num.tanh(xn)
            delx = (p - center)

            fl = self.Al/H * 1 / (1 + self.C1*delx**2 / H**2)
            fg = self.Ag/H * num.exp(-self.C0*delx**2 / H**2)

            dfl_dp = -(2.0 * self.C1) * delx * fl**2 / (self.Al*H)
            dfg_dp = -(2.0 * self.C0) * delx * fg / H**2

            retval += A * ( n*dfl_dp + (1.0 - n)*dfg_dp )
        return retval
    def d_dCenters(self, xVec, p):
        ps = p.shape
        retval = num.zeros((p.size,len(self.centers)))

        for iPeak, center in enumerate(self.centers):
            A  = xVec[2+iPeak*3]
            H  = xVec[2+iPeak*3+1]
            xn = xVec[2+iPeak*3+2]
            n  = 0.5 + 0.5 * num.tanh(xn)
            delx = (p - center).flatten()

            fl = self.Al/H * 1 / (1 + self.C1*delx**2 / H**2)
            fg = self.Ag/H * num.exp(-self.C0*delx**2 / H**2)

            dfl_dc = (2.0 * self.C1) * delx * fl**2 / (self.Al*H)
            dfg_dc = (2.0 * self.C0) * delx * fg / H**2

            retval[:,iPeak] = A * ( n*dfl_dc + (1.0 - n)*dfg_dc )

        retval = retval.reshape(num.hstack((ps,len(self.centers))))
        return retval

class PeakLorentzian1DAtLoc(Peak1DAtLoc):
    def __init__(self, *args, **kwargs):
        Peak1DAtLoc.__init__(self, *args, **kwargs)
        return
    def getNParams(self):
        """2 parameters for background, 2 for intensity and width of each peak"""
        nParams = 2 + len(self.centers) * 2
        return nParams
    def guessXVec(self, xs, vals, w=None):
        # guessXVec(self, widths, mins, maxs)
        xVec = num.empty(self.getNParams())
        xVec[0] = num.min(vals) # num.mean(mins)
        xVec[1] = 0.0e0
        width = (xs.max()-xs.min())/len(self.centers)
        for iPeak in range(len(self.centers)):
            xVec[2+iPeak*2]   = num.max(vals) # maxs[iPeak]
            xVec[2+iPeak*2+1] = 0.24 * width # 0.25 * widths[iPeak]
        return xVec
    def eval(self, xVec, p):
        B  = xVec[0]
        dB = xVec[1]
        retval =  num.tile(B, num.shape(p))
        retval += dB * (p - self.bRef)
        for iPeak, center in enumerate(self.centers):
            A = xVec[2+iPeak*2]
            w = xVec[2+iPeak*2+1]
            gamma=w/2.
            dist = (p - center)
            retval += A * (gamma)**2 / ( dist*dist + (gamma)**2  )#Reformulated so the amplitude given is the function max, DP 2/19/16
        return retval
    def d_dx(self, xVec, p):
        'allow p to be general shape'
        ps = p.shape
        p = p.reshape(p.size)
        retval = num.zeros((p.size, len(xVec)))

        # B  = xVec[0]
        # retval =  num.tile(B, num.shape(p))
        retval[:,0] = 1.0e0

        # dB = xVec[1]
        # retval += dB * (p - self.bRef)
        retval[:,1] = p - self.bRef

        for iPeak, center in enumerate(self.centers):
            A = xVec[2+iPeak*2]
            w = xVec[2+iPeak*2+1]
            gamma=w/2.
            dist = (p - center)

            # retval += A * num.exp( self.expFact * ( dist * dist) )
            evalLtz = (w)**2 / ( dist*dist + (w)**2 )
            retval[:,2+iPeak*2]   =  evalLtz
            retval[:,2+iPeak*2+1] =0.5* A * (2.*dist*dist*w)/(dist*dist + (w)**2 )**2 #Reformulated so the amplitude given is the function max, DP 2/19/16

        retval = retval.reshape(num.hstack((ps,len(xVec))))
        return retval
    def d_dp(self, xVec, p):
        retval = num.zeros(p.shape)

        # B  = xVec[0]
        dB = xVec[1]
        # retval =  num.tile(B, num.shape(p))
        # retval += dB * (p - self.bRef)
        retval = dB

        for iPeak, center in enumerate(self.centers):
            A = xVec[2+iPeak*2]
            w = xVec[2+iPeak*2+1]
            gamma=w/2.
            dist = (p - center)
            #evalLtz = (0.5*w) / ( dist*dist + (0.5*w)**2  )
            # retval += A * num.exp( self.expFact * ( dist * dist) )
            retval += -1.0 * A * 4.0 * gamma**2 * dist*dist /(dist*dist+gamma**2)**2  #Reformulated so the amplitude given is the function max, DP 2/19/16
        return retval
    def d_dCenters(self, xVec, p):
        ps = p.shape
        retval = num.zeros((p.size,len(self.centers)))

        for iPeak, center in enumerate(self.centers):
            A = xVec[2+iPeak*2]
            w = xVec[2+iPeak*2+1]
            gamma=w/2
            dist = (p - center).flatten()
            #evalLtz = (0.5*w) / ( dist*dist + (0.5*w)**2  )
            #retval[:,iPeak] = A * 4.0 * evalLtz * evalLtz * dist * (1/w)
            retval[:,iPeak] += -1.0 * A * 4.0 * gamma**2 * dist*dist /(dist*dist+gamma**2)**2  #Reformulated so the amplitude given is the function max, DP 2/19/16            
            
        retval = retval.reshape(num.hstack((ps,len(self.centers))))
        return retval

class PeakGauss1DAtLoc(Peak1DAtLoc):

    def __init__(self, *args, **kwargs):
        linBG=True
        if kwargs.has_key('linBG'):
            linBG = kwargs.pop('linBG')
        self.linBG = linBG
        self.nPBG = 1
        if self.linBG: self.nPBG += 1
        Peak1DAtLoc.__init__(self, *args, **kwargs)
        self.expFact  = -4.0 * math.log(2.0)
        self.widthFact  = 2.0 *math.sqrt(2.0*math.log(2.0))
        return
    def getNParams(self):
        """2 parameters for background, 2 for intensity and width of each peak"""
        nParams = self.nPBG + len(self.centers) * 2
        return nParams
    def guessXVec(self, xs, vals, w=None):
        xVec = num.empty(self.getNParams())
        if len(self.centers) == 1:
            xGauss = getGaussNDParams([xs], v=vals, w=w)
            xVec[0] = xGauss[3]
            width = (xs.max()-xs.min())/len(self.centers)
            if self.linBG :
                xVec[1] = 0.0e0
            iPeak = 0
            xVec[self.nPBG+iPeak*2]   = xGauss[2]
            xVec[self.nPBG+iPeak*2+1] = xGauss[1]
        else:
            xVec[0] = num.min(vals)
            width = (xs.max()-xs.min())/len(self.centers)
            if self.linBG :
                xVec[1] = 0.0e0
            for iPeak in range(len(self.centers)):
                xVec[self.nPBG+iPeak*2]   = num.max(vals) # maxs[iPeak]
                xVec[self.nPBG+iPeak*2+1] = 0.25 * width # 0.25 * widths[iPeak]
        return xVec
    def eval(self, xVec, p):
        B  = xVec[0]
        if self.linBG:
            dB = xVec[1]
        retval =  num.tile(B, num.shape(p))
        if self.linBG:
            retval += dB * (p - self.bRef)
        for iPeak, center in enumerate(self.centers):
            A = xVec[self.nPBG+iPeak*2]
            w = xVec[self.nPBG+iPeak*2+1]
            sigma=1/self.widthFact*w
            dist = (p - center)
            retval += A * num.exp( -1.0*( dist * dist) /(2. * sigma**2) )
        return retval
    def d_dx(self, xVec, p):
        'allow p to be general shape'
        ps = p.shape
        p = p.reshape(p.size)
        retval = num.zeros((p.size, len(xVec)))

        # B  = xVec[0]
        # retval =  num.tile(B, num.shape(p))
        iXBase = 0
        retval[:,iXBase] = 1.0e0
        iXBase += 1

        if self.linBG:
            # dB = xVec[1]
            # retval += dB * (p - self.bRef)
            retval[:,iXBase] = p - self.bRef
            iXBase += 1

        for iPeak, center in enumerate(self.centers):
            A = xVec[self.nPBG+iPeak*2]
            w = xVec[self.nPBG+iPeak*2+1]
            sigma=1/self.widthFact*w
            dist = (p - center)

            # retval += A * num.exp( self.expFact * ( dist * dist) )
            evalExp = num.exp( self.expFact * ( dist * dist) )
            retval[:,iXBase+iPeak*2]   = evalExp
            retval[:,iXBase+iPeak*2+1] = A * evalExp * self.expFact * 2.0e0 * dist * (-1.0*dist/w)
        iXBase += len(self.centers)*2
        assert iXBase == len(xVec), 'bookkeeping error'

        retval = retval.reshape(num.hstack((ps,len(xVec))))
        return retval
    def d_dp(self, xVec, p):
        retval = num.zeros(p.shape)

        # B  = xVec[0]
        # retval =  num.tile(B, num.shape(p))
        # retval += dB * (p - self.bRef)
        if self.linBG:
            dB = xVec[1]
            retval = dB

        for iPeak, center in enumerate(self.centers):
            A = xVec[self.nPBG+iPeak*2]
            w = xVec[self.nPBG+iPeak*2+1]
            dist = (p - center) / w
            # retval += A * num.exp( self.expFact * ( dist * dist) )
            retval += A * num.exp( self.expFact * ( dist * dist) ) * self.expFact * 2.0e0 * dist * (1/w)
        return retval
    def d_dCenters(self, xVec, p):
        'allow p to be general shape'
        ps = p.shape
        #retval = num.zeros(num.hstack((p.shape,len(self.centers))))
        retval = num.zeros((p.size,len(self.centers)))

        for iPeak, center in enumerate(self.centers):
            A = xVec[self.nPBG+iPeak*2]
            w = xVec[self.nPBG+iPeak*2+1]
            dist = ((p - center) / w).flatten()
            retval[:,iPeak] = -A * num.exp( self.expFact * ( dist * dist) ) * self.expFact * 2.0e0 * dist * (1/w)
        retval = retval.reshape(num.hstack((ps,len(self.centers))))
        return retval






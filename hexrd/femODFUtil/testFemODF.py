# ============================================================
# Copyright (c) 2007-2012, Lawrence Livermore National Security, LLC. 
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
import sys, os

try:
  import sidl.RuntimeException
  import femODF.FemMesh
  import femODF.FemRodrigues
  import femODF.FemHemisphere
  import femODF.DiscreteAggregate
  import utilIO.InMemOrFSPack
except:
  print >> sys.stderr, 'error doing imports; try sourcing /usr/apps/dlsmm/bin/setup.sh or the like'
  sys.exit(1)

import arrayUtil
from arrayUtil import getMem
from arrayUtil import num
if not arrayUtil.haveSidl: 
    raise 'error: haveSidl is false in arrayUtil'

from numpy import random
import math
nSampleVectors = 12
#
sampleNVectors = getMem([3,nSampleVectors])
sampleNVectors[:,:] = random.rand(3*nSampleVectors).reshape(3,nSampleVectors)
for iSV in range(nSampleVectors):
  mag = math.sqrt((sampleNVectors[:,iSV]*sampleNVectors[:,iSV]).sum())
  sampleNVectors[:,iSV] = sampleNVectors[:,iSV] / mag

# sys.path.append('./data')

testStrgDir = 'test_storage'
from femODFUtil import mRodrStorage
mRodrStorage.cleanStorage(storageDir=testStrgDir)

from testFemODF.mRodrUtil import convertFemODFCoo

def testMRodr():
  
  #import meshODFCub
  #mRodr = meshODFCub.getMRodr(0)
  from femODFUtil.mRodrStorage import getMesh, symmGroupString_cub
  mRodr = getMesh(symmGroupString_cub, 0, storageDir=testStrgDir)

  opM = mRodr.getM()
  Mcsr = M = convertFemODFCoo(opM)
  
  crystalVector = getMem([3])
  crystalVector[:] = [1.,0.,0.]
  #
  #
  try:
    odfPfProj = mRodr.getPoleProjection(crystalVector, sampleNVectors, True)
  except sidl.RuntimeException._Exception:
      (etype, eobj, etb) = sys.exc_info()
      # eobj is the SIDL exception object
      print >> sys.stderr, "Exception note: %s" % (eobj.exception.getNote())
      print >> sys.stderr, "Exception trace: %s" % (eobj.exception.getTrace())
      print >> sys.stderr, "About to call sys.exit due to exception"
      sys.exit(1)
  
  Pcsr = P = convertFemODFCoo(odfPfProj)
  # Pcsr.T * Pcsr # seems to work fine

  vProj = mRodr.getVolProj()
  vProjSum = vProj.sum()
  print >> sys.stdout, 'vProjSum is %g' % (vProjSum)

  (inElem, elemXi) = mRodr.find(num.array([0.1,0.1,-0.2]))
  print >> sys.stdout, 'found in element %d' % (inElem)

  return 0

def testMRodrPack():
  try: 

    from femODFUtil.mRodrStorage import getMesh, symmGroupString_cub

    mRodrA = getMesh(symmGroupString_cub, 0, storageDir=testStrgDir)
    vProjA = mRodrA.getVolProj()  
    
    filename = 'mRodr_cub.p'
    import fileUtil
    fileUtil.rmWild(filename)
  
    p = utilIO.InMemOrFSPack.InMemOrFSPack()
    p.initFS(filename)
    packer = utilIO.Packer.Packer(p)
    mRodrA.pack(packer)
    p.closeFile()

    # cheat and load the shelf directly
    import shelve
    s = shelve.open(filename)
    print >> sys.stdout, 'shelf file has keys : %s' % (str(s.keys()))
    del s
  
    p = utilIO.InMemOrFSPack.InMemOrFSPack()
    p.initFS(filename)
    unp = utilIO.Unpacker.Unpacker(p)
    mRodrB = femODF.FemRodrigues.FemRodrigues()
    mRodrB.unpack(unp)
    vProjB = mRodrB.getVolProj()  
    if sum(abs(vProjA-vProjB)) > 1.0e-8:
      print >> sys.stdout, 'ERROR: vProj mismatch'
      return 1
    # call find just to check that things work
    (inElem, elemXi) = mRodrB.find(num.array([0.1,0.1,-0.2]))
    print >> sys.stdout, 'found in element %d' % (inElem)
    p.closeFile()

  except sidl.RuntimeException._Exception:
    (etype, eobj, etb) = sys.exc_info()
      # eobj is the SIDL exception object
    print >> sys.stderr, "Exception note: %s" % (eobj.exception.getNote())
    print >> sys.stderr, "Exception trace: %s" % (eobj.exception.getTrace())
    #print >> sys.stderr, "Ignoring this not-implemented exception for now -- I know it is not implemented!"
    return 1

  return 0

def testPlot():
  try:

    from femODFUtil.mRodrStorage import getMesh, symmGroupString_cub

    symmGroupString = symmGroupString_cub

    mRodr = getMesh(symmGroupString, 3, storageDir=testStrgDir)
    vProj = mRodr.getVolProj() # just to have something to plot
    from femODFUtil import dxWrap
    dxFS = dxWrap.DxFS(prefix='testPlotODF', overWrite=True) 
    print >> sys.stdout, 'handing mesh to dxFS'
    dxFS.setMesh(mRodr)
    print >> sys.stdout, 'attaching nodal data'
    dxFS.attachNodalData(vProj, 'A')
    dxFS.close()
    #
    # ... want to code edges of fundamental region and mechanism for putting in DX file?
    #
    dxCL = dxWrap.DxCL()
    dxCL.odf(dxFS.getFName(), symmGroupString, imageNamePrefix='test_mRodr')

  except sidl.RuntimeException._Exception:
    (etype, eobj, etb) = sys.exc_info()
      # eobj is the SIDL exception object
    print >> sys.stderr, "Exception note: %s" % (eobj.exception.getNote())
    print >> sys.stderr, "Exception trace: %s" % (eobj.exception.getTrace())
    #print >> sys.stderr, "Ignoring this not-implemented exception for now -- I know it is not implemented!"
    return 1
  
  return 0
def testDA():
  from femODFUtil.mRodrStorage import symmGroupString_cub
  symmGroupString = symmGroupString_cub # may ultimately be a better place to get this

  dagg = femODF.DiscreteAggregate.DiscreteAggregate()
  dagg.createFromUniform(512, 0)
  #
  crystalVectors = getMem([3,1])
  crystalVectors[:,:] = num.array([[1.,0.,0.]]).T
  #
  pFigVals = dagg.calculatePoleValues(crystalVectors, sampleNVectors, symmGroupString, 0.1)

  return 0

def testPoleFigCalc():
  from femODFUtil.mRodrStorage import symmGroupString_hex
  import orientations as ors
  import math
  
  mPfig = femODF.FemHemisphere.FemHemisphere()
  mPfig.autoMesh(0.01)

  #import fileUtil, os
  #quatsData = num.array(fileUtil.readFloatData(os.path.join('data','rolling.quat.txt'))).T
  nGrain = 128
  quatsBall = num.ones([4,nGrain])
  quatsBall[1:,:] = num.random.normal(loc=0.0, scale=0.04, size=[3,nGrain])
  for iQ in range(quatsBall.shape[1]):
    mag = num.sum(quatsBall[:,iQ]*quatsBall[:,iQ])
    quatsBall[:,iQ] = quatsBall[:,iQ] 
  quatsData = num.ones([4,nGrain])
  for iQ in range(quatsBall.shape[1]):
      qAboutC = ors.Quat(ors.RotInv(2.*math.pi*num.random.normal(),  (0.,0.,1.)))
      qBall   = ors.Quat(quatsBall[:,iQ])
      qFlip   = ors.Quat(ors.RotInv(math.pi/2., (num.random.normal(),num.random.normal(),0.)))
      q = qFlip * qBall * qAboutC
      quatsData[:,iQ] = q.q[:]
  symmGroupString = symmGroupString_hex

  dagg = femODF.DiscreteAggregate.DiscreteAggregate()
  dagg.setOrientations(quatsData)
  
  crystalVectors = getMem([3,1])
  crystalVectors[:,:] = num.array([[0.,0.,1.]]).T # c-axis
  #
  coords = mPfig.getCoords()
  from femODFUtil import pfig
  nVectors = pfig.sph2n(coords)
  pFigVals = dagg.calculatePoleValues(crystalVectors, nVectors, symmGroupString, 0.1)

  from femODFUtil import dxWrap
  dxFS = dxWrap.DxFS() #  prefix='testPfig', overWrite=True
  print >> sys.stdout, 'handing mesh to dxFS'
  dxFS.setMesh(mPfig)
  print >> sys.stdout, 'attaching nodal data'
  dxFS.attachNodalData(pFigVals[:,0], 'p')
  dxFS.close()
  #
  dxCL = dxWrap.DxCL()
  dxCL.pfig(dxFS.getFName(), imageNamePrefix='test_pfig')
  
  return 0

def main():
  if testMRodr():
    print >> sys.stderr, 'failed test testMRodr'
    sys.exit(1)
  else:
    print >> sys.stdout, 'passed test testMRodr'

  if testMRodrPack():
    print >> sys.stderr, 'failed test testMRodrPack'
    sys.exit(1)
  else:
    print >> sys.stdout, 'passed test testMRodrPack'

  if testPlot():
    print >> sys.stderr, 'failed test testPlot'
    sys.exit(1)
  else:
    print >> sys.stdout, 'passed test testPlot'

  if testDA():
    print >> sys.stderr, 'failed test testDA'
    sys.exit(1)
  else:
    print >> sys.stdout, 'passed test testDA'

  if testPoleFigCalc():
    print >> sys.stderr, 'failed test testPoleFigCalc'
    sys.exit(1)
  else:
    print >> sys.stdout, 'passed test testPoleFigCalc'

  return

if __name__ == '__main__':
  main()

######################################################################

'''

app-open texture-r $texin
fc-list $texin
include $mesh_dir/in.symm_group
#
fc-lop solve-odfdds_system
#
odf_np-full_form
write-sfc-odf_np fmt $output
'''





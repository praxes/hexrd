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
    from arrayUtil import num
    from arrayUtil import getMem
except:
  print >> sys.stderr, 'error doing imports; need python tools in PYTHONPATH'
  sys.exit(1)

try:
  import sidl.RuntimeException
  import femODF.FemMesh
  import femODF.FemRodrigues
  import femODF.ElemType
except:
  print >> sys.stderr, 'error doing imports; try sourcing /usr/apps/dlsmm/bin/setup.sh or the like'
  sys.exit(1)

from femODFUtil.mRodrUtil import makeRodrMeshBall, getStorageDir
import femODFUtil, arrayUtil

shelfFilename = 'meshODFBall.py_shelf'
def cleanStorage(storageDir=None):
    import fileUtil, os
    storageDirThis = getStorageDir(storageDir)
    storageFile = os.path.join(storageDirThis, shelfFilename)
    os.remove(storageFile)
    return

##############################

def getM(inFix='base'):
    from femODFUtil.mRodrUtil import mFromCC
    
    path     = femODFUtil.getPath()[0]
    rsFile   = os.path.join(path, 'p3_%s_unitRodr.rodr' % (inFix))
    connFile = os.path.join(path, 'p3_%s.con' % (inFix))
    
    m = mFromCC(connFile, rsFile)
    
    return m

######################################################################

class MBallStorage:
    '''
    unlike MRodrStorage, just stores the mesh, not femODF.FemRodrigues objects
    '''
    __elemType = femODF.ElemType.stdE4 # 4-point tet
    # __quadRule = num.array([32011, 22004]) # 11-pt quadrature in volume, 4-pt quadrature on surfaces
    __quadRule = arrayUtil.toArray([32015, 22004]) # 15-pt quadrature in volume, 4-pt quadrature on surfaces
    __baseMeshList = ['base','10x']
    def __init__(self, storageDir=None, m0Name='base'):
        import shelve
        self.storageDir = getStorageDir(storageDir)
        self.meshFilename = os.path.join(self.storageDir, 
                                         shelfFilename)
        self.s = shelve.open(self.meshFilename)
        self.m0Name = m0Name
        meshNum = 0
        meshDir = self.__meshDir(meshNum)
        if not self.s.has_key(meshDir):
            self.s[meshDir] = self.__retrieveMeshName(meshDir)
        self.s.close()
    def __meshDir(self, meshName):
        if isinstance(meshName,int):
            assert meshName >= 0, 'meshNum too small'
            retval = 'm_%d' % (meshName)
        elif isinstance(meshName,str):
            retval = 'm_%s' % (meshName)
        return retval
    def __call__(self, meshName):
        return self.getMesh(meshName)
    def __retrieveMesh(self, meshDir):
        import shelve
        self.s = shelve.open(self.meshFilename)
        m = self.s[meshDir]
        self.s.close()
        return m
    def __retrieveMeshName(self, meshName):
        if meshName in self.__baseMeshList:
            m = getM(inFix=meshName)
        elif meshName == 'm_0':
            m = getM(inFix=self.m0Name)
        else:
            raise RuntimeError, 'do not know how to get mesh named : '+str(meshName)
        print >> sys.stdout, 'pulled mesh designated '+str(meshName)
        return m
    def __storeMesh(self, m, meshDir):
        import shelve
        self.s = shelve.open(self.meshFilename)
        self.s[meshDir] = m
        self.s.close()
        return
    def getMesh(self, meshName):
        '''get a mesh with a specified level of mesh refinement; 
        may call itself recursively'''
        import shelve
        from femODFUtil.mRodrUtil import refineTetMesh
        from arrayUtil import num
        
        self.s = shelve.open(self.meshFilename)
        
        meshDir = self.__meshDir(meshName)
        if self.s.has_key(meshDir):
            print 'getMesh : retrieving meshName %s' % (meshName)
            m = self.__retrieveMesh(meshDir)
        elif isinstance(meshName, int):
            meshNum=meshName
            mRef = self.getMesh(meshNum-1) # always have in self.unp meshNum=0
            print 'getMesh : refining to get meshNum %d' % (meshNum)
            (connNew, coordsNew) = refineTetMesh(mRef['conn'], mRef['coords'])
            # mRodr = makeRodrMeshBall() # call elsewhere
            
            'pull new mesh coordinates onto the unit ball'
            mObj = femODF.FemMesh.FemMesh()
            mObj.init(self.__elemType, self.__quadRule, connNew, coordsNew)
            sConn1 = mObj.getSConn()
            sConn = sConn1 - 1
            surfNodes = num.unique(sConn)
            surfNodeMags = num.apply_along_axis(num.linalg.norm, 0, coordsNew[:,surfNodes])
            coordsNew[:,surfNodes] = coordsNew[:,surfNodes] / surfNodeMags
            
            m = {
                'ne' : connNew.shape[1], 
                'nn' : coordsNew.shape[1],
                'conn' : connNew,
                'coords' : coordsNew,
                }
            self.__storeMesh(m, meshDir)
        elif isinstance(meshName, str):
            m = self.__retrieveMeshName(meshName)
        else: 
            raise RuntimeError, 'do not know what to do with meshName : '+str(meshName)
        
        self.s.close()
        return m
    

__meshes = {}
#
def getMeshes(**keyArgs):
    import sys
    storageDir = None
    if keyArgs.has_key('storageDir'):
        storageDir = keyArgs['storageDir']
    if __meshes.has_key(storageDir):
        meshesThis = __meshes[storageDir]
    else:
        meshesThis = MBallStorage(**keyArgs)
        __meshes[storageDir] = meshesThis
    return meshesThis
def getBallMesh(meshNum, **keyArgs):
    meshes = getMeshes(**keyArgs)
    m      = meshes(meshNum)
    return m
def getMesh(symmGroupString, ballSize, meshNum, cRF=None, forceSlave=None, **keyArgs):
    '''
    cRF set to the identity if None is passed, use resetBallCRF to change it later
    '''
    
    m = getBallMesh(meshNum, **keyArgs)
    
    if cRF is None:
        cRF = arrayUtil.getMem((3,3))
        cRF[:,:] = num.eye(3)
    
    mRodr = makeRodrMeshBall(m['conn'], m['coords'], 
                             symmGroupString, ballSize, cRF, forceSlave=forceSlave)
    
    return mRodr

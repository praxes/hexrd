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
'''Wrapper for interacting with dx
'''

import tempfile
import os

from arrayUtil import getMem
from arrayUtil import num

dxLines = '"lines"'
dxCubes = '"cubes"'
dxTets  = '"tetrahedra"'
dxTris  = '"triangles"'
dxQuads = '"quads"'

def getElMap(elemType, quality=1):
    import femODF.ElemType
    if elemType == femODF.ElemType.stdL2:
        npMap     = num.array([1,2]) - 1
        dxEtype   = dxLines
        nelemMult = 1
        dxNnpe    = 2
    elif elemType == femODF.ElemType.stdL3:
        npMap = num.array([1,2,2,3]) - 1
        dxEtype   = dxLines
        nelemMult = 2
        dxNnpe    = 2
        
    elif elemType == femODF.ElemType.stdB8:
        npMap = num.array([1, 4, 2, 3, 5, 8, 6, 7]) - 1
        dxEtype   = dxCubes
        nelemMult = 1
        dxNnpe    = 8

    elif elemType == femODF.ElemType.stdE4:
        npMap = num.array([1, 2, 3, 4]) - 1
        dxEtype   = dxTets
        nelemMult = 1
        dxNnpe    = 4
    
    elif elemType == femODF.ElemType.stdE10:
        if quality > 0:
            npMap = num.array([
                    2, 3, 4, 8, 
                    6, 4, 5, 9, 
                    1, 2, 6, 7, 
                    7, 8, 9, 10, 
                    6, 9, 8, 4, 
                    6, 8, 2, 4, 
                    8, 9, 6, 7, 
                    8, 6, 2, 7]) - 1
            dxEtype   = dxTets
            nelemMult = 8
            dxNnpe    = 4
        else:
            npMap = num.array([1, 3, 5, 10]) - 1
            dxEtype   = dxTets
            nelemMult = 1
            dxNnpe    = 4

    elif elemType == femODF.ElemType.stdT3:
            npMap = num.array([1, 2, 3]) - 1
            dxEtype   = dxTris
            nelemMult = 1
            dxNnpe    = 3
    elif elemType == femODF.ElemType.stdT6:
        if quality > 0:
            npMap = num.array([
                    1, 2, 6,
                    2, 3, 4,
                    6, 4, 5,
                    2, 4, 6]) - 1
            dxEtype   = dxTris
            nelemMult = 4
            dxNnpe    = 3
        else:
            npMap = num.array([1, 3, 5]) - 1
            dxEtype   = dxTris
            nelemMult = 4
            dxNnpe    = 3

    elif elemType == femODF.ElemType.stdQ4:
        npMap = num.array([1, 4, 2, 3]) - 1
        dxEtype   = dxQuads
        nelemMult = 1
        dxNnpe    = 4
    elif elemType == femODF.ElemType.stdQ9:
        if quality > 0:
            npMap = num.array([
                    1, 8, 2, 9,
                    2, 9, 3, 4,
                    8, 7, 9, 6,
                    9, 6, 4, 5 ]) - 1
            dxEtype   = dxQuads
            nelemMult = 4
            dxNnpe    = 4
        else:
            npMap = num.array([1, 7, 3, 5]) - 1
            dxEtype   = dxQuads
            nelemMult = 1
            dxNnpe    = 4
    else:
        raise RuntimeError, 'unknown elemType'
    assert len(npMap) == nelemMult*dxNnpe, 'bad getElMap'
    return (dxEtype, npMap, nelemMult, dxNnpe)

class DxFS:
    'for writing dx files to the file system'
    def __init__(self, prefix=None, overWrite=False, clobber=None):
        'set clobber to True or False to override default'
        from fileUtil import rmDirF
        self.fo = None
        self.fnameList = []
        self.prefix = prefix
        if self.prefix is None:
            self.clobber = True
            self.dataDir = tempfile.mkdtemp()
            (fd, name) = tempfile.mkstemp(suffix='.dx', text=True)
            os.close(fd)
            fo = open(name, 'w')
        else:
            self.clobber = False
            name    = prefix+'.dx'
            dataDir = prefix+'_dxdata'
            if overWrite:
                if os.path.lexists(name):    os.remove(name)
                if os.path.lexists(dataDir): rmDirF(dataDir)
            fo = open(name, 'w')
            self.dataDir = dataDir
            os.mkdir(self.dataDir)
        if clobber is not None:
            self.clobber = clobber
        self.fo = fo
        self.fnameList.append(name)
        self.fieldBuffer = 'object "mesh" class field\n'
        self.dVals = {}
        self.haveData = False
        if self.prefix is None:
            # 'binary' not working yet; even when do float->double
            self.dVals['dataFormat'] = 'ascii' 
        else:
            self.dVals['dataFormat'] = 'ascii'
        return
    def __dumpData(self, data, dataName):
        dataFileName = os.path.join(self.dataDir,dataName+'.dat')
        #if self.prefix is not None:
        if self.dVals['dataFormat'] == 'ascii':
            data.tofile(dataFileName, sep="\t")
        elif self.dVals['dataFormat'] == 'binary':
            data.tofile(dataFileName) 
        else:
            raise RuntimeError, 'unknown data format'
        self.fnameList.append(dataFileName)
        return dataFileName
    def setMesh(self, mesh):
        'mesh should be an femODF.FemMesh (or extension of FemMesh, or implement an interface like it at least)'

        elemType = mesh.getElemType()
        (dxEType, npMap, nelemMult, dxNnpe) = getElMap(elemType)
        
        self.fieldBuffer += '  component "connections" "connectivity"\n'
        conn  = mesh.getConn()
        dxConn = conn[npMap,:] - 1
        if nelemMult > 1:
            connMut = dxConn.reshape(dxNnpe, conn.shape[1]*nelemMult, order='f')
            dxConn = connMut
        self.dVals['nnpe']         = dxConn.shape[0]
        self.dVals['nelm']         = dxConn.shape[1]
        self.dVals['nnpeBase']     = conn.shape[0]
        self.dVals['nelmBase']     = conn.shape[1]
        self.dVals['connDataFile'] = self.__dumpData(dxConn.T, 'connData') # transpose to get in correct order
        self.dVals['dxElemType']   = dxEType
        self.dVals['nelemMult']    = nelemMult
        print >> self.fo, '''
object "connectivity"  class array type int rank 1 shape %(nnpe)d items %(nelm)d %(dataFormat)s
data file %(connDataFile)s
attribute "element type" string %(dxElemType)s
attribute "dep" string "connections"
attribute "ref" string "positions"
''' % self.dVals

        self.fieldBuffer += '  component "positions" "nodal points"\n'
        coords = mesh.getCoords()
        self.dVals['nDim'] = coords.shape[0]
        self.dVals['nPoints'] = coords.shape[1]
        self.dVals['coordsDataFile'] = self.__dumpData(coords.T, 'coordsData') # transpose to get in correct order
        print >> self.fo, '''
object "nodal points" class array type float rank 1 shape %(nDim)s items %(nPoints)s %(dataFormat)s
data file %(coordsDataFile)s
''' % self.dVals
        
        return
    def attachElemData(self, edata, name, asData=True):

        nelemMult = self.dVals['nelemMult']
        if nelemMult > 1:
            edataMult = getMem([len(edata)*nelemMult])
            for iMult in range(nelemMult):
                edataMult[iMult::nelemMult] = edata[:]
            dataFileName = self.__dumpData(edataMult, 'eData'+'_'+name) 
        else:
            dataFileName = self.__dumpData(edata, 'eData'+'_'+name) 
        self.dVals['curData'] = name
        self.dVals['curFile'] = dataFileName

        print >> self.fo, '''
object "%(curData)s" class array type float rank 0  items %(nelm)d %(dataFormat)s
data file %(curFile)s
attribute "dep" string "connections"
''' % self.dVals
        if asData and not self.haveData:
            self.fieldBuffer += '  component "data" "%(curData)s"\n' % self.dVals
        else:
            self.fieldBuffer += '  component "%(curData)" "%(curData)s"\n' % self.dVals
        if asData: self.haveData = True

        return
    def attachNodalData(self, ndata, name, asData=True):

        dataFileName = self.__dumpData(ndata, 'nData'+'_'+name) 
        self.dVals['curData'] = name
        self.dVals['curFile'] = dataFileName

        print >> self.fo, '''
object "%(curData)s" class array type float rank 0  items %(nPoints)d %(dataFormat)s
data file %(curFile)s
attribute "dep" string "positions"
''' % self.dVals
        if asData and not self.haveData:
            self.fieldBuffer += '  component "data" "%(curData)s"\n' % self.dVals
        else:
            self.fieldBuffer += '  component "%(curData)s" "%(curData)s"\n' % self.dVals
        if asData: self.haveData = True

        return
    def close(self):
        self.fieldBuffer += '  attribute "name" string "mesh"\nend\n'
        print >> self.fo, self.fieldBuffer
        self.fo.close()
        return
    def getFName(self):
        return self.fnameList[0]
    def preserve(self):
        'return top-level file descriptor, name, and data directory name; do not then clobber files by default'
        self.clobber = False
        return (self.fo, self.fnameList[0], self.dataDir)
    def __del__(self):
        "cleanup temporary file(s)"
        from fileUtil import rmDirF
        if self.clobber:
            self.fo.close()
            for fn in self.fnameList:
                os.remove(fn)
            rmDirF(self.dataDir)
            
def getFromPipe(command):
    import os
    pipe = os.popen(command)
    val = pipe.readline()[:-1] # [:-1] drops end-of-line character
    pipe.close()
    return val

class DxCL:
    '''interface through the command line and file system'''
    def __init__(self):
        if (os.system('which dx') != 0):
            import sys
            print >> sys.stderr, "need dx to be in the path"
            raise RuntimeError, "unrecoverable error"
        self.dxExec = getFromPipe('which dx')
        import femODFUtil
        self.path = femODFUtil.getPath()[0]
        return
    def pfig(self, dxInFile, 
             imageNamePrefix='pfig', imageFormat='tiff', 
             greyScale=True, logThresh=None,
             dataName='data',
             cameraResolution=1024,
             nautoNetName = 'pfig_nauto.net',
             overSampleFactor=2):
        '''
        can call with nautoNetName = 'pfigInv_nauto.net' for inverse pole figures
        '''
        import os
        (fd, scriptFName) = tempfile.mkstemp(suffix='.dx_script', text=True)
        #os.write(fd, 'include "%s/%s"' % (self.path, nautoNetName))
        imageFileName = imageNamePrefix + '.' + imageFormat
        os.write(fd, 'include "%s"' % (nautoNetName)) # trust DXINCLUDE
        os.write(fd, 'main_WriteImage_1_in_2 = "%s";' % (imageNamePrefix))
        os.write(fd, 'main_WriteImage_1_in_3 = "%s";' % (imageFormat))
        os.write(fd, 'main_FileSelector_1_out_1 = "%s";' % (dxInFile))
        os.write(fd, 'main_Integer_2_out_1 = %d ;' % (cameraResolution))
        os.write(fd, 'main_Integer_1_out_1 = %d ;' % (overSampleFactor))
        os.write(fd, 'main_String_1_out_1 = "%s";' % (dataName))
        if greyScale:
            os.write(fd, 'main_Toggle_8_out_1 =  2  ;')
        else:
            os.write(fd, 'main_Toggle_8_out_1 =  1  ;')
        if logThresh is not None:
            os.write(fd, 'main_Toggle_9_out_1 =  2  ;')
            os.write(fd, 'main_Scalar_1_out_1 = %g ;' % (logThresh))
        else:
            os.write(fd, 'main_Toggle_9_out_1 =  1  ;')
        os.write(fd, 'main();')
        os.write(fd, 'quit;')
        os.close(fd)
        command = self.dxExec+' -script '+scriptFName
        print 'running: %s' % (command)
        os.system(command)
        os.remove(scriptFName)
        return imageFileName
    def odf(self, dxInFile, symmGroupString,
            dataName=None,
            cameraResolution = 400,
            imageNamePrefix='odf', imageFormat='tiff',
            minVal=None,
            maxVal=None,
            ):
        import os
        nautoNetName = 'odf_nauto.net'
        (fd, scriptFName) = tempfile.mkstemp(suffix='.dx_script', text=True)
        #os.write(fd, 'include "%s/%s"' % (self.path, nautoNetName))
        imageFileNames = []
        imageFileNames.append(imageNamePrefix + '_surface.' + imageFormat)
        imageFileNames.append(imageNamePrefix + '_slice.' + imageFormat)
        os.write(fd, 'include "%s"' % (nautoNetName)) # trust DXINCLUDE
        os.write(fd, 'main_String_1_out_1 = "%s";' % (imageNamePrefix))
        if dataName is not None:
            os.write(fd, 'main_String_2_out_1 = "%s";' % (dataName))
        os.write(fd, 'main_WriteImage_1_in_3 = "%s";' % (imageFormat))
        os.write(fd, 'main_WriteImage_2_in_3 = "%s";' % (imageFormat))
        os.write(fd, 'main_FileSelector_1_out_1 = "%s";' % (dxInFile))
        os.write(fd, 'main_MyAutoCamera_1_in_4 = %d ;' % (cameraResolution))
        os.write(fd, 'main_MyAutoCamera_2_in_4 = %d ;' % (cameraResolution))
        if not minVal is None:
            os.write(fd, 'main_AutoColor_1_in_7 = %g ;' % (float(minVal)))
        if not maxVal is None:
            os.write(fd, 'main_AutoColor_1_in_8 = %g ;' % (float(maxVal)))
        os.write(fd, 'main();')
        os.write(fd, 'quit;')
        os.close(fd)
        command = self.dxExec+' -script '+scriptFName
        print 'running: %s' % (command)
        os.system(command)
        os.remove(scriptFName)
        return imageFileNames

def plotPfigs(vals, mPfig, plotList=None, dataFilePrefix=None, names=None, **args):
    'vals can be list of values or [nPoleFigures, nMeshNodes] array'
    from femODFUtil import dxWrap
    import sys

    nPF = len(vals)
    
    if dataFilePrefix is None:
        dxFS = dxWrap.DxFS() #  prefix='testPfig', overWrite=True
    else:
        dxFS = dxWrap.DxFS(prefix=dataFilePrefix, overWrite=True)
    print >> sys.stdout, 'handing mesh to dxFS'
    dxFS.setMesh(mPfig)
    print >> sys.stdout, 'attaching nodal data'
    pList = []
    for iPF in range(nPF):
        v = vals[iPF][:]
        if names is not None:
            dataName = 'p_'+names[iPF]
        else:
            dataName = 'p_%03d' % (iPF)
        pList.append(dataName)
        dxFS.attachNodalData(v, dataName, asData=False)
    dxFS.close()
    #
    dxCL = dxWrap.DxCL()
    imageFileNames = []
    for iPF in range(nPF):
        if names is not None:
            imageNamePrefix = names[iPF]
        else:
            imageNamePrefix = 'pfig-%03d_of_%03d' % (iPF,nPF)
        imageFileName = dxCL.pfig(dxFS.getFName(), imageNamePrefix=imageNamePrefix, 
                                  dataName=pList[iPF],
                                  **args)
        imageFileNames.append(imageFileName)
        if plotList is not None:
            plotList[iPF](imageFileName)
    
    return imageFileNames

def plotODF(v, mRodr, dataFilePrefix=None, name=None, **args):
    'all-in-one utility function'
    from femODFUtil import dxWrap
    import sys
    #
    if dataFilePrefix is None:
        dxFS = dxWrap.DxFS() #  prefix='testPfig', overWrite=True
    else:
        dxFS = dxWrap.DxFS(prefix=dataFilePrefix, overWrite=True)
    print >> sys.stdout, 'handing mesh to dxFS'
    dxFS.setMesh(mRodr)
    print >> sys.stdout, 'attaching nodal data'
    dxFS.attachNodalData(v, 'A')
    dxFS.close()
    #
    # ... want to code edges of fundamental region and mechanism for putting in DX file?
    #
    dxCL = dxWrap.DxCL()
    symmGroupString = mRodr.getSymmGroup()
    prefix = name
    if prefix is None: prefix = 'odf'
    imageFileName = dxCL.odf(dxFS.getFName(), symmGroupString, imageNamePrefix=prefix, **args)
    return imageFileName




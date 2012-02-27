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
    import femODF.FemMesh
    import femODF.FemRodrigues
    import femODF.ElemType
    import sidl.RuntimeException
except:
    print >> sys.stderr, 'error doing imports; try sourcing /usr/apps/dlsmm/bin/setup.sh or the like'
    sys.exit(1)

haveSciPy = True
try:
    from scipy import sparse
except:
    haveSciPy = False

def convertFemODFCoo(A):
    'convert to format that can be used in python for operations'
    if haveSciPy:
        try:
            (AI1, AJ1, AV) = A.getIJV()
        except sidl.RuntimeException._Exception:
            import sys
            (etype, eobj, etb) = sys.exc_info()
            # eobj is the SIDL exception object
            print >> sys.stderr, "Exception note: %s" % (eobj.exception.getNote())
            print >> sys.stderr, "Exception trace: %s" % (eobj.exception.getTrace())
            print >> sys.stderr, "About to call sys.exit due to exception"
            sys.exit(1)
        AI = AI1 - 1 # convert to number from 0
        AJ = AJ1 - 1 # convert to number from 0
        mShape = tuple(A.getShape())
        A_conv = sparse.coo_matrix((AV,(AI,AJ)),shape=mShape).tocsr()
    else:
        A_conv = A.getDense()
    return A_conv

def makeRodrMesh(conn, coords, elemType, quadRule, symmGroupString):
    try:    
        mRodr = femODF.FemRodrigues.FemRodrigues()
        #mBase = femODF.FemMesh.FemMesh(mRodr) # cast; not necessary
        #mBase.init(elemType, quadRule, meshes[meshNum]['conn'], meshes[meshNum]['coords'])
        mRodr.init(elemType, quadRule, conn, coords)
        
        mRodr.setupRodr(symmGroupString)
        mRodr.correctBoundary() # only because got coords from file with some possible loss of symmetry; do not anticipate needing this if get from unpack() method call
    except sidl.RuntimeException._Exception:
        (etype, eobj, etb) = sys.exc_info()
        # eobj is the SIDL exception object
        print >> sys.stderr, "Exception note: %s" % (eobj.exception.getNote())
        print >> sys.stderr, "Exception trace: %s" % (eobj.exception.getTrace())
        print >> sys.stderr, "About to call sys.exit due to exception"
        sys.exit(1)

    return mRodr

def makeRodrMeshBall(conn, coords, symmGroupString, ballSize, cRF, 
                     forceSlave=None,
                     elemType=femODF.ElemType.stdE4, 
                     quadRule=[32015, 22004]):
    import arrayUtil
    quadRule = arrayUtil.toArray(quadRule)
    try:    
       mRodr = femODF.FemRodrigues.FemRodrigues()
       #mBase = femODF.FemMesh.FemMesh(mRodr) # cast; not necessary
       #mBase.init(elemType, quadRule, meshes[meshNum]['conn'], meshes[meshNum]['coords'])
       mRodr.init(elemType, quadRule, conn, coords)

       if forceSlave is None:
           mRodr.setupRodrBall(symmGroupString, ballSize, cRF)
       else:
           mRodr.setupRodrBallWithSlaved(symmGroupString, ballSize, cRF, forceSlave)
    except sidl.RuntimeException._Exception:
       (etype, eobj, etb) = sys.exc_info()
       # eobj is the SIDL exception object
       print >> sys.stderr, "Exception note: %s" % (eobj.exception.getNote())
       print >> sys.stderr, "Exception trace: %s" % (eobj.exception.getTrace())
       print >> sys.stderr, "About to call sys.exit due to exception"
       sys.exit(1)

    return mRodr

def getStorageDir(storageDir=None):
    import femODFUtil
    import stat, os

    if storageDir is None:
        femODFPath0 = femODFUtil.getPath()[0]
        # mode = stat.S_IMODE(os.stat(storageDir)[stat.ST_MODE])
        # mode & getattr(stat,"S_IWUSR")
        # mode & getattr(stat,"S_IWGRP")
        # mode & getattr(stat,"S_IWOTH")
        if os.access(femODFPath0, os.W_OK):
            storageDir = femODFPath0
        elif os.access(os.curdir, os.W_OK):
            subDirName = 'femODFUtil_storage'
            storageDir = os.path.join(os.curdir, subDirName)
        else:
            raise RuntimeError, 'do not know where to put storage directory'
    
    if os.path.exists(storageDir):
        if not stat.S_ISDIR(os.stat(storageDir)[stat.ST_MODE]):
            raise RuntimeError, 'exists but is not a directory : '+str(storageDir)
    else:
        os.mkdir(storageDir)

    return storageDir

def refineTetMesh(*args):
    '''
    Note: works with conn numbered from 1 (converts internally to 0-based numbering)
    '''
    import femODF.ElemType
    from scipy import sparse

    if len(args) == 1:
        mesh = args[0]
        elemType = mesh.getElemType()
        assert elemType == femODF.ElemType.stdE4, 'wrong elemnt type'
        conn1     = mesh.getConn()
        coords    = mesh.getCoords()
    elif len(args) == 2:
        conn1   = args[0]
        coords  = args[1]
    else:
        raise RuntimeError, 'wrong number of arguments'
    conn      = conn1-1

    nnpe = 4
    nd = 3
    assert conn.shape[0] == nnpe, 'nnpe mismatch'
    assert coords.shape[0] == nd, 'nd mismatch'

    nelemMult = 8
    npMap = num.array([
            2, 3, 4, 8, #
            6, 4, 5, 9, #
            1, 2, 6, 7,  # 
            7, 8, 9, 10, #
            6, 4, 9, 8, # 6, 9, 8, 4, 
            4, 6, 2, 8, # 6, 8, 2, 4, 
            7, 9, 8, 6, # 8, 9, 6, 7, 
            8, 2, 7, 6  # 8, 6, 2, 7
            ]).reshape(4,nelemMult,order='f') - 1
    
    nNodesOld = coords.shape[1]
    nElemsOld = conn.shape[1]
    midNodeConn = sparse.dok_matrix((nNodesOld, nNodesOld), dtype='int')

    connNew = getMem([nnpe,nElemsOld*nelemMult], typeInt=True)
    # num.zeros((4,nElemsOld*nelemMult),dtype='int')

    edgeConnList = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]

    iNodeNew = nNodesOld
    for iElem in range(nElemsOld):
        for ni,nj in edgeConnList:
            if midNodeConn[conn[ni,iElem], conn[nj,iElem]] == 0:
                midNodeConn[conn[ni,iElem], conn[nj,iElem]] = iNodeNew
                midNodeConn[conn[nj,iElem], conn[ni,iElem]] = iNodeNew
                iNodeNew += 1
        connE4  = conn[:,iElem]
        connE10 = num.array([
            connE4[0], 
            midNodeConn[connE4[0], connE4[1]],
            connE4[1],
            midNodeConn[connE4[1], connE4[2]],
            connE4[2],
            midNodeConn[connE4[0], connE4[2]],
            midNodeConn[connE4[0], connE4[3]],
            midNodeConn[connE4[1], connE4[3]],
            midNodeConn[connE4[2], connE4[3]],
            connE4[3]
            ])
        connNew[:,iElem*nelemMult:(iElem+1)*nelemMult] = connE10[npMap]

    nNodesNew = iNodeNew
    nNewNodes = nNodesNew - nNodesOld
    assert midNodeConn.getnnz() == nNewNodes*2, 'mismatch in new number of nodes'
    
    coordsNew = getMem([nd,nNodesNew])
    # num.zeros((3,nNodesNew),dtype='float64')
    coordsNew[:,0:nNodesOld] = coords[:,:]
    for ij, val in midNodeConn.iteritems():
        #print 'ij %s val %s' % (str(ij),str(val))
        coordsNew[:,val] = 0.5 * (coords[:,ij[0]] + coords[:,ij[1]])
    
    connNew1 = connNew+1
    return (connNew1, coordsNew)

def setBallRMat(mRodr, rMat):
    '''
    set cRF for a ball mesh to corresponde to a lattice orientation
    
    see mesh_ops_odf_mod.F90:form_poleProj_system_ball and rot_mod.F90:c_rf_rot_refvec
    '''
    import arrayUtil
    cRF = arrayUtil.getMem((3,3))
    cRF[:,:] = rMat.T[:,:]
    mRodr.resetBallCRF(cRF)
    return

def threeSliceODF(mRodr, x, nP=401, prefix=None, cmap=None):
    import plotWrap
    from scipy.interpolate import LinearNDInterpolator
    from matplotlib import cm
    # from detector import getCMap
    
    cmap = cmap or cm.jet
    
    rCoords = mRodr.getCoords()
    xF = mRodr.toFull(x)
    # sliceSize = rBallSizePad
    sliceSize = num.abs(rCoords).max()*1.05
    ntrpSoln = LinearNDInterpolator(rCoords.T, xF)
    
    vmax = xF.max()
    
    # from matplotlib import colors
    # norm = colors.Normalize(vmin=0, vmax=vmax)
    
    if prefix is not None:
        titlePrefix = '(%s) ' % (prefix)
    else:
        titlePrefix = ''
    
    grid1D = num.linspace(-sliceSize, sliceSize, nP)
    gridX, gridY = num.meshgrid(grid1D, grid1D)
    
    pwList = []
    
    pwSx = plotWrap.PlotWrap(title=titlePrefix+'x slice')
    #
    vals = ntrpSoln( num.zeros(gridX.shape), gridX, gridY )
    # clrs = cmap(norm(vals))[:,:,0:3]
    im = pwSx( gridX, gridY, vals, cmap=cmap, vmax=vmax, interpolation='bilinear', aspect='equal')
    pwSx.colorbar(thing=im)
    if prefix is not None:
        pwSx.save(filename=prefix+'_ODFslice-sx.png')
    pwList.append(pwSx)
    
    pwSy = plotWrap.PlotWrap(title=titlePrefix+'y slice')
    #
    vals = ntrpSoln( gridY, num.zeros(gridX.shape), gridX )
    #clrs = cmap(norm(vals))[:,:,0:3]
    im = pwSy( gridX, gridY, vals, cmap=cmap, vmax=vmax, interpolation='bilinear', aspect='equal')
    pwSy.colorbar(thing=im)
    if prefix is not None:
        pwSy.save(filename=prefix+'_ODFslice-sy.png')
    pwList.append(pwSy)
    
    pwSz = plotWrap.PlotWrap(title=titlePrefix+'z slice')
    #
    vals = ntrpSoln( gridX, gridY, num.zeros(gridX.shape) )
    # clrs = cmap(norm(vals))[:,:,0:3]
    im = pwSz( gridX, gridY, vals, cmap=cmap, vmax=vmax, interpolation='bilinear', aspect='equal')
    pwSz.colorbar(thing=im)
    if prefix is not None:
        pwSz.save(filename=prefix+'_ODFslice-sz.png')
    pwList.append(pwSz)
    
    return pwList, ntrpSoln

def multiSliceODF(mRodr, x, nP=401, 
                  nSlices=5, sliceNormal=[0.,0.,1.], 
                  prefix=None, cmap=None):
    import plotWrap
    from scipy.interpolate import LinearNDInterpolator
    from matplotlib import cm
    import orientations as ors
    # from detector import getCMap
    
    sliceNormal = num.asarray(sliceNormal)
    cmap = cmap or cm.jet
    
    rCoords = mRodr.getCoords()
    xF = mRodr.toFull(x)
    # sliceSize = rBallSizePad
    sliceSize = num.abs(rCoords).max()*1.05
    ntrpSoln = LinearNDInterpolator(rCoords.T, xF)
    
    vmax = xF.max()
    
    # from matplotlib import colors
    # norm = colors.Normalize(vmin=0, vmax=vmax)
    
    if prefix is not None:
        titlePrefix = '(%s) ' % (prefix)
    else:
        titlePrefix = ''
    
    grid1D = num.linspace(-sliceSize, sliceSize, nP)
    gridX, gridY = num.meshgrid(grid1D, grid1D)
    gridZ = num.zeros(gridX.shape)
    
    'rotate coordinates'
    rMat = ors.RotInv('align', num.array([0.,0.,1.]), sliceNormal).toMatrix()
    crds = num.vstack((gridX.flatten(), gridY.flatten(), gridZ.flatten()))
    rotCrds = num.dot(rMat, crds)
    gridX = rotCrds[0,:].reshape((nP,nP))
    gridY = rotCrds[1,:].reshape((nP,nP))
    gridZ = rotCrds[2,:].reshape((nP,nP))    

    pwList = []
    
    stepsCrd = num.linspace(-sliceSize, sliceSize, nSlices+2)[1:-1]
    
    for iStep, stepCrd in enumerate(stepsCrd):
        pw = plotWrap.PlotWrap(title=titlePrefix+'slice %s'%(iStep))
        offset = stepCrd * sliceNormal
        gridXThis = gridX + offset[0]
        gridYThis = gridY + offset[1]
        gridZThis = gridZ + offset[2]
        vals = ntrpSoln( gridXThis, gridYThis, gridZThis )
        im = pw( gridX, gridY, vals, cmap=cmap, vmax=vmax, interpolation='bilinear', aspect='equal')
        pw.a.set_xlim(-sliceSize, sliceSize)
        pw.a.set_ylim(-sliceSize, sliceSize)
        pw.colorbar(thing=im)
        pw.show()
        if prefix is not None:
            pw.save(filename=prefix+'_ODFslice-%s.png'%(iStep))
        pwList.append(pw)
    
    return pwList, ntrpSoln

def mFromCC(conn, crds, from1 = True):
    """
    get a mesh from connectivity conn and coordinates crds;
    if from1 then incoming conn is numbered from 1, otherwise from 0;
    the result is always numbered from 1
    """
    import arrayUtil
    
    if isinstance(crds, str):
        rs   = num.loadtxt(crds)
        rsA  = arrayUtil.toArray(rs.T)
    else:
        rsA = arrayUtil.toArray(crds)
    
    if isinstance(conn, str):
        conn1 = num.loadtxt(conn, dtype=arrayUtil.dtypeI)
        connA = arrayUtil.toArray(conn1.T)
    else:
        connA = arrayUtil.toArray(conn)
    if not from1:
        connA = connA + 1
    
    'flip elements right-side-out, otherwise volumes are negative (differing connectivity conventions)'
    'must do on an element-by-element basis -- starting mesh is not uniform in this regard'
    # import copy
    # tmp = copy.copy(connA[2,:])
    # connA[2,:] = connA[1,:]
    # connA[1,:] = tmp[:]
    
    def _calcVol(iEl0):
        p1 = rsA[:,connA[0,iEl0]-1]
        p2 = rsA[:,connA[1,iEl0]-1]
        p3 = rsA[:,connA[2,iEl0]-1]
        p4 = rsA[:,connA[3,iEl0]-1]
        v1 = p2 - p1
        v2 = p3 - p1
        v3 = p4 - p1
        'flipped connectivity:' 
        # v1 = p3 - p1
        # v2 = p2 - p1
        # v3 = p4 - p1
        v12 = num.cross(v1,v2)
        vol = num.dot(v12, v3)
        return vol
    # print >> sys.stdout, 'volume of first element is '+str(vol)
    
    for iEl0 in range(len(connA.T)):
        vol = _calcVol(iEl0)
        if vol < 0:
            # flip element
            tmp = connA[1,iEl0]
            connA[1,iEl0] = connA[2,iEl0]
            connA[2,iEl0] = tmp
    
    m = {
        'ne' : connA.shape[1],
        'nn' : rsA.shape[1],
        'conn' : connA,
        'coords' : rsA,
        }
    
    return m

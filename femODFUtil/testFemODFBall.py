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
'''
. /usr/gapps/mdef/dev/chaos_4_x86_64_ib-bbabel-1.4.0-r6691.91-c0.2.7/dbg/bin/setup.sh
'''
import sys, os

import femODFUtil
import arrayUtil
from arrayUtil import num
from femODFUtil import mRodrUtil, mRodrBallStorage

from mRodrUtil import convertFemODFCoo

testOpH1 = True

testScalings = True

calcCondition = False

######################################################################

import matplotlib
from matplotlib import cm
# import pylab as p
# import mpl_toolkits.mplot3d.axes3d as p3

# fig = p.figure()
# ax = p3.Axes3D(fig)
import plotWrap

######################################################################

cRF = arrayUtil.getMem((3,3))
cRF[:,:] = num.eye(3)

ballSize = 0.1
rBallSize = num.tan(ballSize*0.5)
rBallSizePad = rBallSize*1.05

#mRodr = mRodrBallStorage.getMesh('cubic', ballSize, '10x', cRF)
mRodr = mRodrBallStorage.getMesh('cubic', ballSize, 3, cRF)

vProj = mRodr.getVolProj()
nnr = mRodr.getNNodesReduced()
print >> sys.stdout, 'nnr is %d' % (nnr)
vProjSum = vProj.sum()
print >> sys.stdout, 'vProjSum is %g' % (vProjSum)

##############################

irnp1_master = mRodr.getMasterRNP()
if irnp1_master < 0:
    raise RuntimeError, 'no master node found'
irnp_master = irnp1_master - 1

##############################

(inElem, elemXi) = mRodr.find(num.array([0.01,0.01,-0.02]))
print >> sys.stdout, 'found in element %d' % (inElem)

##############################

opM  = convertFemODFCoo(mRodr.getM()).tocsr()
opMSparsity = float(opM.nnz) / float(num.product(opM.shape))
print >> sys.stdout, 'non-zero fill of opM is %5.2f%%' % (opMSparsity*100)
#
if testOpH1:
    opH1 = convertFemODFCoo(mRodr.getH1()).tocsr()
    opH1Sparsity = float(opH1.nnz) / float(num.product(opH1.shape))
    print >> sys.stdout, 'non-zero fill of opH1 is %5.2f%%' % (opH1Sparsity*100)

##############################

from numpy import random
num.random.seed(seed=0)

import math

nSampleVectors = 12
#
sampleNVectors = arrayUtil.getMem([3,nSampleVectors])
'need sample vectors in a specific area to have them have contributions!'
# sampleNVectors[:,:] = random.rand(3*nSampleVectors).reshape(3,nSampleVectors)
sampleNVectors[0,:] = num.ones(nSampleVectors)
sampleNVectors[1,:] = random.randn(nSampleVectors)*ballSize*0.125 # 0.05
sampleNVectors[2,:] = random.randn(nSampleVectors)*ballSize*0.125 # 0.05

for iSV in range(nSampleVectors):
  mag = math.sqrt((sampleNVectors[:,iSV]*sampleNVectors[:,iSV]).sum())
  sampleNVectors[:,iSV] = sampleNVectors[:,iSV] / mag

crystalVector = arrayUtil.getMem([3])
crystalVector[:] = [1.,0.,0.]
odfPfProj = convertFemODFCoo(mRodr.getPoleProjection(crystalVector, sampleNVectors, True))
odfPfProjSparsity = float(odfPfProj.nnz) / float(num.product(odfPfProj.shape))
print >> sys.stdout, 'non-zero fill of odfPfProj is %5.2f%%' % (odfPfProjSparsity*100)
#dok = odfPfProj.todok()
#dok.values()
'to look at values for parts of paths falling outside of the ball:'
print >> sys.stdout, odfPfProj[:,irnp_master].todok().values()

'these sums should all be 1'
print >> sys.stdout, odfPfProj.sum(axis=1)

rsA = mRodr.getCoords()
conn1 = mRodr.getConn()
conn = conn1 - 1
nnf = rsA.shape[1]

sConn1 = mRodr.getSConn()
sConn = sConn1 - 1
nnps = sConn.shape[0]

p1 = odfPfProj.getrow(0)
colInd = p1.nonzero()[1]
p1Dense = num.array(p1.todense()).reshape(nnr)

p1FullB = mRodr.toFull(p1Dense)
p1TouchedB = p1FullB > 0
#
p1Dense[irnp_master] = 0
p1Full = mRodr.toFull(p1Dense)
p1Touched = p1Full > 0

p1FullNZ = p1Full[p1Touched]

rsA_p1 = rsA[:,p1Full > 0]

verts = ( [ 
        num.array( [ rsA[:,triElem[0]] , rsA[:,triElem[1]] , rsA[:,triElem[2]] ]  )
        for triElem in sConn.T  ] )
x = num.ones(nnf)
# average from nodes to elements to get colors for faces
nnps_inv = 1. / float(nnps)
xSElem = num.array( [
        num.sum(x[triElem]) * nnps_inv
        for triElem in sConn.T ] )

#win = plotWrap.PlotWinP(as3D=True)
#pw = plotWrap.PlotWrap(window=win, axes=0)
pw = plotWrap.PlotWrap(as3D=True, title='pole path nodes')
fig = pw.f
ax = pw.a

import mpl_toolkits.mplot3d.art3d  as ar3
triCol = ar3.Poly3DCollection( verts, cmap=cm.jet , alpha=0.05) # color='blue'
triCol.set_array( xSElem )
#a
#triCol.set_edgecolor(None)
#triCol.set_linewidth( 2.0 )
triCol.set_edgecolors([0., 0., 0., 0.])
triCol.set_linewidth( 0.0 )

p1ElTouched = num.all(p1TouchedB[conn], axis=0)
lines = ( [ [
        num.array( [ rsA[:,tetElem[0]] , rsA[:,tetElem[1]] ]  ),
        num.array( [ rsA[:,tetElem[0]] , rsA[:,tetElem[2]] ]  ),
        num.array( [ rsA[:,tetElem[0]] , rsA[:,tetElem[3]] ]  ),
        num.array( [ rsA[:,tetElem[1]] , rsA[:,tetElem[2]] ]  ),
        num.array( [ rsA[:,tetElem[1]] , rsA[:,tetElem[3]] ]  ),
        num.array( [ rsA[:,tetElem[2]] , rsA[:,tetElem[3]] ]  ),
        ] for tetElem in conn[:,p1ElTouched].T  ] )
lines = reduce(lambda x,y:x+y, lines) # flatten one level
linCol = ar3.Line3DCollection( lines )

ax.add_collection( triCol )
ax.add_collection( linCol )

# ax.plot(rsA_p1[0,:], rsA_p1[1,:], rsA_p1[2,:], 'ro', alpha=1.0)
# ax.scatter(rsA_p1[0,:], rsA_p1[1,:], rsA_p1[2,:], c=p1FullNZ, s=40)
pSize = p1FullNZ * (40/p1FullNZ.mean())
ax.scatter(rsA_p1[0,:], rsA_p1[1,:], rsA_p1[2,:], c=pSize, s=pSize)
ax.set_autoscale_on(True)
# ax.clear()

# plotting in the reference ball, so set limits to +-1
ax.set_xlim3d(-rBallSizePad, rBallSizePad)
ax.set_ylim3d(-rBallSizePad, rBallSizePad)
ax.set_zlim3d(-rBallSizePad, rBallSizePad)
ax.set_aspect('equal')
ax.mouse_init()
# plt.ioff()
# plt.show()
pw.show()

ax.axis('off')
ax.grid(on=False)
pw.show()

# raise RuntimeError, 'stop on purpose'

##############################

# make a pole figure patch around the +z axis
patchSize = ballSize*1.2
dxP = patchSize*2/40.
nP = 41
xP = num.linspace(-patchSize, patchSize, nP)
yP = num.linspace(-patchSize, patchSize, nP)
xGrid, yGrid = num.meshgrid(xP, yP)
zGrid = num.ones(xGrid.shape)
nSampleVectors = xGrid.size # per patch!
sampleNVectors = num.vstack((xGrid.flatten(), yGrid.flatten(), zGrid.flatten()))
for iSV in range(nSampleVectors):
  mag = math.sqrt((sampleNVectors[:,iSV]*sampleNVectors[:,iSV]).sum())
  sampleNVectors[:,iSV] = sampleNVectors[:,iSV] / mag

pfigs = [
    {'crystalVector' : arrayUtil.toArray([1., 0., 0.])},
    {'crystalVector' : arrayUtil.toArray([1., 1., 0.])},
    {'crystalVector' : arrayUtil.toArray([1., 2., 3.])},
    ]
if testScalings:
    pfigArtificialScalings = [1.0, 1.2, 3.0]
else:
    pfigArtificialScalings = [1.0, 1.0, 1.0]
'make polePoints'
import orientations as ors
symmGroupString = mRodr.getSymmGroup()
symm = ors.SymmGroup(symmGroupString)
qArray = symm.getQArray()
for pfig in pfigs:
    cV = pfig['crystalVector']
    from XRD import Symmetry
    cvecs = Symmetry.applySym(cV.reshape(3,1), qArray.T, cullPM=True)
    theseNVecs = []
    for cvec in cvecs.T:
        rMat = ors.RotInv('align', num.array([1.,0.,0.]), cvec).toMatrix()
        theseNVecs.append(num.dot(rMat, sampleNVectors))
    pfig['polePoints'] = num.hstack(theseNVecs)

'the following works for 100 only:'
# sNVx = num.vstack((sampleNVectors[2,:],sampleNVectors[1,:],sampleNVectors[0,:]))
# sNVy = num.vstack((sampleNVectors[0,:],sampleNVectors[2,:],sampleNVectors[1,:]))
# sNVz = sampleNVectors
# sNVAll = num.hstack((sNVx, sNVy, sNVz))

patchSize1 = patchSize + dxP*0.5
xyP1 = num.arange(-patchSize1,patchSize1+dxP*0.5,dxP)
xGrid1, yGrid1 = num.meshgrid(xyP1, xyP1)

from femODFUtil import pfig
# hFactor=1e-3 is too little?
h1Factor=1e-2
# h1Factor = 1e-1 # is okay for under-detemined system
# h1Factor = 1e0 # is too much?
sysBlocks = pfig.calcPfigProjOps(pfigs, mRodr, h1Factor=h1Factor)

if calcCondition:
    from scipy.sparse import linalg
    projOps    = sysBlocks['projOps']
    projOpFull = (sparse.vstack(projOps)).tocsr()
    projOp = projOpFull[:,:-1]
    u, s, vh = linalg.svds(projOp, k=projOp.shape[1]-1)
    conditionNumber = s[-1]/s[0]
    print 'conditionNumber is : %g' % (conditionNumber)


'''
make ODF using normal distribution, Euclidean in Rodrigues space -- 
not taking account of the metric;
this should be fine because this is just a test problem and the metric is pretty
uniform over such a small extent anyway
'''
# get the real coords
rCoords = mRodr.getCoords()
# num.max(rCoords) = num.tan(ballSize/2.)

def distrib(x, width, scale):
    mags = num.apply_along_axis(num.linalg.norm, 0, x)
    retval = scale * num.exp(-(1./(width*width)*mags*mags))
    return retval
rCoordsB = rCoords - num.tile(num.array([num.tan(0.25*ballSize*0.5), 0., 0.]), (nnf,1)).T
aValF = distrib(rCoords, ballSize/6., 1.) + distrib(rCoordsB, ballSize/6., 0.75)
#
aValR = mRodr.toReduced(aValF, True)
aValR[irnp_master] = 0. # force the master boundary node to have zero value

sysBlocks['scaleFactors'] = pfigArtificialScalings
rhsBlocks = pfig.calcPfigProj(sysBlocks, aValR)

pfigData100 = rhsBlocks[0]

pwPatch = plotWrap.PlotWrap(title='pole figure patch')
#iPF = 0
iPF = 1
#iPF = 2
data = pfigData100[nSampleVectors*iPF:nSampleVectors*(iPF+1)].reshape(nP,nP)
cont = pwPatch(xGrid1, yGrid1, data, interpolation='nearest', aspect='equal')
pwPatch.colorbar(thing=cont)

'pretend that we do not know anything about the scale factors'
sysBlocks['scaleFactors'] = num.ones(len(sysBlocks['projOps']))

if testScalings:
    x, scaleFactors = pfig.invert(sysBlocks, rhsBlocks, verbose=True, scaleIntensities=True, solverTol=1e-9)
    print 'got scaleFactors : '+str(scaleFactors)
    sysBlocks['scaleFactors'] = scaleFactors
else:
    x = pfig.invert(sysBlocks, rhsBlocks, verbose=True)

relDiff = num.linalg.norm(aValR - x)/num.linalg.norm(aValR)
print >> sys.stdout, 'relative x difference : '+str(relDiff)

print 'min, max, and mean of solution : %g %g %g' % (x.min(), x.max(), x.mean())
print 'min, max, and mean of true ODF : %g %g %g' % (aValR.min(), aValR.max(), aValR.mean())

pfigDataX = pfig.calcPfigProj(sysBlocks, x)
pfigData100X = pfigDataX[0]

pwPatchX = plotWrap.PlotWrap(title='reconstructed')
data = pfigData100X[nSampleVectors*iPF:nSampleVectors*(iPF+1)].reshape(nP,nP)
cont = pwPatchX(xGrid1, yGrid1, data, interpolation='nearest', aspect='equal')
pwPatchX.colorbar(thing=cont)

import matplotlib.collections as collections

# Note: projOpsAM is done _without_ scale factors
mPfig, nVectors = pfig.getAutoMesh(0.06)
projOpsAM  = pfig.calcPfigProjOps(pfigs, mRodr, polePoints=nVectors)
projValsAM = pfig.calcPfigProj(projOpsAM, x)

'''
pole figure mesh is quads with irregular connectivity;
have not yet spotted a nice mesh plotter for irregulat connectivity;
matplotlib.collections.QuadMesh is for regular connectivity -- maybe not!;
just do interpolate onto a regular grid for now and plot that -- not as pretty, especially around the edges, but it does the trick;
'''

import pfigUtil
eaProj = pfigUtil.n2eap(nVectors, flip=False)
from scipy.interpolate import LinearNDInterpolator
ntrp = LinearNDInterpolator(eaProj.T, projValsAM[0])
xi = num.linspace(-1.05, 1.05, 401)
yi = num.linspace(-1.05, 1.05, 401)
xiGrid, yiGrid = num.meshgrid(xi,yi)
ziGrid = ntrp(xiGrid, yiGrid)
pwAMX = plotWrap.PlotWrap(title='coarse global pole figure')
pwAMX(xiGrid, yiGrid, ziGrid, interpolation='bilinear', aspect='equal')

##############################

"""
# the following sort of works, but matplotlib does not figure out the intersections of the slices -- it just draws the most recent thing on top, so the 3D nature is ruined; do not yet know whether thers is a work-aroudn
from scipy.interpolate import LinearNDInterpolator
xF = mRodr.toFull(x)
rBallSize = num.tan(ballSize*0.5)
sliceSize = rBallSize*1.05
ntrpSoln = LinearNDInterpolator(rCoords.T, xF, fill_value=0) # fill_value can go away?
nP = 101

from matplotlib import colors
norm = colors.Normalize(vmin=0, vmax=xF.max())

pwSlice = plotWrap.PlotWrap(as3D=True, title='ODF solution')
ax = pwSlice.a

xiSoln = num.linspace(-sliceSize, sliceSize, nP)
yiSoln = num.linspace(-sliceSize, sliceSize, nP)
ziSoln = num.linspace(-sliceSize, sliceSize, nP)

xiSolnGrid, yiSolnGrid = num.meshgrid(xiSoln,yiSoln)
ziSolnGrid = num.zeros(xiSolnGrid.shape) # slice in z=0 plane
valsSolnGrid = ntrpSoln(xiSolnGrid, yiSolnGrid, ziSolnGrid) 
# ax.set_autoscale_on(True)
# ax.clear()
# ... may be better to plot using patches if get to the point of being able to slice the mesh
# colors = valsSolnGrid.flatten()
facecolors = cm.jet(norm(valsSolnGrid))[:,:,0:3]
surf = ax.plot_surface(xiSolnGrid, yiSolnGrid, ziSolnGrid, facecolors=facecolors,
                       cstride=1, rstride=1, 
                       linewidth=0, antialiased=False, shade=False)

xiSolnGrid, ziSolnGrid = num.meshgrid(xiSoln,ziSoln)
yiSolnGrid = num.zeros(xiSolnGrid.shape) # slice in y=0 plane
#
valsSolnGrid = ntrpSoln(xiSolnGrid, yiSolnGrid, ziSolnGrid) 
facecolors = cm.jet(norm(valsSolnGrid))[:,:,0:3]
surf = ax.plot_surface(xiSolnGrid, yiSolnGrid, ziSolnGrid, facecolors=facecolors,
                       cstride=1, rstride=1, 
                       linewidth=0, antialiased=False, shade=False)

yiSolnGrid, ziSolnGrid = num.meshgrid(yiSoln,ziSoln)
xiSolnGrid = num.zeros(yiSolnGrid.shape) # slice in x=0 plane
#
valsSolnGrid = ntrpSoln(xiSolnGrid, yiSolnGrid, ziSolnGrid) 
facecolors = cm.jet(norm(valsSolnGrid))[:,:,0:3]
surf = ax.plot_surface(xiSolnGrid, yiSolnGrid, ziSolnGrid, facecolors=facecolors,
                       cstride=1, rstride=1, 
                       linewidth=0, antialiased=False, shade=False)

ax.set_xlim3d(-sliceSize, sliceSize)
ax.set_ylim3d(-sliceSize, sliceSize)
ax.set_zlim3d(-sliceSize, sliceSize)
ax.set_aspect('equal')
ax.mouse_init()
# just example junk:
# ax.w_zaxis.set_major_locator(LinearLocator(10))
# ax.w_zaxis.set_major_formatter(FormatStrFormatter('%.03f'))
# pwSlice.f.colorbar(surf, shrink=0.5, aspect=5)
pwSlice.show()
"""

##############################
# do slices of ODFs
# in three plotWrap instances for backend compatibility reasons (plotWrap deficiency)
from scipy.interpolate import LinearNDInterpolator
xF = mRodr.toFull(x)
sliceSize = rBallSizePad
ntrpSoln = LinearNDInterpolator(rCoords.T, xF)
nP = 401

from matplotlib import colors
vmax = xF.max()
norm = colors.Normalize(vmin=0, vmax=vmax)

grid1D = num.linspace(-sliceSize, sliceSize, nP)
gridX, gridY = num.meshgrid(grid1D, grid1D)

pwSz = plotWrap.PlotWrap(title='z slice (fit)')
#
vals = ntrpSoln( gridX, gridY, num.zeros(gridX.shape) )
# clrs = cm.jet(norm(vals))[:,:,0:3]
im = pwSz( gridX, gridY, vals, cmap=cm.jet, vmax=vmax, interpolation='bilinear', aspect='equal')
pwSz.colorbar(thing=im)
pwSz.save(filename='fit_sz.png')

pwSx = plotWrap.PlotWrap(title='x slice (fit)')
#
vals = ntrpSoln( num.zeros(gridX.shape), gridX, gridY )
# clrs = cm.jet(norm(vals))[:,:,0:3]
im = pwSx( gridX, gridY, vals, cmap=cm.jet, vmax=vmax, interpolation='bilinear', aspect='equal')
pwSx.colorbar(thing=im)
pwSx.save(filename='fit_sx.png')

pwSy = plotWrap.PlotWrap(title='y slice (fit)')
#
vals = ntrpSoln( gridY, num.zeros(gridX.shape), gridX )
#clrs = cm.jet(norm(vals))[:,:,0:3]
im = pwSy( gridX, gridY, vals, cmap=cm.jet, vmax=vmax, interpolation='bilinear', aspect='equal')
pwSy.colorbar(thing=im)
pwSy.save(filename='fit_sy.png')

##############################
# same as above, but replay xF with aValF

from scipy.interpolate import LinearNDInterpolator
rBallSize = num.tan(ballSize*0.5)
sliceSize = rBallSize*1.05
ntrpSoln = LinearNDInterpolator(rCoords.T, aValF)
nPODF = 401

# from matplotlib import colors
vmax = aValF.max()
# norm = colors.Normalize(vmin=0, vmax=vmax)

grid1D = num.linspace(-sliceSize, sliceSize, nPODF)
gridX, gridY = num.meshgrid(grid1D, grid1D)

pwSz = plotWrap.PlotWrap(title='z slice (original)')
#
vals = ntrpSoln( gridX, gridY, num.zeros(gridX.shape) )
# clrs = cm.jet(norm(vals))[:,:,0:3]
im = pwSz( gridX, gridY, vals, cmap=cm.jet, vmax=vmax, interpolation='bilinear', aspect='equal')
pwSz.colorbar(thing=im)
pwSz.save(filename='original_sz.png')

pwSx = plotWrap.PlotWrap(title='x slice (original)')
#
vals = ntrpSoln( num.zeros(gridX.shape), gridX, gridY )
#clrs = cm.jet(norm(vals))[:,:,0:3]
im = pwSx( gridX, gridY, vals, cmap=cm.jet, vmax=vmax, interpolation='bilinear', aspect='equal')
pwSx.colorbar(thing=im)
pwSx.save(filename='original_sx.png')

pwSy = plotWrap.PlotWrap(title='y slice (original)')
#
vals = ntrpSoln( gridY, num.zeros(gridX.shape), gridX )
#clrs = cm.jet(norm(vals))[:,:,0:3]
im = pwSy( gridX, gridY, vals, cmap=cm.jet, vmax=vmax, interpolation='bilinear', aspect='equal')
pwSy.colorbar(thing=im)
pwSy.save(filename='original_sy.png')

'''
... play with h1 factor?

... if go to some mRodrStorage type thing, need to watch for ballSize in mesh storage -- H1 and M operators depend on what ballSize is!

... # plot 3D slices through the mRodr ball

scipy.optimize.slsqp?
	# Minimize a function using Sequential Least SQuares Programming
	# (at least) box and linear equality constraints
	takes function for scalar f and its gradient?
scipy.optimize.lbfgsb?
	# Minimize a function func using the L-BFGS-B algorithm.
	# handles box constraints
scipy.optimize.tnc?
	# Minimize a function with variables subject to bounds, using gradient information.

...



'''

'''
cV001 = toArray([0.,0.,1.])
poleInfo = [('001',cV001)]
pfigsInit = []
for inFix,cV in poleInfo:
    dataFile = os.path.join(dataPath, '8378.%s.data' % inFix)
    data = num.array(fileUtil.readDataFlat(dataFile))
    pfigsInit.append({
            'crystalVector' : cV,
            'polePoints'    : polePoints,
            'poleVals'      : data,
            })
pfigs = pfigsInit

... # for iron data:
(*) initial grain location 
	figure out how to set cRF from orientation
(*) swag of grain spread
(*) spot locations
	including eta and omega ranges
(*) make pfig data around those spot locations, by hkl type
	(*) like in arrayUtilhistoFit?
	(*) modify xrdUtils.CollapseOmeEta to be like that too?
	(*) mask any elements that do not get any hits in radial rebinning?
	(*) cann apply in 3D to bit in 2-theta as well? -- num.histogramdd
	(*) can look over frames and sub-loop over spots, pulling out only what is needed?
	(*) make sure do background subtraction
		(*) clean up background to despecle? -- JVB has code to do this?
	(*) and then apply Lorenz factor type stuff
(*) and zero any negative values?
'''

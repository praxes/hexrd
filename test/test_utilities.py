import sys, os, time, random
import numpy as np

from hexrd.xrd import transforms as xf
from hexrd.xrd import transforms_CAPI as xfcapi

epsf = 2.2e-16

vec = np.array([[random.uniform(-np.pi,np.pi),random.uniform(-np.pi,np.pi),random.uniform(-np.pi,np.pi)]])
vHat1 = xf.unitVector(vec.T)
vHat2 = xfcapi.unitRowVector(vec)
print "unitVector results match:             ",np.linalg.norm(vHat1.T-vHat2)/np.linalg.norm(vHat1) < epsf

tAng = np.array([0.0011546340766314521,-0.0040527538387122993,-0.0026221336905160211])
rMat1 = xf.makeDetectorRotMat(tAng)
rMat2 = xfcapi.makeDetectorRotMat(tAng)
print "makeDetectorRotMat results match:     ",np.linalg.norm(rMat1-rMat2)/np.linalg.norm(rMat1) < epsf

oAng = np.array([-0.0011591608938627839,0.0011546340766314521])
rMat1 = xf.makeOscillRotMat(oAng)
rMat2 = xfcapi.makeOscillRotMat(oAng)
print "makeOscillRotMat results match:       ",np.linalg.norm(rMat1-rMat2)/np.linalg.norm(rMat1) < epsf

eMap = np.array([ 0.66931818,-0.98578066,0.73593251])
rMat1 = xf.makeRotMatOfExpMap(eMap)
rMat2 = xfcapi.makeRotMatOfExpMap(eMap)
print "makeRotMatOfExpMap results match:     ",np.linalg.norm(rMat1-rMat2)/np.linalg.norm(rMat1) < epsf

axis = np.array([ 0.66931818,-0.98578066,0.73593251])
rMat1 = xf.makeBinaryRotMat(axis)
rMat2 = xfcapi.makeBinaryRotMat(axis)
print "makeBinaryRotMat results match:       ",np.linalg.norm(rMat1-rMat2)/np.linalg.norm(rMat1) < epsf

bHat = np.array([0.0,0.0,-1.0])
eta  = np.array([1.0,0.0,0.0])
rMat1 = xf.makeEtaFrameRotMat(bHat,eta)
rMat2 = xfcapi.makeEtaFrameRotMat(bHat,eta)
print "makeEtaFrameRotMat results match:     ",np.linalg.norm(rMat1-rMat2)/np.linalg.norm(rMat1) < epsf

angles = np.array([random.uniform(-np.pi,np.pi),random.uniform(-np.pi,np.pi),random.uniform(-np.pi,np.pi)])
aMin = np.array([random.uniform(-np.pi,np.pi),random.uniform(-np.pi,np.pi)])
aMax = np.array([random.uniform(-np.pi,np.pi),random.uniform(-np.pi,np.pi)])
va1 = xf.validateAngleRanges(angles,aMin,aMax)
va2 = xfcapi.validateAngleRanges(angles,aMin,aMax)

print angles
print aMin,aMax
print va1,va2
print "validateAngleRanges results match:    ",np.array_equiv(va1,va2)

angle = np.array([0.1])
axis = np.array([[ 0.47792,-0.703887,0.525486]])
axis /= np.linalg.norm(axis)
vecs = np.array([[ 1.0, 2.0, 3.0]])
rVec1 = xf.rotate_vecs_about_axis(angle,axis.T,vecs.T)
rVec2 = xfcapi.rotate_vecs_about_axis(angle,axis,vecs)
print "rotate_vecs_about_axis results match: ",np.linalg.norm(rVec1.T-rVec2)/np.linalg.norm(rVec1) < epsf

### Timing Results ###

vec = np.array([[random.uniform(-np.pi,np.pi),random.uniform(-np.pi,np.pi),random.uniform(-np.pi,np.pi)]])
n = 26000
start1 = time.time()                      # time this
for i in range(n):
    vHat1 = xf.unitVector(vec.T)
elapsed1 = (time.time() - start1)
start2 = time.time()                      # time this
for i in range(n):
    vHat2 = xfcapi.unitRowVector(vec)
elapsed2 = (time.time() - start2)
print "Time for %d unitVector:            %g v. %g (%f)"%(n,elapsed1/n,elapsed2/n,elapsed1/elapsed2)

tAng = np.array([0.0011546340766314521,-0.0040527538387122993,-0.0026221336905160211])
n = 13000
start1 = time.time()                      # time this
for i in range(n):
    rMat1 = xf.makeDetectorRotMat(tAng)
elapsed1 = (time.time() - start1)
start2 = time.time()                      # time this
for i in range(n):
    rMat2 = xfcapi.makeDetectorRotMat(tAng)
elapsed2 = (time.time() - start2)
print "Time for %d makeDetectorRotMat:    %g v. %g (%f)"%(n,elapsed1/n,elapsed2/n,elapsed1/elapsed2)

oAng = np.array([-0.0011591608938627839,0.0011546340766314521])
n = 20000
start1 = time.time()                      # time this
for i in range(n):
    rMat1 = xf.makeOscillRotMat(oAng)
elapsed1 = (time.time() - start1)
start2 = time.time()                      # time this
for i in range(n):
    rMat2 = xfcapi.makeOscillRotMat(oAng)
elapsed2 = (time.time() - start2)
print "Time for %d makeOscillRotMat:      %g v. %g (%f)"%(n,elapsed1/n,elapsed2/n,elapsed1/elapsed2)

eMap = np.array([ 0.66931818,-0.98578066,0.73593251])
n = 20000
start1 = time.time()                      # time this
for i in range(n):
    rMat1 = xf.makeRotMatOfExpMap(eMap)
elapsed1 = (time.time() - start1)
start2 = time.time()                      # time this
for i in range(n):
    rMat2 = xfcapi.makeRotMatOfExpMap(eMap)
elapsed2 = (time.time() - start2)
print "Time for %d makeRotMatOfExpMap:    %g v. %g (%f)"%(n,elapsed1/n,elapsed2/n,elapsed1/elapsed2)

axis = np.array([ 0.66931818,-0.98578066,0.73593251])
n = 78000
start1 = time.time()                      # time this
for i in range(n):
    rMat1 = xf.makeBinaryRotMat(axis)
elapsed1 = (time.time() - start1)
start2 = time.time()                      # time this
for i in range(n):
    rMat2 = xfcapi.makeBinaryRotMat(axis)
elapsed2 = (time.time() - start2)
print "Time for %d makeBinaryRotMat:      %g v. %g (%f)"%(n,elapsed1/n,elapsed2/n,elapsed1/elapsed2)

bHat = np.array([0.0,0.0,-1.0])
eta  = np.array([1.0,0.0,0.0])
n = 7500
start1 = time.time()                      # time this
for i in range(n):
    rMat1 = xf.makeEtaFrameRotMat(bHat,eta)
elapsed1 = (time.time() - start1)
start2 = time.time()                      # time this
for i in range(n):
    rMat2 = xfcapi.makeEtaFrameRotMat(bHat,eta)
elapsed2 = (time.time() - start2)
print "Time for %d makeEtaFrameRotMat:     %g v. %g (%f)"%(n,elapsed1/n,elapsed2/n,elapsed1/elapsed2)

angles = np.array([-np.pi/3.0,np.pi/6.0,2.0*np.pi/3.0])
aMin = np.array([-0.75*np.pi,0.25*np.pi])
aMax = np.array([-0.25*np.pi,0.75*np.pi])
n = 3500
start1 = time.time()                      # time this
for i in range(n):
    va1 = xf.validateAngleRanges(angles,aMin,aMax)
elapsed1 = (time.time() - start1)
start2 = time.time()                      # time this
for i in range(n):
    va2 = xfcapi.validateAngleRanges(angles,aMin,aMax)
elapsed2 = (time.time() - start2)
print "Time for %d validateAngleRanges:    %g v. %g (%f)"%(n,elapsed1/n,elapsed2/n,elapsed1/elapsed2)

angle = np.array([0.1])
axis = np.array([[ 0.47792,-0.703887,0.525486]])
axis /= np.linalg.norm(axis)
vecs = np.array([[ 1.0, 2.0, 3.0]])
n = 4500
start1 = time.time()                      # time this
for i in range(n):
    rVec1 = xf.rotate_vecs_about_axis(angle,axis.T,vecs.T)
elapsed1 = (time.time() - start1)
start2 = time.time()                      # time this
for i in range(n):
    rVec2 = xfcapi.rotate_vecs_about_axis(angle,axis,vecs)
elapsed2 = (time.time() - start2)
print "Time for %d rotate_vecs_about_axis: %g v. %g (%f)"%(n,elapsed1/n,elapsed2/n,elapsed1/elapsed2)

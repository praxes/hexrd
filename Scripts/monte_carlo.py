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
import detector
import spotFinder
import grain_wrap_precession as gw
import grain as G
import valUnits
import grain_wrap as gw
import numpy as num
import numpy.random
import scipy.stats
import Rotations
import grain
import detector
import grain_wrap_precession as gw
import matrixUtils as mUtil
from numpy import dot
from scipy.linalg import inv
from Vector_Data_Structures_LLNL import Get_LLNL_Angles,mapOme
from Vector_Data_Structures_LLNL import Get_Pixel_Polar_Coords, alphabeta_polar_to_etatheta, Get_LLNL_Angles,mapOme,Get_LLNL_Pixel_Sols_From_Recip_Vector
from scipy import optimize
from Vector_funcs import *
import uncertainty_analysis
try:
    from XRD import readGE
    from XRD import spotFinder
    from XRD import crystallography
    from XRD.examples.rubyNIST import planeData
except:
    # probably running from XRD directory; put ".." in path just in case
    import sys, os
    sys.path.append(os.path.join(os.getcwd(),os.path.pardir))
    import readGE
    import spotFinder
    import crystallography
    from examples.rubyNIST import planeData,wavelength

def MC_make_angles(R,U, gIs, sigmas,wmin,wmax):
    assert len(gIs) == len(sigmas)
    rIs = []
    F = dot(R,U)
    FinvT = inv(F).T
    angles = []
    for i in range(len(gIs)):
        gI = gIs[i]
        rI = dot(FinvT,gI)
        rI_angle_sol = Get_LLNL_Angles(rI, wavelength)
        ct = 0
        s_tth, s_eta,s_w = sigmas[i]
        r_tth = numpy.random.normal(0, s_tth, 1)
        r_eta = numpy.random.normal(0, s_eta, 1)
        r_w = numpy.random.normal(0, s_w, 1)
        for key in rI_angle_sol.keys():
            ct+=1
            tth,eta,w = rI_angle_sol[key]
            w = mapOme(w,wmin,wmax)
            eta = mapOme(eta,-num.pi,num.pi)
            if w is not None:
                hit = [tth+r_tth,eta+r_eta,w + r_w]
                angles.append(hit)
                if ct==2:
                    print 'double hit', gI, i, rI
    return angles
def MC_fitC_wrap_order(Uparams,Rparams,gIs,angles,wavelength,weighting = False, report = True, evaluate = False):
    return MC_fitC_angle_weights(Rparams,Uparams,gIs,angles,wavelength,weighting,report,evaluate)

def MC_fitC_angle_weights(Rparams, Uparams, gIs,angles0, wavelength, weighting = False, report = True, evaluate = False):
    """
    alternate fitting approach for U, uses no Rotation information.
    Cparams: [C11,C22,C33,C12,C23,C13]
    """
    from Vector_funcs import Mag
    from numpy import dot,sqrt
    import uncertainty_analysis
    'internal objective function'
    def _fitC_residual(angles, r0, invCmat, wavelength):
        rI = gw.makeARecipVector(angles, wavelength)
        #print rI
        df = 1./Mag(rI)
        d0 = 1./Mag(r0)
        gii_0 = dot(r0,r0)
        gii = dot(dot(invCmat,r0),r0)
        r_i = df/d0 - sqrt(gii_0)/(sqrt(gii))
        return r_i
    def _fitC(Cparams,gIs,angles,wavelength):
        Cmat = gw.vec_to_sym(Cparams)
        invCmat = inv(Cmat)
        out = []
        for i in range(len(gIs)):
            gI = gIs[i]
            angCOM = angles[i] #, spot in gI_spots:
            r_i = _fitC_residual(angCOM,gI, invCmat, wavelength)
            if weighting:
                weight = uncertainty_analysis.propagateUncertainty(_fitC_residual, weighting[i], 1e-8,angCOM,gI,invCmat,wavelength)
            else:
                weight = 1
            out.append(r_i/weight)
        return num.array(out)
    
    if Rparams == 'default':
        angl,axxx = Rotations.angleAxisOfRotMat(num.eye(3))
        Rparams = (angl*axxx).flatten()
    else:
        Rparams = Rparams
    if Uparams=='default':
        Uparams = sym_to_vec(num.eye(3))
    else:
        Uparams = Uparams
        
    U = gw.vec_to_sym(Uparams)
    R = Rotations.rotMatOfExpMap(Rparams)
    F = num.dot(R,U)
    C_ = num.dot(F.T,F)
    Cparams = gw.sym_to_vec(C_)

    if evaluate:
        return _fitC(Cparams,gIs,angles0,wavelength)
    
    optResults = optimize.leastsq(_fitC, x0 = Cparams, args = (gIs,angles0,wavelength), full_output = 1)
    Cparams = optResults[0]    
    C = gw.vec_to_sym(Cparams)
    from scipy.linalg import eig
    'spectral decomposition to get U'
    eval,evec = eig(C)
    l1_sqr,l2_sqr,l3_sqr = num.asarray(eval,dtype = 'float')
    u1,u2,u3 = evec[:,0],evec[:,1],evec[:,2]
    U = sqrt(l1_sqr)*num.outer(u1,u1) + sqrt(l2_sqr)*num.outer(u2,u2) + sqrt(l3_sqr)*num.outer(u3,u3)
    if report:
        print "refined stretch matrix: \n%g, %g, %g\n%g, %g, %g\n%g, %g, %g\n" \
              % (U[0, 0], U[0, 1], U[0, 2], \
                 U[1, 0], U[1, 1], U[1, 2], \
                 U[2, 0], U[2, 1], U[2, 2], )
        

    return optResults
def _fitR_residual(angles, r0, Fmat, wavelength):
    Cmat = dot(Fmat.T,Fmat)
    
    rI = gw.makeARecipVector(angles,wavelength)
    n = Unit(rI)
    N = Unit(r0)
    #print Fmat,Cmat
    
    alpha_ = sqrt(Inner_Prod(Star(Cmat),N,N))
    #print alpha_
    r_i = Inner_Prod(Star(Fmat),n,N) - alpha_
    #print 'r_i'
    return r_i

def MC_fitRotation_angles(Rparams, Uparams, gIs,angles0, wavelength, weighting = False, report = True,evaluate = False):
    assert len(gIs)==len(angles0), 'different lengths'
    def _fitR_residual(angles, r0, Fmat, wavelength):
        Cmat = dot(Fmat.T,Fmat)
        
        rI = gw.makeARecipVector(angles,wavelength)
        n = Unit(rI)
        N = Unit(r0)
        #print Fmat,Cmat
        
        alpha_ = sqrt(Inner_Prod(Star(Cmat),N,N))
        #print alpha_
        r_i = Inner_Prod(Star(Fmat),n,N) - alpha_
        #print 'r_i'
        return r_i        
    def _fitRotation_lsq(Rparams, gIs,angles0, U,wavelength):
        out = []
        Rmat = Rotations.rotMatOfExpMap(Rparams)
        Fmat = dot(Rmat,U)
        for i in range(len(gIs)):
            gI,angles = gIs[i], angles0[i]
            r_i = _fitR_residual(angles,gI, Fmat, wavelength)
            if weighting:
                weight = uncertainty_analysis.propagateUncertainty(_fitR_residual, weighting[i], 1e-8,angles,gI,Fmat,wavelength)
                maxiter = 100
                ct = 0
                while weight==0:
                    weight = uncertainty_analysis.propagateUncertainty(_fitR_residual, weighting[i], 1e-8,angles,gI,Fmat,wavelength)
                    ct+=1
                    if ct>=maxiter:
                        print 'zero weight error, i = ',i
                        weight = 1
                        break
                    
            else:
                weight = 1.
            out.append(r_i/weight)
            if weight == 0:
                print 'wefiht0',i,weighting[i],Fmat
            

        out = num.array(out)

    #    print out
        return out
    if Rparams == 'default':
        angl,axxx = Rotations.angleAxisOfRotMat(num.eye(3))
        Rparams = (angl*axxx).flatten()
    else:
        Rparams = Rparams
    if Uparams=='default':
        Uparams = sym_to_vec(num.eye(3))
    else:
        Uparams = Uparams
        
    U = gw.vec_to_sym(Uparams)
    if evaluate:
        return _fitRotation_lsq(Rparams,gIs,angles0,U,wavelength)
    
    optResults = optimize.leastsq(_fitRotation_lsq, x0 = Rparams, args = (gIs,angles0, U, wavelength), full_output = 1)
    r1 = optResults[0]
    U1 = Rotations.rotMatOfExpMap(r1)
    if report:
        print "refined orientation matrix: \n%g, %g, %g\n%g, %g, %g\n%g, %g, %g\n" \
              % (U1[0, 0], U1[0, 1], U1[0, 2], \
                 U1[1, 0], U1[1, 1], U1[1, 2], \
                 U1[2, 0], U1[2, 1], U1[2, 2], )
        
    return optResults

def MC_fitAll_Angles_wrap(allparams,agrain,wmin,wmax,angles, weighting = True):
    pVec = allparams[0:3]
    Rparams = allparams[3:6]
    Uparams = allparams[6:12]
    U = gw.vec_to_sym(Uparams)
    R = Rotations.rotMatOfExpMap(Rparams)
    #angles = MC_make_angles_from_spots_grain(R,U,agrain,pVec,wmin,wmax,tweak = False)
    return MC_Angle_Fit(R,U,agrain,pVec,wmin,wmax,angles,weighting)
def MC_fitU_Angle_wrap(Uparams,R,agrain,pVec,wmin,wmax,angles,weighting = True):
    uMat = gw.vec_to_sym(Uparams)
    #angles = MC_make_angles_from_spots_grain(R,U,agrain,pVec,wmin,wmax,tweak = False)
    return MC_Angle_Fit(R,uMat,agrain,pVec,wmin,wmax,angles, weighting)

def MC_fitR_Angle_wrap(Rparams,U,agrain,pVec,wmin,wmax,angles,weighting = True):
    R = Rotations.rotMatOfExpMap(Rparams)
    #angles = MC_make_angles_from_spots_grain(R,U,agrain,pVec,wmin,wmax,tweak = False)
    return MC_Angle_Fit(R,U,agrain,pVec,wmin,wmax,angles,weighting)
def MC_fitPrecession_Angle_wrap(pVec,R,U,agrain,wmin,wmax,angles,weighting = True):
    #angles = MC_make_angles_from_spots_grain(R,U,agrain,pVec,wmin,wmax,tweak = False)
    return MC_Angle_Fit(R,U,agrain,pVec,wmin,wmax,angles,weighting)

def MC_Angle_Fit(R,U,pVec,detectorGeom, planeData, omeMin,omeMax,angles,tweak = False,weighting = False,stdev = rad(.1)):
    angles0 = MC_make_angles_(R,U,pVec,detectorGeom,planeData,omeMin,omeMax,tweak = False)
    out = []
    key1 = angles.keys()
    key2 = angles0.keys()
    key1.sort()
    key2.sort()
    print "equalkeys?", key1==key2, len(angles.keys()), len(angles0.keys())

    for key in key1:
        if angles0.has_key(key):
            meas1,unc1 = angles[key]
            pred1,delm = angles0[key]
            diff = meas1 - pred1
            out.append(diff[0]/w1)
            out.append(diff[1]/w2)
            out.append(diff[2]/w3)
            
    for sol in angles:
        tth,eta,w = sol
        
        
    return num.array(out)


def MC_make_vectors_from_spots_grain(R,U,gwInstance,pVec,wmin,wmax, tweak = False):
    angles = MC_make_angles_from_spots_grain(R,U,gwInstance,pVec,wmin,wmax,tweak)
    wavelength = gwInstance.planeData.wavelength
    from data_class import inv_dict
    vecs = inv_dict()
    for key in angles.keys():
        if angles.triggerdict[key]==1:
            for sol in angles[key]:
                (tth,eta,w), mu_uncs = sol
                rI = makeARecipVector((tth,eta,w), wavelength)
                vecs.Add(key,[rI,mu_rI])
        else:
            (tth,eta,w), mu_uncs = sol
            rI = makeARecipVector((tth,eta,w), wavelength)
            vecs.Add(key,[rI,mu_rI])
    return vecs

def mc_replicate_jvb(R0,U0,pVec0,detectorGeom,planeData,omeMin,omeMax,stdev = rad(.1),tweak = False):
    import detector
    import spotFinder
    import grain_wrap_precession as gw
    import grain as G
    import valUnits
    
    angles = MC_make_angles_jvb(R0,U0,pVec0,detectorGeom, planeData,omeMin,omeMax,stdev = stdev, tweak = tweak)
    spotsForFit = spotFinder.Spots(planeData, angles.T, detectorGeom, (omeMin, omeMax))
    grain = G.Grain(spotsForFit,
                    rMat=R0,
                    pVec=pVec0,
                    etaTol=valUnits.valWUnit('etaTol', 'angle', 1.0, 'degrees'),
                    omeTol=valUnits.valWUnit('omeTol', 'angle', 1.0, 'degrees'))
    grain_wrap = gw.grainStrainAnalysis(grain)
    grain_wrap.assignUncertainties_angles(1e-4,1e-4,1e-2,1e-5,1e-5,1e-3)
    return grain_wrap, grain, angles
    
def mc_replicate(R0,U0,pVec0,detectorGeom,planeData,omeMin,omeMax,stdev = rad(.1),tweak = False):
    import detector
    import spotFinder
    import grain_wrap_precession as gw
    import grain as G
    import valUnits
    
    angles = MC_make_angles_(R0,U0,pVec0,detectorGeom, planeData,omeMin,omeMax,stdev = stdev, tweak = tweak)
    all_angs = []
    for key in angles.keys():
        angsol, unc = angles[key]
        all_angs.append(angsol)

    all_angs_ = num.array(all_angs)
    spotsForFit = spotFinder.Spots(planeData, all_angs_, detectorGeom, (omeMin, omeMax))
    grain = G.Grain(spotsForFit,
                    rMat=R0,
                    pVec=pVec0,
                    etaTol=valUnits.valWUnit('etaTol', 'angle', 1.0, 'degrees'),
                    omeTol=valUnits.valWUnit('omeTol', 'angle', 1.0, 'degrees'))
    grain_wrap = gw.grainStrainAnalysis(grain)
    grain_wrap.assignUncertainties_angles(1e-4,1e-4,1e-2,1e-5,1e-5,1e-3)
    return grain_wrap, grain, all_angs_
def MC_make_angles_jvb(R,U,pVec,detectorGeom,planeData,omeMin,omeMax,tweak = False, stdev = rad(0.1)):
    
    wlen    = planeData.get_wavelength()
    refCell = planeData.getLatticeOperators()['dparms']
    refBmat = planeData.getLatticeOperators()['B']
    laueGroup = 'd3d'
    tThMax = valUnits.valWUnit('tThMax', 'angle',  detectorGeom.getTThMax(), 'radians')
    
    # set exclusions and strain mag in planeData
    planeData.exclusions = num.where(planeData.getTTh() > tThMax.getVal('radians'))[0]
    planeData.strainMag  = 1e-2
    
    # full set of symmetric HKLs subject to exclusions
    fHKLs = num.hstack(planeData.getSymHKLs())
    bMat = refBmat
    qVec, qAng0, qAng1 = planeData.makeScatteringVectors(fHKLs, R, bMat, wlen)
    validAng0 = num.zeros(qVec.shape[1], dtype='bool')
    validAng1 = num.zeros(qVec.shape[1], dtype='bool')
    if hasattr(omeMin, '__len__'):
        for j in range(len(omeMin)):
            validAng0 = validAng0 | (
                ( qAng0[2, :] >= omeMin[j].getVal('radians') ) & ( qAng0[2, :] <= omeMax[j].getVal('radians') ) )
            validAng1 = validAng1 | (
                ( qAng1[2, :] >= omeMin[j].getVal('radians') ) & ( qAng1[2, :] <= omeMax[j].getVal('radians') ) )
    else:
        validAng0 = ( qAng0[2, :] >= omeMin.getVal('radians') ) & ( qAng0[2, :] <= omeMax.getVal('radians') )
        validAng1 = ( qAng1[2, :] >= omeMin.getVal('radians') ) & ( qAng1[2, :] <= omeMax.getVal('radians') )    
    tmpAngs = num.hstack( [
        qAng0[:, validAng0], qAng1[:, validAng1] ] )
    tmpDG = detector.DetectorGeomGE(detectorGeom, pVec=pVec)    
    tmpX, tmpY, tmpO = tmpDG.angToXYO(tmpAngs[0, :], tmpAngs[1, :], tmpAngs[2, :])
    tmpTTh, tmpEta, tmpOme = detectorGeom.xyoToAng(tmpX, tmpY, tmpO)
    
    tmpAngs = num.vstack([tmpTTh, tmpEta, tmpOme])
    if tweak:
        pert = numpy.random.normal(num.zeros(tmpAngs.shape), num.zeros(tmpAngs.shape)+stdev)
        tmpAngs = tmpAngs+pert
        
    spotAngs = tmpAngs
    
    return spotAngs
def mc_replicate_nrb(R0,U0,pVec0,detectorGeom,planeData,omeMin,omeMax,tweak = False, stdev = rad(.1),angle_stdev_fact = 10):
    import detector
    import spotFinder
    import grain_wrap_precession as gw
    import grain as G
    import valUnits
    print 'making angles'
    angles,ang_dev = MC_make_angles_nrb(R0,U0,pVec0,detectorGeom, planeData,omeMin,omeMax,stdev = stdev, tweak = tweak)
    spotsForFit = spotFinder.Spots(planeData, angles, detectorGeom, (omeMin, omeMax))
    grain = G.Grain(spotsForFit,
                    rMat=R0,
                    pVec=pVec0,
                    etaTol=valUnits.valWUnit('etaTol', 'angle', 1.0, 'degrees'),
                    omeTol=valUnits.valWUnit('omeTol', 'angle', 1.0, 'degrees'))
    grain_wrap = gw.grainStrainAnalysis(grain)
    #grain_wrap.assignUncertainties_angles(1e-4,1e-4,1e-2,1e-5,1e-5,1e-3)
    an1 = angle_stdev_fact*stdev
    an2 = angle_stdev_fact*stdev
    an3 = angle_stdev_fact*10*stdev
    st1 = .1*an1
    st2 = .1*an2
    st3 = .1*an3
    grain_wrap.assignUncertainties_angles(an1,an2,an3,st1,st2,st3)
    return grain_wrap, grain, angles
def MC_make_angles_nrb(R,U,pVec,detectorGeom,planeData,omeMin,omeMax,tweak = False, stdev = rad(.1)):
    angl,axxx = Rotations.angleAxisOfRotMat(R)
    Rparams0 = (angl*axxx).flatten()
    #import uncertainty_analysis
    
    tmpDG = detector.DetectorGeomGE(detectorGeom, pVec=pVec)
    #print 'recip vects...',
    recips = planeData.makeRecipVectors(R = R, U = U)
    #print 'made'
    wavelength = planeData.wavelength
    nvecs = recips.shape[1]
    from data_class import inv_dict
    #allangs = inv_dict()
    allangs = []
    #import sys
    #print "nvecs",nvecs
    for j in range(10):
        print j,
    #print "done"
    for i in range(nvecs):
        #print i
        #sys.stdout.write(str(i)+ " ")
        rI = recips[:,i]
        #print "rI",rI,tweak
        if tweak:
            Rparams_tweak = numpy.random.normal([0,0,0],[stdev,stdev,stdev])
            #Rparams = Rparams0 + Rparams_tweak
            Rtweak = Rotations.rotMatOfExpMap(Rparams_tweak)
            #print Rparams0, Rparams
            
            rI = num.dot(Rtweak,rI)
        
        angle_sols  = Get_LLNL_Angles(rI,wavelength)
        #print rI, angle_sols
        for key in angle_sols.keys():
            tth,eta,w = angle_sols[key]
            w = mapOme(w,-num.pi,num.pi)
            if eta!=num.pi and eta!=-num.pi:
                eta = mapOme(eta,-num.pi,num.pi)
            if w<omeMin.getVal('radians') or w>omeMax.getVal('radians'):
                continue
            if w is not None:
                #pred_angle_vector = num.array(detectorGeom.xyoToAng(*(px,py,w)))
                tmpX, tmpY, tmpO = tmpDG.angToXYO([tth,tth],[eta,eta],[w,w])
                tmpX = tmpX[0]
                tmpY = tmpY[0]
                tmpO = tmpO[0]
                
                tmpTTh, tmpEta, tmpOme = detectorGeom.xyoToAng(tmpX, tmpY, tmpO)
                tmpAngs = num.vstack([tmpTTh, tmpEta, tmpOme])
                #allangs.Add((i,key),(num.array((tmpTTh,tmpEta,tmpOme)),(stdev,stdev,stdev)))
                allangs.append([tmpTTh,tmpEta,tmpOme])
    return num.array(allangs)

    
def MC_make_angles_(R,U,pVec,detectorGeom,planeData,omeMin,omeMax,tweak = False, stdev = rad(.1)):
    
    tmpDG = detector.DetectorGeomGE(detectorGeom, pVec=pVec)    
    recips = planeData.makeRecipVectors(R = R, U = U)
    wavelength = planeData.wavelength
    nvecs = recips.shape[1]
    from data_class import inv_dict
    allangs = inv_dict()
    for i in range(nvecs):
        rI = recips[:,i]
        angle_sols  = Get_LLNL_Angles(rI,wavelength)
        for key in angle_sols.keys():
            tth,eta,w = angle_sols[key]
            w = mapOme(w,-num.pi,num.pi)
            if eta!=num.pi and eta!=-num.pi:
                eta = mapOme(eta,-num.pi,num.pi)
            if w<omeMin.getVal('radians') or w>omeMax.getVal('radians'):
                continue
            if w is not None:
                #pred_angle_vector = num.array(detectorGeom.xyoToAng(*(px,py,w)))
                tmpX, tmpY, tmpO = tmpDG.angToXYO([tth,tth],[eta,eta],[w,w])
                tmpX = tmpX[0]
                tmpY = tmpY[0]
                tmpO = tmpO[0]
                if tweak:
                    tmpX = tmpX + num.random.normal(0,stdev)
                    tmpY = tmpY + num.random.normal(0,stdev)
                    tmpO = tmpO + num.random.normal(0,stdev)
                tmpTTh, tmpEta, tmpOme = detectorGeom.xyoToAng(tmpX, tmpY, tmpO)
                tmpAngs = num.vstack([tmpTTh, tmpEta, tmpOme])
                allangs.Add((i,key),(num.array((tmpTTh,tmpEta,tmpOme)),(stdev,stdev,stdev)))
                #fangs.append([tmpTTh,tmpEta,tmpOme])
    return allangs

def MC_make_angles_from_spots_grain(R,U,gwInstance,pVec,wmin,wmax,tweak = False):
    pxmin,pxmax = 0,2048
    pymin,pymax = 0,2048
    tmpDG = gwInstance.detectorGeom.__class__(gwInstance.detectorGeom)
    tmpDG.pVec = pVec
    wavelength = gwInstance.planeData.wavelength
    F = num.dot(R,U)
    FinvT = inv(F).T
    masterReflInfo = gwInstance.grainSpots
    hitRelfId = num.where(masterReflInfo['iRefl'] >= 0)[0]
    measHKLs = masterReflInfo['hkl'][hitRelfId, :]
    gIs = num.dot(gwInstance.bMat, measHKLs.T)
    pred_rIs = num.dot(FinvT, gIs)
    import data_class
    angles = data_class.inv_dict()
    hits =0
    for i in range(len(hitRelfId)):
        try:
            Id = hitRelfId[i]
            spotId = masterReflInfo['iRefl'][Id]
            rI = pred_rIs[:,i]

            'measured angles and uncertainties'
            spot = gwInstance.spots._Spots__spots[spotId]
            angCOM, angCOM_unc = spot.angCOM(useFit=True, getUncertainties=True)
            xyoCOM = num.array(spot.detectorGeom.angToXYO(*angCOM)).flatten()
            new_angCOM = tmpDG.xyoToAng(*xyoCOM)
            #print 'newang old',new_angCOM,angCOM
            mus     = num.hstack(new_angCOM).flatten()
            mu_uncs = num.hstack(angCOM_unc).flatten()
            
            meas_tth, meas_eta, meas_ome = mus
            meas_angle_vector = mus
            
            'predicted angles computation'
           
            rI_xyosol = Get_LLNL_Pixel_Sols_From_Recip_Vector(rI,tmpDG,pVec,wavelength)
            # rI_xyosol = Get_LLNL_Angles(rI,wavelength)
            hits =0
            for key in rI_xyosol.keys():
                w,px,py=  rI_xyosol[key]

                if px<pxmin or px>pxmax:
                    continue
                if py<pymin or py>pymax:
                    continue
                #if eta==num.pi or w ==num.pi:
                #    continue
                'here, eta, w are defined in the range [0,2pi] constrained to be counterclockwise, to match'
                w = mapOme(w,-num.pi,num.pi)
                #eta = mapOme(eta,-num.pi,num.pi)
                if w<wmin or w>wmax:
                    continue
                #print 'meas sol', xyoCOM, (px,py,w)
                if w is not None:
                    pred_angle_vector = num.array(tmpDG.xyoToAng(*(px,py,w)))
                    #pred_angle_vector = num.array([tth,eta,w])
                   
                    
                    'found predicted solution matching measured'
                    print 'pred',pred_angle_vector, mus
                    pred_sol = pred_angle_vector
                    if tweak:
                        m1,m2,m3 = pred_angle_vector
                        u1,u2,u3 = mu_uncs
                        #m1 = m1+num.random.normal(0,abs(u1)*tweak)
                        #m2 = m2+num.random.normal(0,abs(u2)*tweak)
                        #m3 = m3+num.random.normal(0,abs(u3)*tweak)
                        m1 = m1+num.random.normal(0,tweak)
                        m2 = m2+num.random.normal(0,tweak)
                        m3 = m3+num.random.normal(0,tweak)
                        pred_angle_vector = num.array([m1,m2,m3])
                    angles.Add(Id, [pred_angle_vector,mu_uncs])
                    hits=1

                    
                    
        except RuntimeError:
#            print 'runtimeerror'
            continue
        if hits==0:
            print 'missing hit',Id,spotId
            
    return angles



    
def MC_make_angles_grain(R,U, grainInstance, sigmas,wmin,wmax):
    #assert len(gIs) == len(sigmas)
    F = dot(R,U)
    FinvT = inv(F).T
    
    masterReflInfo = grainInstance.grainSpots
    hitRelfId = num.where(masterReflInfo['iRefl'] >= 0)[0]
    measHKLs = masterReflInfo['hkl'][hitRelfId, :]
    
    pred_rIs = num.dot( FinvT, num.dot(grainInstance.bMat, measHKLs.T) )
    pred_rIs = UnMatlab(pred_rIs)
    rIs = []
    wavelength = grainInstance.planeData.wavelength
    angles = []

    for i in range(len(pred_rIs)):
        rI = pred_rIs[i]
        rI_angle_sol = Get_LLNL_Angles(rI, wavelength)
        ct = 0
        (s_tth), (s_eta),(s_w) = sigmas[i]        
        r_tth = float(numpy.random.normal(0, abs(s_tth), 1))
        r_eta = float(numpy.random.normal(0, abs(s_eta), 1))
        r_w = float(numpy.random.normal(0, abs(s_w), 1))
        
        for key in rI_angle_sol.keys():
            
            tth,eta,w = rI_angle_sol[key]
            
            w = mapOme(w,wmin,wmax)
            eta = mapOme(eta,-num.pi,num.pi)
            if eta==num.pi:
                print 'etapi'
            if w is not None and eta is not None:
                ct+=1
                
                hit = [tth+r_tth,eta+r_eta,w + r_w]
                angles.append(hit)
                
                if ct==2:
                    print 'double hit', i, rI
        if ct==0:
            print 'no solution',i,rI,rI_angle_sol
    return angles
       

def MC_make_rIs(Rparams, gIs, sigma):
    rIs = []
    R = Rotations.rotMatOfExpMap(Rparams)
    for i in range(len(gIs)):
        gI = gIs[i]
        rI_mus = num.dot(R,gI)
        errors = numpy.random.normal(0, sigma, 3)
        rI = rI_mus + errors
        rIs.append(rI)
    return rIs
def MC_make_rIs_grain(Rparams, grainInstance, scale = 1, sigma = 0):
    rIs = []
    R = Rotations.rotMatOfExpMap(Rparams)
    masterReflInfo = grainInstance.grainSpots
    hitRelfId = num.where(masterReflInfo['iRefl'] >= 0)[0]
    measHKLs = masterReflInfo['hkl'][hitRelfId, :]
    #measAngs = masterReflInfo['measAngles'][hitRelfId, :]
    #measQvec = grain.makeScatteringVectors(measAngs[:, 0], measAngs[:, 1], measAngs[:, 2])
    #predQvec = num.dot( R, mUtil.unitVector( num.dot(grainInstance.bMat, measHKLs.T) ) )
    pred_rIs = num.dot( R, num.dot(grainInstance.bMat, measHKLs.T) )
    pred_rIs = UnMatlab(pred_rIs)
    err_rIs = []
    for i in range(len(pred_rIs)):
        rI = pred_rIs[i]
        if sigma == 0:
            errors = num.zeros(3)
        else:
            errors = numpy.random.normal(0, sigma, 3)
        rI = rI + errors
        err_rIs.append(rI/scale)
    #rIS_2 = grainInstance.planeData.makeRecipVectors(R = R)
    return err_rIs


def fitOrientation_func(Rparams, gIs, rIs):
    R = Rotations.rotMatOfExpMap(Rparams)
    residual = []
    for i in range(len(gIs)):
        gI = gIs[i]
        rI = rIs[i]
        diff = num.dot(R,gI) - rI
        residual.append(diff[0])
        residual.append(diff[1])
        residual.append(diff[2])
    residual = num.array(residual)
    return residual
def UnMatlab(col_vecs):
    out = []
    for i in range(col_vecs.shape[1]):
        out.append(col_vecs[:,i])
    return out
"""
from Vector_funcs import *
angstrom = 1e-10
lam = planeData.wavelength
gIs = (1./lam)*planeData.makeRecipVectors(R = num.eye(3))
gIs = UnMatlab(gIs)
axis = erhofunc(.4,.5)
angle = .3
R = Rotations.rotMatOfExpMap(angle*axis)
sigma = .1 * lam

rIs_err = MC_make_rIs(angle*axis, gIs, sigma)
rIs_0= MC_make_rIs_grain(angle*axis, rubyGrain, sigma)
from scipy.optimize import leastsq
#out = leastsq(fitOrientation_func, x0 = [.1,.1,.1], args = (gIs, rIs_err), full_output = 1) 



masterReflInfo = rubyGrain.grainSpots
hitRelfId = num.where(masterReflInfo['iRefl'] >= 0)[0]
measHKLs = masterReflInfo['hkl'][hitRelfId, :]
measAngs = masterReflInfo['measAngles'][hitRelfId, :]
measQvec = grain.makeScatteringVectors(measAngs[:, 0], measAngs[:, 1], measAngs[:, 2])

zeroangle = 0
gIs= MC_make_rIs_grain(zeroangle*axis, rubyGrain, scale = 1e-10, sigma = 0)
Rmeasurements = []
n_replications = 100
sigma = .1 * lam
cov_xs = []
for i in range(n_replications):
    rIs_err = MC_make_rIs_grain(angle*axis, rubyGrain, scale = 1e-10, sigma = sigma)
    out = leastsq(fitOrientation_func, x0 = [.1,.1,.1], args = (gIs, rIs_err), full_output = 1)
    Rmeasurements.append(out[0])
    cov_xs.append(out[1])

    
Rmeasurements = num.array(Rmeasurements)
r1 = Rmeasurements[:,0]
r2 = Rmeasurements[:,1]
r3 = Rmeasurements[:,2]
N = len(r1)
cov = num.zeros((len(out[0]),len(out[0])))
rMean = num.mean(Rmeasurements,0)

for i in range(len(r1)):
    Ri = Rmeasurements[i]
    cov += num.outer((Ri - rMean),(Ri - rMean))
    
    

cov_matrix = (1./(N-1)) * cov
'or'
cov_ = numpy.cov(Rmeasurements,rowvar = 0)

z = num.zeros(cov_.shape)
for i in range(len(cov_xs)):
    z+=cov_xs[i]
z/=len(cov_xs)
    
bias = rMean - angle*axis
"""
"""
#getting uncertainties
sigmas = []
rubyGrain.fitGrainSpotUncertainties()
spots_ = rubyGrain.spots._Spots__spots
for spot in spots_:
    try:
        
        
        angCOM, angCOM_unc = spot.angCOM(useFit=True, getUncertainties=True)
        sigmas.append(angCOM_unc)
    except RuntimeError:
        sigmas.append([1e-3,1e-3,1e-3])
        

angs = MC.MC_make_angles_grain(rubyGrain.rMat,rubyGrain.uMat,rubyGrain,sigmas, wmin,wmax)

spotsForFit2 = spotFinder.Spots(rubyGrain.planeData, num.array(angs), detectorGeom, [wmin, wmax])


rubyGrain_sim = G.Grain(spotsForFit2,
                        rMat=rMat,
                        etaTol=valUnits.valWUnit("etaTol","ANGLE",0.25,"degrees"),
                        omeTol=valUnits.valWUnit("omeTol","ANGLE",0.50,"degrees") )
rubyGrain_sim.refineGrain()

"""

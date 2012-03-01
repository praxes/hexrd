# ============================================================
# Copyright (c) 2012, Lawrence Livermore National Security, LLC. 
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
"""friend functions for grain.py's Grain class"""

import numpy as num
from numpy import dot,sqrt
from scipy import optimize
from scipy.linalg import inv
from scipy.linalg import eig

from hexrd.XRD import Rotations 
from hexrd.Vector_Data_Structures_LLNL import Get_Pixel_Polar_Coords,\
    alphabeta_polar_to_etatheta, Get_LLNL_Angles, mapOme
from hexrd.Vector_funcs import Mag
from hexrd.Vector_funcs import Unit,Star,Inner_Prod
from hexrd.Vector_funcs import polarDecomposition
from hexrd.XRD import uncertainty_analysis

def vec_to_sym(UV_vec):
    UV_vec = UV_vec.flatten()
    UV_11,UV_22,UV_33,UV_12,UV_23,UV_13 = UV_vec
    UV = num.array([[UV_11,UV_12,UV_13],
                    [UV_12,UV_22,UV_23],
                    [UV_13,UV_23,UV_33]])
    return UV
def sym_to_vec(A):
    A11,A22,A33 = A[0,0],A[1,1],A[2,2]
    A12,A23,A13 = A[0,1],A[1,2],A[0,2]
    return num.array([A11,A22,A33,A12,A23,A13])

def getUfromVR(V,R):
    U = num.dot(num.dot(R.T, V),R)
    return U
def getVfromUR(U,R):
    V = num.dot(num.dot(R, U),R.T)
    return V

def makeARecipVector((tTh,eta,ome),wavelength):
    return makeRecipVectors([tTh],[eta],[ome],wavelength)[0]

def makeRecipVectors(tTh, eta, ome, wavelength):
    """
    """
    rIs = []
    for i in range(len(tTh)):
        t,e,w = tTh[i]/2., eta[i], ome[i]
        rI_hat = num.array([num.cos(w)*num.cos(e)*num.cos(t) - num.sin(w)*num.sin(t),
                            num.sin(e)*num.cos(t),
                            num.sin(w)*num.cos(e)*num.cos(t) + num.cos(w)*num.sin(t) ])
        d = wavelength/(2.*num.sin(t))
        rI = (1./d)*rI_hat
        rIs.append(rI)
    return num.array(rIs)

class grainStrainAnalysis:
    """
    essentially a friend class
    """
    def __init__(self, agrain):
        self.grain = agrain
        self.haveGrainSpotUncertainties = False
        self.grainSpots = agrain.grainSpots
        self.spots = agrain.spots
        self.planeData = agrain.planeData
        self.vecs_spots_associated = False
        self.vecs_associated = False
        self.gI_spots = None
        self.gI_rIs_u_rIs = None
        self.gI_rIs = None
        self.bMat = agrain.bMat
        self.rMat = agrain.rMat
        self.uMat = agrain.uMat
        self.confidence_level = agrain.confidence_level
    def reset(self):
        self.rMat = self.grain.rMat
        self.uMat = self.grain.uMat
        self.vecs_spots_associated = False
        self.vecs_associated = False
        self.gI_spots = None
        self.gI_rIs_u_rIs = None
        self.gI_rIs = None
       
        
    def fitGrainSpotUncertainties(self):
        if self.haveGrainSpotUncertainties==True:
            return
        masterReflInfo = self.grain.grainSpots
        hitRelfId = num.where(masterReflInfo['iRefl'] >= 0)[0]
        for i in range(len(hitRelfId)):
            Id = hitRelfId[i]
            spotId = masterReflInfo['iRefl'][Id]
            spot = self.spots._Spots__spots[spotId]
            spot.cleanFit()
            spot.fitWrap(uncertainties = True, confidence_level = self.confidence_level);
            try:
                angCOM, angCOM_unc = spot.angCOM(useFit=True, getUncertainties=True)
            except RuntimeError:
                continue
            #print angCOM,angCOM_unc
        self.haveGrainSpotUncertainties = True
        
    def associateRecipVectors(self, forceReRead = False):
        """
        pairs up gI from B.[h,k,l] to the measured rI for that [h,k,l].
        For use with rotation/strain refinement functions.
        stores results in self.gI_rIs
        """
        if self.vecs_associated == True and forceReRead == False:
            return self.gI_rIs
        'else do the association'
        
        masterReflInfo = self.grainSpots
        hitRelfId = num.where(masterReflInfo['iRefl'] >= 0)[0]
        measHKLs = masterReflInfo['hkl'][hitRelfId, :]
        measAngs = masterReflInfo['measAngles'][hitRelfId, :]
        'measAngs = [tTh, eta, ome]'
        meas_rIs = makeRecipVectors(measAngs[:, 0], measAngs[:, 1], measAngs[:, 2],self.planeData.wavelength)
        'make predicted gI vectors'
        'these are organized by columns'
        gIs = num.dot(self.bMat, measHKLs.T)
        'make associations'
        gI_rIs=[]
        for i in range(len(meas_rIs)):
            rI = meas_rIs[i]
            gI = gIs[:,i]
            gI_rIs.append([gI,rI])

        self.gI_rIs = gI_rIs
        self.vecs_associated = True
        return gI_rIs
    
    def associateRecipVector_Spots(self,forceReRead = False):
        """
        pairs up gI from B.[h,k,l] to the measured rI for that [h,k,l].
        For use with rotation/strain refinement functions.
        stores results in self.gI_rIs
        """
        if self.vecs_spots_associated == True and forceReRead == False:
            return self.gI_spots
        'else do the association'
        
        masterReflInfo = self.grainSpots
        hitRelfId = num.where(masterReflInfo['iRefl'] >= 0)[0]
        measHKLs = masterReflInfo['hkl'][hitRelfId, :]
        measAngs = masterReflInfo['measAngles'][hitRelfId, :]
        'make predicted gI vectors'
        'these are organized by columns'
        gIs = num.dot(self.bMat, measHKLs.T)
        'make associations'
        gI_spots=[]
        for i in range(len(hitRelfId)):
            Id = hitRelfId[i]
            spotId = masterReflInfo['iRefl'][Id]
            spot = self.spots._Spots__spots[spotId]
            gI = gIs[:,i]
            gI_spots.append([gI,spot])

        self.gI_spots = gI_spots
        self.vecs_spots_associated = True
        return gI_spots
    def associateRecipVectorAngles(self, R, U):
        """
        """
        self.fitGrainSpotUncertainties()

        wavelength = self.planeData.wavelength

        F = num.dot(R,U)
        FinvT = inv(F).T
        
        masterReflInfo = self.grainSpots
        hitRelfId = num.where(masterReflInfo['iRefl'] >= 0)[0]
        measHKLs = masterReflInfo['hkl'][hitRelfId, :]

        gIs = num.dot(self.bMat, measHKLs.T)
        pred_rIs = num.dot(FinvT, gIs)
        
        angles = []
        for i in range(len(hitRelfId)):
            try:
                Id = hitRelfId[i]
                spotId = masterReflInfo['iRefl'][Id]
                rI = pred_rIs[:,i]

                'measured angles and uncertainties'
                spot = self.spots._Spots__spots[spotId]
                angCOM, angCOM_unc = spot.angCOM(useFit=True, getUncertainties=True)
                mus     = num.hstack(angCOM).flatten()
                mu_uncs = num.hstack(angCOM_unc).flatten()

                meas_tth, meas_eta, meas_ome = mus
                meas_angle_vector = mus

                'predicted angles computation'
                rI_angle_sol = Get_LLNL_Angles(rI, wavelength)
                for key in rI_angle_sol.keys():
                    tth,eta,w = rI_angle_sol[key]
                    if eta==num.pi or w ==num.pi:
                        continue
                    'here, eta, w are defined in the range [0,2pi] constrained to be counterclockwise, to match'
                    w = mapOme(w,-num.pi,num.pi)
                    eta = mapOme(eta,-num.pi,num.pi)
                    
                    if w is not None:
                        pred_angle_vector = num.array([tth,eta,w])
                        #print pred_angle_vector, meas_angle_vector
                        diff = Mag(meas_angle_vector - pred_angle_vector)
                        if diff<.01:
                            #print 'diff angle 0 ', Mag(meas_angle_vector - pred_angle_vector)
                            'found predicted solution matching measured'
                            pred_sol = pred_angle_vector
                            angles.append([meas_angle_vector,pred_angle_vector,mu_uncs])
                            break

            except RuntimeError:
                #print 'runtimeerror'
                continue
            

        return angles
    def associateRecipVectors_w_Uncertainty(self, forceReRead = False):
        """
        pairs up gI from B.[h,k,l] to the measured rI for that [h,k,l].
        For use with rotation/strain refinement functions.
        stores results in self.gI_rIs
        """
        if self.vecs_associated_w_uncertainty == True and forceReRead == False:
            return self.gI_rIs_u_rIs
        'else do the association'
        self.fitGrainSpotUncertainties()
        wavelength = self.planeData.wavelength
        
        masterReflInfo = self.grainSpots
        hitRelfId = num.where(masterReflInfo['iRefl'] >= 0)[0]
        measHKLs = masterReflInfo['hkl'][hitRelfId, :]
        gIs = num.dot(self.bMat, measHKLs.T)
        
        gI_rIs_u_rIs=[]
        for i in range(len(hitRelfId)):
            Id = hitRelfId[i]
            spotId = masterReflInfo['iRefl'][Id]
            spot = self.spots._Spots__spots[spotId]
            try:
                angCOM, angCOM_unc = spot.angCOM(useFit=True, getUncertainties=True)
            except RuntimeError:
                continue
            mus     = num.hstack(angCOM).flatten()
            mu_uncs = num.hstack(angCOM_unc).flatten()
            
            meas_tth, meas_eta, meas_ome = mus
            rI = makeARecipVector([meas_tth,meas_eta,meas_ome], wavelength)
            u_rI = uncertainty_analysis.propagateUncertainty(makeARecipVector, mu_uncs, 1e-8, mus,wavelength )
            gI = gIs[:,i]
            gI_rIs_u_rIs.append([gI, rI, u_rI])
        
        self.vecs_associated_w_uncertainty = True
        self.gI_rIs_u_rIs = gI_rIs_u_rIs
        return gI_rIs_u_rIs
    def computeF(self):
        fMat = num.dot(self.rMat,self.uMat)
        self.fMat = fMat
        return fMat

    def fitAngles_lsq(self,Rparams,Umat, weighting = True):
        R_ = Rotations.rotMatOfExpMap(Rparams)
        angles = self.associateRecipVectorAngles(R_,Umat)
        out = []
        for angleset in angles:
            measangles, predangles, unc = angleset
            diff = measangles - predangles
            #print measangles,predangles
            if weighting:
                w1,w2,w3 = unc
            else:
                w1=w2=w3=1.
            out.append(diff[0]/w1)
            out.append(diff[1]/w2)
            out.append(diff[2]/w3)
        return num.array(out)
    def fitU_Angles_(self, Rparams = 'default', Uparams = 'default', weighting = True, report = True,evaluate = False):
        if Uparams=='default':
            Uparams = sym_to_vec(self.uMat)
        else:
            Uparams = Uparams
        if Rparams == 'default':
            angl, axxx = Rotations.angleAxisOfRotMat(self.rMat)
            Rparams = (angl*axxx).flatten()
        else:
            Rparams = Rparams

        def _lsq_U(Uparams, Rparams):
            U = vec_to_sym(Uparams)
            return self.fitAngles_lsq(Rparams, U, weighting)
        if evaluate:
            return _lsq_U(Uparams, Rparams, weighting)
        optResults = optimize.leastsq(_lsq_U, x0 = Uparams, args = Rparams, full_output = 1)
        r1 = optResults[0]
        U1 = vec_to_sym(r1)
        if report:
            print "refined U matrix: \n%g, %g, %g\n%g, %g, %g\n%g, %g, %g\n" \
                  % (U1[0, 0], U1[0, 1], U1[0, 2], \
                     U1[1, 0], U1[1, 1], U1[1, 2], \
                     U1[2, 0], U1[2, 1], U1[2, 2], )
        self.uMat = U1
        return optResults

    def fitRotation_Angles_(self, Rparams = 'default', Uparams = 'default', report = True, weighting = True, evaluate = False):
        """
        """
        if Rparams == 'default':
            angl,axxx = Rotations.angleAxisOfRotMat(self.rMat)
            Rparams = (angl*axxx).flatten()
        else:
            Rparams = Rparams
        if Uparams=='default':
            Uparams = sym_to_vec(self.uMat)
        else:
            Uparams = Uparams
            
        U = vec_to_sym(Uparams)

        def _lsq_Rotation(Rparams, U, weighting = True):
            return self.fitAngles_lsq(Rparams, U, weighting)
        if evaluate:
            return _lsq_Rotation(Rparams,U,weighting)
        optResults = optimize.leastsq(_lsq_Rotation, x0 = Rparams, args = (U,weighting), full_output = 1)
        r1 = optResults[0]
        U1 = Rotations.rotMatOfExpMap(r1)
        if report:
            print "refined orientation matrix: \n%g, %g, %g\n%g, %g, %g\n%g, %g, %g\n" \
                  % (U1[0, 0], U1[0, 1], U1[0, 2], \
                     U1[1, 0], U1[1, 1], U1[1, 2], \
                     U1[2, 0], U1[2, 1], U1[2, 2], )
        
        # reset pVec in self
        self.rMat = U1
        return optResults
    def fitRotation(self, Rparams = 'default', Uparams = 'default', report = True,evaluate = False):
        """
        uses F*N = \alpha n to refine rotation
        """
        if Rparams == 'default':
            angl,axxx = Rotations.angleAxisOfRotMat(self.rMat)
            Rparams = (angl*axxx).flatten()
        else:
            Rparams = Rparams
        if Uparams=='default':
            Uparams = sym_to_vec(self.uMat)
        else:
            Uparams = Uparams

        def _fitRotation_lsq(Rparams, vecs, U):
            out = []
            Rmat = Rotations.rotMatOfExpMap(Rparams)
            Fmat = dot(Rmat,U)
            Cmat = dot(Fmat.T,Fmat)
            for r0,rd in vecs:
                N = Unit(r0)
                n = Unit(rd)
                alpha_ = sqrt(Inner_Prod(Star(Cmat),N,N))
                r_i = Inner_Prod(Star(Fmat),n,N) - alpha_
                out.append(r_i)
            return num.array(out)

        gI_rIs = self.associateRecipVectors()
        U = vec_to_sym(Uparams)
        if evaluate:
            return _fitRotation_lsq(Rparams,gI_rIs,U)
        
        optResults = optimize.leastsq(_fitRotation_lsq, x0 = Rparams, args = (gI_rIs, U), full_output = 1)
        
        r1 = optResults[0]
        U1 = Rotations.rotMatOfExpMap(r1)
        if report:
            print "refined orientation matrix: \n%g, %g, %g\n%g, %g, %g\n%g, %g, %g\n" \
                  % (U1[0, 0], U1[0, 1], U1[0, 2], \
                     U1[1, 0], U1[1, 1], U1[1, 2], \
                     U1[2, 0], U1[2, 1], U1[2, 2], )
        
        # reset pVec in self
        self.rMat = U1
        return optResults
    def fitRotationStar_angle_weights(self, Rparams = 'default', Uparams = 'default', report = True, weighting = True, evaluate = False):
        """
        uses F*N = \alpha n to refine rotation
        """
        def _fitR_residual(angles, r0, Fmat, wavelength):
            Cmat = dot(Fmat.T,Fmat)
            rI = makeARecipVector(angles,wavelength)
            n = Unit(rI)
            N = Unit(r0)
            alpha_ = sqrt(Inner_Prod(Star(Cmat),N,N))
            r_i = Inner_Prod(Star(Fmat),n,N) - alpha_
            return r_i
        
        def _fitRotation_lsq(Rparams, gI_spots, U,wavelength):
            
            out = []
            Rmat = Rotations.rotMatOfExpMap(Rparams)
            Fmat = dot(Rmat,U)
            for gI,spot in gI_spots:
                try:
                    angCOM, angCOM_unc = spot.angCOM(useFit=True, getUncertainties=True)
                except RuntimeError:
                    continue
                r_i = _fitR_residual(angCOM.flatten(),gI, Fmat, wavelength)
                if weighting:
                    weight = uncertainty_analysis.propagateUncertainty(_fitR_residual, angCOM_unc, 1e-8,angCOM.flatten(),gI,Fmat,wavelength)
                    #print 'Rotation using weight', weight, r_i, r_i/weight, angCOM_unc
                else:
                    weight = 1.
                out.append(r_i/weight)
            return num.array(out)
        
        if Rparams == 'default':
            angl,axxx = Rotations.angleAxisOfRotMat(self.rMat)
            Rparams = (angl*axxx).flatten()
        else:
            Rparams = Rparams
        if Uparams=='default':
            Uparams = sym_to_vec(self.uMat)
        else:
            Uparams = Uparams

        self.fitGrainSpotUncertainties()
        gI_spots = self.associateRecipVector_Spots()
        U = vec_to_sym(Uparams)
        wavelength = self.planeData.wavelength
        if evaluate:
            return _fitRotation_lsq(Rparams,gI_spots,U,wavelength)
        
        optResults = optimize.leastsq(_fitRotation_lsq, x0 = Rparams, args = (gI_spots, U, wavelength), full_output = 1)
        r1 = optResults[0]
        U1 = Rotations.rotMatOfExpMap(r1)
        if report:
            print "refined orientation matrix: \n%g, %g, %g\n%g, %g, %g\n%g, %g, %g\n" \
                  % (U1[0, 0], U1[0, 1], U1[0, 2], \
                     U1[1, 0], U1[1, 1], U1[1, 2], \
                     U1[2, 0], U1[2, 1], U1[2, 2], )
        
        self.rMat = U1
        return optResults
    def refineF(self, tol = 1e-9, maxiter = 100, report = True):
        """
        iteratively calls fitRotation(), fitC(), computes the difference in the results F in a two norm sense, and stops when this difference goes below the specified tolerance, or when max iterations is reached.
        """
        gI_rIs = self.associateRecipVectors()
        fMat0 = self.computeF()        
        fMat0 = fMat0.flatten()
        ct = 0
        while 1:
            self.fitRotation(report = False)
            rMat = self.rMat
            self.fitC(report = False)
            uMat = self.uMat
            fMat = num.dot(rMat,uMat)
            fMat = fMat.flatten()
            fDiff = Mag(fMat - fMat0)
            print 'difference from previous iteration', fDiff
            if fDiff < tol:
                break
            if ct>maxiter:
                break
            fMat0 = fMat
            ct+=1
        F = self.computeF()
        if report:
            print 'refined deformation gradient','\n',F

        
    def refineGrain(self, tol = 1e-9, maxiter = 100, report = True):
        """
        iteratively calls fitRotation(), fitC(), computes the difference in the results F in a two norm sense, and stops when this difference goes below the specified tolerance, or when max iterations is reached.
        """
        gI_rIs = self.associateRecipVectors()
        fMat0 = self.computeF()        
        fMat0 = fMat0.flatten()
        ct = 0
        while 1:
            self.fitRotation(report = False)
            rMat = self.rMat
            self.fitC(report = False)
            uMat = self.uMat
            fMat = num.dot(rMat,uMat)
            fMat = fMat.flatten()
            fDiff = Mag(fMat - fMat0)
            print 'difference from previous iteration', fDiff
            if fDiff < tol:
                break
            if ct>maxiter:
                break
            fMat0 = fMat
            ct+=1
        F = self.computeF()
        if report:
            print 'refined deformation gradient','\n',F
            
    def refineGrain_Angles(self, tol = 1e-9, maxiter = 100, report = True):
        """
        iteratively calls fitRotation(), fitC(), computes the difference in the results F in a two norm sense, and stops when this difference goes below the specified tolerance, or when max iterations is reached.
        """
        fMat0 = self.computeF()        
        fMat0 = fMat0.flatten()
        ct = 0
        while 1:
            self.fitRotation_Angles_(report = False)
            rMat = self.rMat
            self.fitU_Angles_(report = False)
            uMat = self.uMat
            fMat = num.dot(rMat,uMat)
            fMat = fMat.flatten()
            fDiff = Mag(fMat - fMat0)
            print 'difference from previous iteration', fDiff
            if fDiff < tol:
                break
            if ct>maxiter:
                break
            fMat0 = fMat
            ct+=1
        F = self.computeF()
        if report:
            print 'refined deformation gradient','\n',F
            
    def refineGrain_Angles_2(self, tol = 1e-9, maxiter = 100, report = True):
        """
        iteratively calls fitRotation(), fitC(), computes the difference in the results F in a two norm sense, and stops when this difference goes below the specified tolerance, or when max iterations is reached.
        """
        fMat0 = self.computeF()        
        fMat0 = fMat0.flatten()
        ct = 0
        while 1:
            self.fitRotationStar_angle_weights(weighting = True, report = False)
            rMat = self.rMat
            self.fitC_angle_weights(weighting = True, report = False)
            uMat = self.uMat
            fMat = num.dot(rMat,uMat)
            fMat = fMat.flatten()
            fDiff = Mag(fMat - fMat0)
            print 'difference from previous iteration', fDiff
            if fDiff < tol:
                break
            if ct>maxiter:
                break
            fMat0 = fMat
            ct+=1
        F = self.computeF()
        if report:
            print 'refined deformation gradient','\n',F
            

    def fitF(self):
        def _fitF(FParams, vecs):
            Rparams = FParams[0:3]
            Uparams = FParams[3:9]
            
            rMat = Rotations.rotMatOfExpMap(Rparams)
            uMat = vec_to_sym(Uparams)
            F = dot(rMat,uMat)
            FinvT = inv(F).T 

            out = []
            for r0,rd in vecs:
                diff = dot(FinvT,r0) - rd
                out.append(diff[0])
                out.append(diff[1])
                out.append(diff[2])
            return num.array(out)
        
        gI_rIs = self.associateRecipVectors()
        angle,axxx = Rotations.angleAxisOfRotMat(self.rMat)
        Rparams = (angle*axxx).flatten()
        Uparams = sym_to_vec(self.uMat)
        Fparams = []
        Fparams.extend(list(Rparams))
        Fparams.extend(list(Uparams))

        optResults = optimize.leastsq(_fitF, x0 = Fparams, args = gI_rIs, full_output = 1)
        Rparams = optResults[0][0:3]
        Uparams = optResults[0][3:9]
        rMat = Rotations.rotMatOfExpMap(Rparams)
        uMat = vec_to_sym(Uparams)
        F = dot(rMat,uMat)
                
        print 'finalF'
        print F
        
        R_,U = polarDecomposition(F)
        print R_,'\n',U
        return optResults
                
    def fitC(self, Cparams = 'default', report = True, evaluate = False):
        """
        alternate fitting approach for U, uses no Rotation information.
        Cparams: [C11,C22,C33,C12,C23,C13]
        """
        'internal objective function'
        def _fitC(Cparams,vecs):
            Cmat = vec_to_sym(Cparams)
            invCmat = inv(Cmat)
            out = []
            for r0,rd in vecs:
                df = 1./Mag(rd)
                d0 = 1./Mag(r0)
                gii_0 = dot(r0,r0)
                gii = dot(dot(invCmat,r0),r0)
                r_i = df/d0 - sqrt(gii_0)/(sqrt(gii))
                out.append(r_i)
            return num.array(out)

        if Cparams == 'default':
            U = self.uMat
            Cparams = sym_to_vec(dot(U,U))
        else:
            Cparams = Cparams

        gI_rIs = self.associateRecipVectors()
        if evaluate:
            return _fitC(Cparams, gI_rIs)
        
        optResults = optimize.leastsq(_fitC, x0 = Cparams, args = gI_rIs, full_output = 1)
        Cparams = optResults[0]

        C = vec_to_sym(Cparams)
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

        self.uMat = U
    def fitC_angle_weights(self, weighting = True,Cparams = 'default', report = True,evaluate = False):
        """
        alternate fitting approach for U, uses no Rotation information.
        Cparams: [C11,C22,C33,C12,C23,C13]
        """
        'internal objective function'
        def _fitC_residual(angles, r0, invCmat, wavelength):
            rI = makeARecipVector(angles, wavelength)
            #print rI
            df = 1./Mag(rI)
            d0 = 1./Mag(r0)
            gii_0 = dot(r0,r0)
            gii = dot(dot(invCmat,r0),r0)
            r_i = df/d0 - sqrt(gii_0)/(sqrt(gii))
            return r_i
        def _fitC(Cparams,gI_spots,wavelength):
            Cmat = vec_to_sym(Cparams)
            invCmat = inv(Cmat)
            out = []
            for gI, spot in gI_spots:
                try:
                    angCOM, angCOM_unc = spot.angCOM(useFit=True, getUncertainties=True)
                except RuntimeError:
                    continue
                r_i = _fitC_residual(angCOM.flatten(),gI, invCmat, wavelength)
                if weighting:

                    weight = uncertainty_analysis.propagateUncertainty(_fitC_residual, angCOM_unc, 1e-8,angCOM.flatten(),gI,invCmat,wavelength)
                    #print 'using weight', weight, r_i, r_i/weight, angCOM_unc
                else:
                    weight = 1
                out.append(r_i/weight)
            return num.array(out)

        if Cparams == 'default':
            U = self.uMat
            Cparams = sym_to_vec(dot(U,U))
        else:
            Cparams = Cparams


        self.fitGrainSpotUncertainties()
        gI_spots = self.associateRecipVector_Spots()
        wavelength = self.planeData.wavelength
        if evaluate:
            return _fitC(Cparams,gI_spots,wavelength)
        
        optResults = optimize.leastsq(_fitC, x0 = Cparams, args = (gI_spots,wavelength), full_output = 1)
        Cparams = optResults[0]

        C = vec_to_sym(Cparams)
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

        self.uMat = U
        return optResults

    def fitU(self, Rparams = 'default', Uparams = 'default'):
        if Uparams=='default':
            Uparams = sym_to_vec(self.uMat)
        else:
            Uparams = Uparams
        if Rparams == 'default':
            angl, axxx = Rotations.angleAxisOfRotMat(self.rMat)
            Rparams = (angl*axxx).flatten()
        else:
            Rparams = Rparams
            
        
        optResults = optimize.leastsq(self._fitU_objFunc_lsq, x0 = Uparams, args = (Rparams, self.planeData.wavelength), xtol=1e-5, ftol=1e-5, full_output = 1)
        
        Uparams = optResults[0]
        U = vec_to_sym(Uparams)
        
                
        print "refined stretch matrix: \n%g, %g, %g\n%g, %g, %g\n%g, %g, %g\n" \
              % (U[0, 0], U[0, 1], U[0, 2], \
                 U[1, 0], U[1, 1], U[1, 2], \
                 U[2, 0], U[2, 1], U[2, 2], )



        self.uMat = U

    

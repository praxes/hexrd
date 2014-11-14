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

import math

from scipy import \
     array, asarray, ndarray, \
     append, compress, \
     sqrt, sum, dot, sparse, \
     concatenate, c_, r_, ix_, hstack, vstack, \
     squeeze, tile, reshape, \
     ones, zeros, int_, \
     setxor1d, sort, union1d, unique, \
     finfo
from scipy.linalg import svd
import numpy as num

# module variables
sqr6i  = 1./sqrt(6.)
sqr3i  = 1./sqrt(3.)
sqr2i  = 1./sqrt(2.)
sqr2   = sqrt(2.)
sqr2b3 = sqrt(2./3.)

fpTol = finfo(float).eps                # ~2.2e-16
vTol  = 1.0e-14


def columnNorm(a):
    """
    normalize array of column vectors (hstacked, axis = 0)
    """
    if len(a.shape) > 2:
        raise RuntimeError, "incorrect shape: arg must be 1-d or 2-d, yours is %d" %(len(a.shape))
    
    cnrma = sqrt(sum(asarray(a)**2, 0))
    
    return cnrma

def rowNorm(a):
    """
    normalize array of row vectors (vstacked, axis = 1)
    """
    if len(a.shape) > 2:
        raise RuntimeError, "incorrect shape: arg must be 1-d or 2-d, yours is %d" %(len(a.shape))
    
    cnrma = sqrt(sum(asarray(a)**2, 1))
    
    return cnrma

def unitVector(a):
    """
    normalize array of column vectors (hstacked, axis = 0)
    """
    assert a.ndim in [1, 2], "incorrect arg shape; must be 1-d or 2-d, yours is %d-d" % (a.ndim)
    
    ztol = 1.0e-14
    
    m = a.shape[0]; n = 1
    
    nrm = tile(sqrt(sum(asarray(a)**2, 0)), (m, n))
    
    # prevent divide by zero
    zchk = nrm <= ztol
    nrm[zchk] = 1.0
    
    nrma = a/nrm
    
    return nrma

def nullSpace(A, tol=vTol):
    """
    computes the null space of the real matrix A
    """
    assert A.ndim == 2, 'input must be 2-d; yours is %d-d' % (A.ndim)
    
    n, m = A.shape

    if n > m :
        return nullSpace(A.T, tol).T

    U, S, V = svd(A)
    
    S = hstack([S, zeros(m-n)])

    null_mask  = (S <= tol)
    null_space = V[null_mask, :]

    return null_space

def blockSparseOfMatArray(matArray):
    """
    blockSparseOfMatArray

    Constructs a block diagonal sparse matrix (csc format) from a
    (p, m, n) ndarray of p (m, n) arrays
    
    ...maybe optional args to pick format type?
    """

    # if isinstance(args[0], str):
    #    a = args[0]
    # if a == 'csc': ...
        
    if len(matArray.shape) != 3:
        raise RuntimeError, "input array is not the correct shape!"
    
    l = matArray.shape[0]
    m = matArray.shape[1]
    n = matArray.shape[2]
    
    mn   = m*n;
    jmax = l*n;
    imax = l*m;
    ntot = l*m*n;
    
    rl    = asarray(range(l), 'int') 
    rm    = asarray(range(m), 'int')
    rjmax = asarray(range(jmax), 'int')
    
    sij = matArray.transpose(0, 2, 1).reshape(1, ntot).squeeze()
    j   = reshape(tile(rjmax, (m, 1)).T, (1, ntot))
    i   = reshape(tile(rm, (1, jmax)), (1, ntot)) + reshape(tile(m*rl, (mn, 1)).T, (1, ntot))
    
    ij  = concatenate((i, j), 0)
    
    smat = sparse.csc_matrix((sij, ij), shape=(imax, jmax)) # syntax as of scipy-0.7.0
    
    return smat

def symmToVecMV(A):
    """
    convert from symmetric matrix to Mandel-Voigt vector
    representation (JVB)
    """ 
    mvvec = zeros(6, dtype='float64')
    mvvec[0] = A[0,0]
    mvvec[1] = A[1,1]
    mvvec[2] = A[2,2]
    mvvec[3] = sqr2 * A[1,2]
    mvvec[4] = sqr2 * A[0,2]
    mvvec[5] = sqr2 * A[0,1] 
    return mvvec

def vecMVToSymm(A):
    """
    convert from Mandel-Voigt vector to symmetric matrix 
    representation (JVB) 
    """ 
    symm = zeros((3, 3), dtype='float64')
    symm[0, 0] = A[0]
    symm[1, 1] = A[1]
    symm[2, 2] = A[2]
    symm[1, 2] = A[3] / sqr2
    symm[0, 2] = A[4] / sqr2
    symm[0, 1] = A[5] / sqr2 
    symm[2, 1] = A[3] / sqr2
    symm[2, 0] = A[4] / sqr2
    symm[1, 0] = A[5] / sqr2 
    return symm

def vecMVCOBMatrix(R):
    """
    GenerateS array of 6 x 6 basis transformation matrices for the
    Mandel-Voigt tensor representation in 3-D given by: 
    
    [A] = [[A_11, A_12, A_13],
           [A_12, A_22, A_23],
           [A_13, A_23, A_33]]
    
    {A} = [A_11, A_22, A_33, sqrt(2)*A_23, sqrt(2)*A_13, sqrt(2)*A_12]
              
    where the operation :math:`R*A*R.T` (in tensor notation) is obtained by
    the matrix-vector product [T]*{A}.
    
    USAGE
    
        T = vecMVCOBMatrix(R)
    
    INPUTS
    
        1) R is (3, 3) an ndarray representing a change of basis matrix
    
    OUTPUTS
    
        1) T is (6, 6), an ndarray of transformation matrices as
           described above
    
    NOTES
    
        1) Compoments of symmetric 4th-rank tensors transform in a
           manner analogous to symmetric 2nd-rank tensors in full
           matrix notation. 

    SEE ALSO

    symmToVecMV, vecMVToSymm, quatToMat
    """
    rdim = len(R.shape)
    if rdim == 2:
        nrot = 1
        R = tile(R, (1, 1, 1))
    elif rdim == 3:
        nrot = R.shape[0]
    else:
        raise RuntimeError, \
            "R array must be (3, 3) or (n, 3, 3); input has dimension %d" \
            % (rdim)
    
    T = zeros((nrot, 6, 6), dtype='float64')
    
    T[:, 0, 0] = R[:, 0, 0]**2    
    T[:, 0, 1] = R[:, 0, 1]**2
    T[:, 0, 2] = R[:, 0, 2]**2
    T[:, 0, 3] = sqr2 * R[:, 0, 1] * R[:, 0, 2]
    T[:, 0, 4] = sqr2 * R[:, 0, 0] * R[:, 0, 2]
    T[:, 0, 5] = sqr2 * R[:, 0, 0] * R[:, 0, 1]
    T[:, 1, 0] = R[:, 1, 0]**2
    T[:, 1, 1] = R[:, 1, 1]**2
    T[:, 1, 2] = R[:, 1, 2]**2
    T[:, 1, 3] = sqr2 * R[:, 1, 1] * R[:, 1, 2]
    T[:, 1, 4] = sqr2 * R[:, 1, 0] * R[:, 1, 2]
    T[:, 1, 5] = sqr2 * R[:, 1, 0] * R[:, 1, 1]
    T[:, 2, 0] = R[:, 2, 0]**2
    T[:, 2, 1] = R[:, 2, 1]**2
    T[:, 2, 2] = R[:, 2, 2]**2
    T[:, 2, 3] = sqr2 * R[:, 2, 1] * R[:, 2, 2]
    T[:, 2, 4] = sqr2 * R[:, 2, 0] * R[:, 2, 2]
    T[:, 2, 5] = sqr2 * R[:, 2, 0] * R[:, 2, 1]
    T[:, 3, 0] = sqr2 * R[:, 1, 0] * R[:, 2, 0]
    T[:, 3, 1] = sqr2 * R[:, 1, 1] * R[:, 2, 1]
    T[:, 3, 2] = sqr2 * R[:, 1, 2] * R[:, 2, 2]
    T[:, 3, 3] = R[:, 1, 2] * R[:, 2, 1] + R[:, 1, 1] * R[:, 2, 2]
    T[:, 3, 4] = R[:, 1, 2] * R[:, 2, 0] + R[:, 1, 0] * R[:, 2, 2]
    T[:, 3, 5] = R[:, 1, 1] * R[:, 2, 0] + R[:, 1, 0] * R[:, 2, 1]
    T[:, 4, 0] = sqr2 * R[:, 0, 0] * R[:, 2, 0]
    T[:, 4, 1] = sqr2 * R[:, 0, 1] * R[:, 2, 1]
    T[:, 4, 2] = sqr2 * R[:, 0, 2] * R[:, 2, 2]
    T[:, 4, 3] = R[:, 0, 2] * R[:, 2, 1] + R[:, 0, 1] * R[:, 2, 2]
    T[:, 4, 4] = R[:, 0, 2] * R[:, 2, 0] + R[:, 0, 0] * R[:, 2, 2]
    T[:, 4, 5] = R[:, 0, 1] * R[:, 2, 0] + R[:, 0, 0] * R[:, 2, 1]
    T[:, 5, 0] = sqr2 * R[:, 0, 0] * R[:, 1, 0]
    T[:, 5, 1] = sqr2 * R[:, 0, 1] * R[:, 1, 1]
    T[:, 5, 2] = sqr2 * R[:, 0, 2] * R[:, 1, 2]
    T[:, 5, 3] = R[:, 0, 2] * R[:, 1, 1] + R[:, 0, 1] * R[:, 1, 2]
    T[:, 5, 4] = R[:, 0, 0] * R[:, 1, 2] + R[:, 0, 2] * R[:, 1, 0]
    T[:, 5, 5] = R[:, 0, 1] * R[:, 1, 0] + R[:, 0, 0] * R[:, 1, 1]
    
    if nrot == 1:
        T = T.squeeze()

    return T

def nrmlProjOfVecMV(vec):
    """
    Gives vstacked p x 6 array to To perform n' * A * n as [N]*{A} for
    p hstacked input 3-vectors using the Mandel-Voigt convention.

    Nvec = normalProjectionOfMV(vec)

    *) the input vector array need not be normalized; it is performed in place
    
    """
    # normalize in place... col vectors!
    n = unitVector(vec)
    
    nmat = array([n[0, :]**2, 
                  n[1, :]**2, 
                  n[2, :]**2, 
                  sqr2 * n[1, :] * n[2, :], 
                  sqr2 * n[0, :] * n[2, :], 
                  sqr2 * n[0, :] * n[1, :]], 
                 dtype='float64')
    
    return nmat.T

def rankOneMatrix(vec1, *args):
    """
    Create rank one matrices (dyadics) from vectors.
      
      r1mat = rankOneMatrix(vec1)
      r1mat = rankOneMatrix(vec1, vec2)
    
      vec1 is m1 x n, an array of n hstacked m1 vectors
      vec2 is m2 x n, (optional) another array of n hstacked m2 vectors
    
      r1mat is n x m1 x m2, an array of n rank one matrices
                   formed as c1*c2' from columns c1 and c2
    
      With one argument, the second vector is taken to
      the same as the first.
    
      Notes:
    
      *)  This routine loops on the dimension m, assuming this 
          is much smaller than the number of points, n.
    """
    if len(vec1.shape) > 2:
        raise RuntimeError, "input vec1 is the wrong shape"
        
    if (len(args) == 0):
        vec2 = vec1.copy()
    else:
        vec2 = args[0]
        if len(vec1.shape) > 2:
            raise RuntimeError, "input vec2 is the wrong shape"
        
    m1, n1 = asmatrix(vec1).shape
    m2, n2 = asmatrix(vec2).shape
    
    if (n1 != n2):
        raise RuntimeError, "Number of vectors differ in arguments."
    
    m1m2 = m1 * m2
    
    r1mat = zeros((m1m2, n1), dtype='float64')
    
    mrange = asarray(range(m1), dtype='int')
    
    for i in range(m2):
        r1mat[mrange, :] = vec1 * tile(vec2[i, :], (m1, 1))
        mrange = mrange + m1
    
    r1mat = reshape(r1mat.T, (n1, m2, m1)).transpose(0, 2, 1)
    return squeeze(r1mat)

def skew(A):
    """
    skew-symmetric decomposition of n square (m, m) ndarrays.  Result
    is a (squeezed) (n, m, m) ndarray
    """
    if not isinstance(A, ndarray):
        raise RuntimeError, "input argument is of incorrect type; should be numpy ndarray."

    if A.ndim == 2:
        m = A.shape[0]
        n = A.shape[1]
        if m != n:
            raise RuntimeError, "this function only works for square arrays; " \
                  + "yours is (%d, %d)" %(m, n)
        A.resize(1, m, n)  
    elif A.ndim == 3:
        m = A.shape[1]
        n = A.shape[2]
        if m != n:
            raise RuntimeError, "this function only works for square arrays"
    else:
        raise RuntimeError, "this function only works for square arrays"
    
    return squeeze(0.5*(A - A.transpose(0, 2, 1)))
    
def symm(A):
    """
    symmetric decomposition of n square (m, m) ndarrays.  Result
    is a (squeezed) (n, m, m) ndarray.  
    """
    if not isinstance(A, ndarray):
        raise RuntimeError, "input argument is of incorrect type; should be numpy ndarray."

    if A.ndim == 2:
        m = A.shape[0]
        n = A.shape[1]
        if m != n:
            raise RuntimeError, "this function only works for square arrays; " \
                  + "yours is (%d, %d)" %(m, n)
        A.resize(1, m, n)  
    elif A.ndim == 3:
        m = A.shape[1]
        n = A.shape[2]
        if m != n:
            raise RuntimeError, "this function only works for square arrays"
    else:
        raise RuntimeError, "this function only works for square arrays"
    
    return squeeze(0.5*(A + A.transpose(0, 2, 1)))

def skewMatrixOfVector(w):
    """
    skewMatrixOfVector(w)

    given a (3, n) ndarray, w,  of n hstacked axial vectors, computes
    the associated skew matrices and stores them in an (n, 3, 3)
    ndarray.  Result is (3, 3) for w.shape = (3, 1) or (3, ).

    See also: vectorOfSkewMatrix
    """
    dims = w.ndim
    stackdim = 0
    if dims == 1:
        if len(w) != 3:
            raise RuntimeError, 'input is not a 3-d vector'
        else:
            w = vstack(w)
            stackdim = 1
    elif dims == 2:
        if w.shape[0] != 3:
            raise RuntimeError, 'input is of incorrect shape; expecting shape[0] = 3'
        else:
            stackdim = w.shape[1]
    else:
        raise RuntimeError, 'input is incorrect shape; expecting ndim = 1 or 2'
    
    zs = zeros((1, stackdim), dtype='float64')
    W = vstack([ zs,
                -w[2, :],
                 w[1, :],
                 w[2, :],
                 zs,
                -w[0, :],
                -w[1, :],
                 w[0, :],
                 zs ])
    
    return squeeze(reshape(W.T, (stackdim, 3, 3)))
    
def vectorOfSkewMatrix(W):
    """
    vectorOfSkewMatrix(W)

    given an (n, 3, 3) or (3, 3) ndarray, W, of n stacked 3x3 skew
    matrices, computes the associated axial vector(s) and stores them
    in an (3, n) ndarray.  Result always has ndim = 2.

    See also: skewMatrixOfVector
    """
    stackdim = 0
    if W.ndim == 2:
        if W.shape[0] != 3 or W.shape[0] != 3:
            raise RuntimeErrorl, 'input is not (3, 3)'
        stackdim = 1
        W.resize(1, 3, 3)
    elif W.ndim == 3:
        if W.shape[1] != 3 or W.shape[2] != 3:
            raise RuntimeError, 'input is not (3, 3)'
        stackdim = W.shape[0]
    else:
        raise RuntimeError, 'input is incorrect shape; expecting (n, 3, 3)'
    
    w = zeros((3, stackdim), dtype='float64')
    for i in range(stackdim):
        w[:, i] = r_[-W[i, 1, 2], W[i, 0, 2], -W[i, 0, 1]]
        
    return w

def multMatArray(ma1, ma2):
    """
    multiply two 3-d arrays of 2-d matrices
    """
    shp1 = ma1.shape
    shp2 = ma2.shape
    
    if len(shp1) != 3 or len(shp2) != 3:
        raise RuntimeError, 'input is incorrect shape; ' \
              + 'expecting len(ma1).shape = len(ma2).shape = 3'

    if shp1[0] != shp2[0]:
        raise RuntimeError, 'mismatch on number of matrices'

    if shp1[2] != shp2[1]:
        raise RuntimeError, 'mismatch on internal matrix dimensions'
    
    prod = zeros((shp1[0], shp1[1], shp2[2]))
    for j in range(shp1[0]):
        prod[j, :, :] = dot( ma1[j, :, :], ma2[j, :, :] )

    return prod

def uniqueVectors(v, tol=1.0e-12):
    """
    Sort vectors and discard duplicates.

      USAGE:
    
          uvec = uniqueVectors(vec, tol=1.0e-12)

    v   -- 
    tol -- (optional) comparison tolerance

    D. E. Boyce 2010-03-18
    """
        
    vdims = v.shape
    
    iv    = zeros(vdims)
    iv2   = zeros(vdims, dtype="bool")
    bsum  = zeros((vdims[1], ), dtype="bool")
    for row in range(vdims[0]):
        tmpord = num.argsort(v[row, :]).tolist()
        tmpsrt = v[ix_([row], tmpord)].squeeze()
        tmpcmp = abs(tmpsrt[1:] - tmpsrt[0:-1])
        indep  = num.hstack([True, tmpcmp > tol]) # independent values 
        rowint = indep.cumsum()
        iv[ix_([row], tmpord)] = rowint
        pass
    
    #
    #  Dictionary sort from bottom up
    #
    iNum = num.lexsort(iv)
    ivSrt = iv[:, iNum]
    vSrt = v[:, iNum]
    
    ivInd = zeros(vdims[1], dtype='int')
    nUniq = 1; ivInd[0] = 0
    for col in range(1, vdims[1]):
        if any(ivSrt[:, col] != ivSrt[:, col -1]):
            ivInd[nUniq] = col
            nUniq += 1
            pass
        pass
     
    return vSrt[:, ivInd[0:nUniq]]

def findDuplicateVectors(vec, tol=vTol, equivPM=False):
    """
    Find vectors in an array that are equivalent to within
    a specified tolerance
      
      USAGE:
    
          eqv = DuplicateVectors(vec, *tol)
    
      INPUT:
    
          1) vec is n x m, a double array of m horizontally concatenated
                           n-dimensional vectors.
         *2) tol is 1 x 1, a scalar tolerance.  If not specified, the default
                           tolerance is 1e-14.
         *3) set equivPM to True if vec and -vec are to be treated as equivalent
    
      OUTPUT:
    
          1) eqv is 1 x p, a list of p equivalence relationships.
    
      NOTES:
    
          Each equivalence relationship is a 1 x q vector of indices that
          represent the locations of duplicate columns/entries in the array
          vec.  For example:
    
                | 1     2     2     2     1     2     7 |
          vec = |                                       |
                | 2     3     5     3     2     3     3 |
    
          eqv = [[1x2 double]    [1x3 double]], where
    
          eqv[0] = [0  4]
          eqv[1] = [1  3  5]
    """
        
    vlen  = vec.shape[1]
    vlen0 = vlen
    orid  = asarray(range(vlen), dtype="int")

    torid = orid.copy()
    tvec  = vec.copy()
    
    eqv    = []
    eqvTot = 0
    uid    = 0
    
    ii = 1
    while vlen > 1 and ii < vlen0:
        dupl = tile(tvec[:, 0], (vlen, 1))
        
        if not equivPM:
            diff  = abs(tvec - dupl.T).sum(0)
            match = abs(diff[1:]) <= tol    # logical to find duplicates
        else:
            diffn  = abs(tvec - dupl.T).sum(0)
            matchn = abs(diffn[1:]) <= tol
            diffp  = abs(tvec + dupl.T).sum(0)
            matchp = abs(diffp[1:]) <= tol
            match = matchn + matchp
    
        kick = hstack([True, match])    # pick self too
        
        if kick.sum() > 1:
            eqv    += [torid[kick].tolist()]
            eqvTot  = hstack( [ eqvTot, torid[kick] ] )
            uid     = hstack( [ uid, torid[kick][0] ] )
        
        cmask       = ones((vlen,))
        cmask[kick] = 0
        cmask       = cmask != 0 
    
        tvec  = tvec[:, cmask]

        torid = torid[cmask]
        
        vlen = tvec.shape[1]

        ii += 1 

    if len(eqv) == 0:
        eqvTot = []
        uid    = []
    else:
        eqvTot = eqvTot[1:].tolist()
        uid    = uid[1:].tolist()
        
    # find all single-instance vectors
    singles = sort( setxor1d( eqvTot, range(vlen0) ) )

    # now construct list of unique vector column indices
    uid = int_( sort( union1d( uid, singles ) ) ).tolist()
    # make sure is a 1D list
    if not hasattr(uid,'__len__'):
        uid = [uid]
    
    return eqv, uid

def normvec(v):
    #mag = math.sqrt(num.sum(v[:]*v[:]))
    mag = num.linalg.norm(v)
    return mag

def normvec3(v):
    mag = math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
    return mag

def normalized(v):
    mag = normvec(v) # normvec3(v)
    n = v / mag
    return n

def cross(v1,v2):
    # return the cross product of v1 with another vector
    # return a vector
    newv3 = zeros(3,dtype='float64')
    newv3[0] = v1[1]*v2[2]-v1[2]*v2[1]
    newv3[1] = v1[2]*v2[0]-v1[0]*v2[2]
    newv3[2] = v1[0]*v2[1]-v1[1]*v2[0]
    return newv3

def determinant3(mat):
    v = cross(mat[0,:],mat[1,:])
    det = sum(mat[2,:] * v[:])
    return det


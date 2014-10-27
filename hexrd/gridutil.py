from numpy        import array, c_, r_, hstack, vstack, tile, \
                         cross, dot, diff, \
                         ceil, floor, \
                         ones, zeros
from numpy        import sum as asum
from numpy.linalg import det
import numpy as np
from hexrd.config import USE_NUMBA
if USE_NUMBA:
    import numba

def cellIndices(edges, points_1d):
    """
    get indices in a 1-d regular grid.

    edges are just that:

    point:            x (2.5)
                      |
    edges:   |1    |2    |3    |4    |5
             -------------------------
    indices: |  0  |  1  |  2  |  3  |
             -------------------------

    above the deltas are + and the index for the point is 1

    point:                  x (2.5)
                            |
    edges:   |5    |4    |3    |2    |1
             -------------------------
    indices: |  0  |  1  |  2  |  3  |
             -------------------------

    here the deltas are - and the index for the point is 2

    * can handle grids with +/- deltas
    * be careful when using with a cyclical angular array!  edges and points
      must be mapped to the same branch cut, and
      abs(edges[0] - edges[-1]) = 2*pi
    """
    ztol = 1e-12

    assert len(edges) >= 2, "must have at least 2 edges"

    points_1d = r_[points_1d].flatten()
    delta     = float(edges[1] - edges[0])

    on_last_rhs = points_1d >= edges[-1] - ztol
    points_1d[on_last_rhs] = points_1d[on_last_rhs] - ztol

    if delta > 0:
        on_last_rhs = points_1d >= edges[-1] - ztol
        points_1d[on_last_rhs] = points_1d[on_last_rhs] - ztol
        idx = floor( (points_1d - edges[0]) / delta )
    elif delta < 0:
        on_last_rhs = points_1d <= edges[-1] + ztol
        points_1d[on_last_rhs] = points_1d[on_last_rhs] + ztol
        idx = ceil( (points_1d - edges[0]) / delta ) - 1
    else:
        raise RuntimeError, "edges array gives delta of 0"
    # ...will catch exceptions elsewhere... 
    # if np.any(np.logical_or(idx < 0, idx > len(edges) - 1)):
    #     raise RuntimeWarning, "some input points are outside the grid"
    return array(idx, dtype=int)


def _fill_connectivity(out, m, n, p):
    i_con = 0
    for k in range(p):
        for j in range(m):
            for i in range(n):
                extra = k*(n+1)*(m+1)
                out[i_con, 0] = i + j*(n + 1) + 1 + extra
                out[i_con, 1] = i + j*(n + 1) + extra
                out[i_con, 2] = i + j + n*(j+1) + 1 + extra
                out[i_con, 3] = i + j + n*(j+1) + 2 + extra
                i_con += 1

if USE_NUMBA:
    _fill_connectivity = numba.njit(_fill_connectivity)

def cellConnectivity(m, n, p=1, origin='ul'):
    """
    p x m x n (layers x rows x cols)

    origin can be upper left -- 'ul' <default> or lower left -- 'll'

    choice will affect handedness (cw or ccw)
    """
    nele = p*m*n
    con  = np.empty((nele, 4), dtype=int)

    _fill_connectivity(con, m, n, p)

    if p > 1:
        nele = m*n*(p-1)
        tmp_con3 = con.reshape(p, m*n, 4)
        hex_con = []
        for layer in range(p - 1):
            hex_con.append(hstack([tmp_con3[layer], tmp_con3[layer + 1]]))
        con = vstack(hex_con)
        pass
    if origin.lower().strip() == 'll':
        con = con[:, ::-1]
    return con


if USE_NUMBA:
    @numba.jit # relies on loop extraction
    def cellCentroids(crd, con):
        nele, conn_count = con.shape
        dim = crd.shape[1]
        out = np.empty((nele, dim))
        inv_conn = 1.0/conn_count
        for i in range(nele):
            for j in range(dim):
                acc = 0.0
                for k in range(conn_count):
                    acc += crd[con[i,k], j]
                out[i,j] = acc * inv_conn
        return out
else:
    def cellCentroids(crd, con):
        """
        con.shape = (nele, 4)
        crd.shape = (ncrd, 2)

        con.shape = (nele, 8)
        crd.shape = (ncrd, 3)
        """
        nele = con.shape[0]
        dim  = crd.shape[1]
        centroid_xy = zeros((nele, dim))
        for i in range(len(con)):
            el_crds = crd[con[i, :], :] # (4, 2)
            centroid_xy[i, :] = (el_crds).mean(axis=0)
        return centroid_xy


def computeArea(polygon):
    """
    must be ORDERED and CONVEX!
    """
    n_vertices = len(polygon)
    polygon = array(polygon)
    
    triv = array([ [ [0, i-1], [0, i] ] for i in range(2, n_vertices) ])
             
    area = 0
    for i in range(len(triv)):
        tvp = diff( hstack([ polygon[triv[i][0], :], 
                             polygon[triv[i][1], :] ]), axis=0).flatten()
        area += 0.5 * cross(tvp[:2], tvp[2:])
    return area

def computeIntersection(line1, line2):
    """
    compute intersection of two-dimensional line intersection

    this is an implementation of two lines:

    line1 = [ [x0, y0], [x1, y1] ]
    line1 = [ [x3, y3], [x4, y4] ]

    
     <http://en.wikipedia.org/wiki/Line-line_intersection>
    """
    intersection = zeros(2)

    l1 = array(line1)
    l2 = array(line2)
     
    det_l1 = det(l1)
    det_l2 = det(l2)

    det_l1_x = det( vstack([ l1[:, 0], ones(2) ]).T )
    det_l1_y = det( vstack([ l1[:, 1], ones(2) ]).T )

    det_l2_x = det( vstack([ l2[:, 0], ones(2) ]).T )
    det_l2_y = det( vstack([ l2[:, 1], ones(2) ]).T )

    denominator = det( vstack([ [det_l1_x, det_l1_y], [det_l2_x, det_l2_y] ]) )

    if denominator == 0:
        intersection = []
    else:
        intersection[0] = det( vstack([ [det_l1, det_l1_x], [det_l2, det_l2_x] ]) ) / denominator
        intersection[1] = det( vstack([ [det_l1, det_l1_y], [det_l2, det_l2_y] ]) ) / denominator
    return intersection

def isinside(point, boundary, ccw=True):
    """
    Assumes CCW boundary ordering
    """
    pointPositionVector = hstack([         point - boundary[0, :], 0.])
    boundaryVector      = hstack([boundary[1, :] - boundary[0, :], 0.])
    
    crossVector = cross(pointPositionVector, boundaryVector)

    inside = False
    if crossVector[2] > 0:
        if ccw:
            inside = True
    elif crossVector[2] < 0:
        if not ccw:
            inside = True
    else:
        inside = True
    
    return inside

def sutherlandHodgman(subjectPolygon, clipPolygon):
    """
    """
    subjectPolygon = array(subjectPolygon)
    clipPolygon    = array(clipPolygon)
    
    numClipEdges = len(clipPolygon)

    prev_clipVertex = clipPolygon[-1, :]
 
    # loop over clipping edges
    outputList = array(subjectPolygon)
    for iClip in range(numClipEdges):
         
        curr_clipVertex = clipPolygon[iClip, :]

        clipBoundary = vstack([ curr_clipVertex, 
                                prev_clipVertex ])
 
        inputList  = array(outputList)
        if len(inputList) > 0:
            prev_subjectVertex = inputList[-1, :]            

        outputList = []
        
        for iInput in range(len(inputList)):

            curr_subjectVertex = inputList[iInput, :]
            
            if isinside(curr_subjectVertex, clipBoundary):
                if not isinside(prev_subjectVertex, clipBoundary):
                    subjectLineSegment = vstack([ curr_subjectVertex,
                                                  prev_subjectVertex ])
                    outputList.append(computeIntersection(subjectLineSegment, clipBoundary))
                    pass
                outputList.append(curr_subjectVertex)
            elif isinside(prev_subjectVertex, clipBoundary):
                subjectLineSegment = vstack([ curr_subjectVertex,
                                              prev_subjectVertex ])
                outputList.append(computeIntersection(subjectLineSegment, clipBoundary))
                pass
            prev_subjectVertex = curr_subjectVertex
            prev_clipVertex = curr_clipVertex
            pass            
        pass
    return outputList


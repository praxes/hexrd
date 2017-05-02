
# Some quaternion operations

from numba import guvectorize

@guvectorize(["(float64[:], float64[:])"], '(n)->(n)')
def conj(op, out):
    "conjugate of a quaterion, assumes (w, i, j, k)"
    out[0] = op[0]
    out[1:] = -op[1:]


@guvectorize(["(float64[:], float64[:], float64[:])"], '(n),(n)->(n)')
def prod(lho, rho, out):
    "product of two quaternions, assumes n==4"
    # based on http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/arithmetic/
    
    a, b, c, d = lho
    e, f, g, h = rho

    out[0] = a*e - b*f - c*g - d*h # w
    out[1] = b*e + a*f + c*h - d*g # i
    out[2] = a*g - b*h + c*e + d*f # j
    out[3] = a*h + b*g - c*f + d*e # k


@guvectorize(["(float64[:], float64[:], float64[:])"], '(n),(n)->(n)')
def quat_diff(lho, rho, out):
    out = prod(lho, conj(rho))



    

    

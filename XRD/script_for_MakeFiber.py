from grainIndex_MakeFiber import MakeFiber
import numpy
from numpy import dot



a1 = numpy.array([4,5,3])
R =  numpy.array([[ 0.95533649, -0.29552021 , 0.        ], [ 0.29552021 , 0.95533649,  0.        ], [ 0.   ,       0.   ,       1.        ]])
a2 = dot(R,a1)
f1 = MakeFiber(a1,a2)

b1 = numpy.array([-6,4,5])
b2 = dot(R,b1)
f2 = MakeFiber(b1,b2)
dist_,Closest_Rotation = f1.distBetweenFibers(f2)
print "original rotation", '\n',R
print "closest_rotation", '\n',Closest_Rotation


R2 = numpy.array([[ 0.69670671, -0.50724736,  0.50724736],       [ 0.50724736,  0.84835335,  0.15164665],       [-0.50724736,  0.15164665,  0.84835335]])

b2_prime = dot(R2,b1)
f2_prime = MakeFiber(b1,b2_prime)
dist_prime, Closest_Rotation_prime = f1.distBetweenFibers(f2_prime)
print "original rotation, different grain", '\n',R2
print "closest_rotation", '\n',Closest_Rotation_prime

#alternate method using rodrigues space
print "alternate method, rodrigues space distance, Closest Rotation"
dist_rodrigues, Closest_Rotation = f1.distBetweenFibers_Rodrigues(f2)
print dist_rodrigues, '\n', Closest_Rotation

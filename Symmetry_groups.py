import numpy
from math import pi,sin,cos
from Vector_funcs import *
from Diffraction_Analysis_10 import quaternion_map,Unit_Quaternion
def Sym_Group(schoenfliesTag):
    
    """ QUATOFLAUEGROUP - Generate quaternion representation for the specified """

    # Laue group.
    #
    # USAGE:
    #
    #      qsym = QuatOfLaueGroup(schoenfliesTag)
    #
    # INPUTS:
    #
    #      1) schoenfliesTag 1 x 1, a case-insensitive string representing the
    #      Schoenflies symbol for the desired Laue group.  The 11 available
    #      choices are:
    #
    #           Class           Symbol      n
    #          -------------------------------
    #           Triclinic       Ci (S2)     1
    #           Monoclinic      C2h         2
    #           Orthorhombic    D2h (Vh)    4
    #           Tetragonal      C4h         4
    #                           D4h         8
    #           Trigonal        C3i (S6)    3
    #                           D3d         6
    #           Hexagonal       C6h         6
    #                           D6h         12
    #           Cubic           Th          12
    #                           Oh          24
    #
    # OUTPUTS:
    #
    #      1) qsym is 4 x n, the quaterions associated with each element of the
    #      chosen symmetry group having n elements (dep. on group -- see INPUTS
    #      list above).
    #
    # NOTES:
    #
    #      *) The conventions used for assigning a RHON basis, {x1, x2, x3}, to
    #      each point group are consistent with those published in Appendix B
    #      of [1].
    #
    #      *) This routine replaces the older (individual) symmetry group
    #      routines, which are still available in the distribution, but made
    #      redundant on 12/1/2006.
    #
    # REFERENCES:
    #
    #      [1] Nye, J. F., "Physical Properties of Crystals: Their
    #      Representation by Tensors and Matrices", Oxford University Press,
    #      1985. ISBN 0198511655
    
    if schoenfliesTag == 'Ci' or schoenfliesTag == 'S2':
        # TRICLINIC
        angleAxis = numpy.array([0.0 ,  1 ,  0 ,  0])			# identity
        
    elif schoenfliesTag == 'C2h':
        # MONOCLINIC
        a1 = numpy.array([	0.0,   1,   0,   0])	# identity
        a2 = numpy.array([pi ,   0  , 1 ,  0])      	# twofold about 010 (x2)
        angleAxis = numpy.array([a1,a2])
    elif schoenfliesTag== 'D2h' or schoenfliesTag == 'Vh':
        # ORTHORHOMBIC
        
        a1 = numpy.array([0.0  ,    1 ,   0,    0])  		# identity
        a2 = numpy.array([pi     ,  1   , 0   , 0]) 		# twofold about 100
        a3 = numpy.array([pi      , 0,    1 ,   0]) 		# twofold about 010
        a4 = numpy.array([pi    ,   0 ,   0   , 1]) 		# twofold about 001
        angleAxis = numpy.array([a1,a2,a3,a4])
    elif schoenfliesTag== 'C4h':
        # TETRAGONAL (LOW)
        
        a1 = numpy.array([0.0  ,    1  ,  0  ,  0])  		# identity
        a2 = 	numpy.array([pi*0.5 ,  0 ,   0 ,   1])  		# fourfold about 001 (x3)
        a3 =numpy.array([pi    ,   0 ,   0  ,  1]) 		#
        a4 =numpy.array([	pi*1.5 ,  0 ,   0  ,  1]) 		#
        angleAxis = numpy.array([a1,a2,a3,a4])
    elif schoenfliesTag== 'D4h':
        # TETRAGONAL (HIGH)
        
        a1 = numpy.array([0.0   ,  1,   0,   0])  		# identity
        a2 = numpy.array([pi*0.5 , 0 ,  0 ,  1])  		# fourfold about  0  0  1 (x3)
        a3 = numpy.array([pi    ,  0  , 0  , 1]) 		#
        a4 = numpy.array([pi*1.5 , 0   ,0   ,1]) 		#
        a5 = numpy.array([pi      ,1,   0 ,  0])  		# twofold about  1  0  0 (x1)
        a6 = numpy.array([pi    ,  0 ,  1  , 0])  		# twofold about  0  1  0 (x2)
        a7 = numpy.array([pi     , 1  , 1   ,0])  		# twofold about  1  1  0
        a8 = numpy.array([pi     ,-1   ,1 ,  0]) 		# twofold about -1  1  0
        angleAxis = numpy.array([a1,a2,a3,a4,a5,a6,a7,a8])
    elif schoenfliesTag=='C3i' or schoenfliesTag == 'S6':
        # TRIGONAL (LOW)
        
        a1 =numpy.array([0.0 ,     1   , 0,    0])  		# identity
        a2 =numpy.array([ pi*2/3,   0 ,   0,    1])  		# threefold about 0001 (x3, c)
        a3 =numpy.array([ pi*4/3 ,  0,    0 ,   1]) 		#
        angleAxis = numpy.array([a1,a2,a3])
    elif schoenfliesTag == 'D3d':
        # TRIGONAL (HIGH)
        
        a1 = numpy.array([0.0,      1,    0,    0])  		# identity
        a2 = numpy.array([pi*2./3,   0 ,   0 ,   1])  		# threefold about 0001 (x3, c)
        a3 = numpy.array([pi*4./3 ,  0  ,  0  ,  1]) 		#
        a4 = numpy.array([pi     ,  1   , 0   , 0]) 		# twofold about  2 -1 -1  0 (x1, a1)
        a5 = numpy.array([pi   ,   -0.5 , sqrt(3)/2.,  0]) 		# twofold about -1  2 -1  0 (a2)
        a6 = numpy.array([pi    ,  -0.5, -sqrt(3)/2. , 0]) 		# twofold about -1 -1  2  0 (a3)
        angleAxis = numpy.array([a1,a2,a3,a4,a5,a6])
    elif schoenfliesTag =='C6h':
        # HEXAGONAL (LOW)
        
        a1 = numpy.array([0.0,      1   , 0,    0])  		# identity
        a2 = numpy.array([pi/3,     0   , 0 ,   1])  		# sixfold about 0001 (x3, c)
        a3 = numpy.array([pi*2./3,   0   , 0   , 1]) 
        a4 = numpy.array([ pi    ,   0  ,  0 ,   1]) 
        a5 = numpy.array([ pi*4./3 ,  0 ,   0  ,  1]) 
        a6 = numpy.array([ pi*5./3  , 0,    0   , 1]) 
        angleAxis = numpy.array([a1,a2,a3,a4,a5,a6])
    elif schoenfliesTag == 'D6h':
        # HEXAGONAL (HIGH)
        
        a1 = numpy.array([0.0 ,     1,    0   , 0])  		# identity
        a2 = numpy.array([pi/3 ,    0 ,   0   , 1])  		# sixfold about 0001 (x3, c)
        a3= numpy.array([pi*2./3 ,  0   , 0    ,1])
        a4 =numpy.array([ pi     ,  0   , 0   , 1]) 
        a5 =numpy.array([ pi*4./3  , 0,    0  ,  1]) 
        a6 =numpy.array([ pi*5./3  , 0 ,   0 ,   1]) 
        a7 =numpy.array([ pi,       1  ,  0 ,   0]) 		# twofold about  2 -1 -1  0 (x1, a1)
        a8 =numpy.array([ pi ,     -0.5 , sqrt(3)/2 , 0]) 		# twofold about -1  2 -1  0 (a2)
        a9 =numpy.array([ pi  ,    -0.5 ,-sqrt(3)/2 , 0]) 		# twofold about -1 -1  2  0 (a3)
        a10 =numpy.array([ pi  ,     sqrt(3)/2,  0.5,  0]) 		# twofold about  1  0 -1  0
        a11 = numpy.array([pi   ,    0   , 1   , 0]) 		# twofold about -1  1  0  0 (x2)
        a12 = numpy.array([pi    ,  -sqrt(3)/2 , 0.5,  0]) 		# twofold about  0 -1  1  0
        angleAxis = numpy.array([a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12])
    elif schoenfliesTag ==  'Th':
        # CUBIC (LOW)
        
        a1 = numpy.array([0.0  ,    1 ,   0,    0  		])# identity
        a2 = numpy.array([pi    ,   1 ,   0,    0  		])# twofold about  1  0  0 (x1)
        a3 = numpy.array([pi     ,  0 ,   1,    0  	])	# twofold about  0  1  0 (x2)
        a4 = numpy.array([pi      , 0 ,   0,    1  ])		# twofold about  0  0  1 (x3)
        a5 = numpy.array([pi*2/3,   1 ,   1,    1 ]) 	       # threefold about  1  1  1
        a6 = numpy.array([pi*4/3 ,  1 ,   1,    1 	])	#
        a7 = numpy.array([pi*2/3 , -1 ,   1,    1  ])		# threefold about -1  1  1
        a8 = numpy.array([pi*4/3 , -1 ,   1,    1 ])		#
        a9 = numpy.array([pi*2/3 , -1 ,  -1,    1 ]) 		# threefold about -1 -1  1
        a10 =numpy.array([ pi*4/3,  -1,   -1,    1]) 		#
        a11 =numpy.array([ pi*2/3,   1,   -1 ,   1])  		# threefold about  1 -1  1
        a12 =numpy.array([ pi*4/3,   1,   -1  ,  1]) 		#
        angleAxis = numpy.array([a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12])
    elif schoenfliesTag == 'Oh':
        # CUBIC (HIGH)
    
        a1=	numpy.array([0.0,      1,    0,    0])  		# identity
        a2 = numpy.array([pi*0.5,   1  ,  0   , 0])  		# fourfold about  1  0  0 (x1)
        a3 = numpy.array([pi    ,   1  ,  0   , 0]) 		#
        a4 = numpy.array([pi*1.5,   1  ,  0   , 0]) 		#
        a5 = numpy.array([pi*0.5,   0  ,  1   , 0])  		# fourfold about  0  1  0 (x2)
        a6 = numpy.array([pi    ,   0  ,  1   , 0]) 		#
        a7 = numpy.array([pi*1.5,   0  ,  1   , 0]) 		#
        a8 = numpy.array([pi*0.5,   0  ,  0   , 1])  		# fourfold about  0  0  1 (x3)
        a9 = numpy.array([pi    ,   0  ,  0   , 1]) 		#
        a10 = numpy.array([pi*1.5,   0 ,   0  ,  1]) 		#
        a11 = numpy.array([pi*2/3,   1 ,   1  ,  1])  		# threefold about  1  1  1
        a12 = numpy.array([pi*4/3,   1 ,   1  ,  1]) 		#
        a13 = numpy.array([pi*2/3,  -1 ,   1  ,  1])  		# threefold about -1  1  1
        a14 = numpy.array([pi*4/3,  -1 ,   1  ,  1]) 		#
        a15 = numpy.array([pi*2/3,  -1 ,  -1  ,  1])  		# threefold about -1 -1  1
        a16 = numpy.array([pi*4/3,  -1 ,  -1  ,  1]) 		#
        a17 = numpy.array([pi*2/3,   1 ,  -1  ,  1])  		# threefold about  1 -1  1
        a18 = numpy.array([pi*4/3,   1 ,  -1  ,  1]) 		#
        a19 = numpy.array([pi   ,    1 ,   1  ,  0])  		# twofold about  1  1  0
        a20 = numpy.array([pi    ,  -1 ,   1  ,  0]) 		# twofold about -1  1  0
        a21 = numpy.array([pi    ,    1,    0 ,   1]) 		# twofold about  1  0  1
        a22 = numpy.array([pi    ,   0,    1  ,  1]) 		# twofold about  0  1  1
        a23 = numpy.array([pi    ,  -1,    0  ,  1]) 		# twofold about -1  0  1
        a24 = numpy.array([pi    ,   0,   -1  ,  1]) 		# twofold about  0 -1  1
        angleAxis = numpy.array([a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24])
        
    else:
        error('unrecognized symmetry group.  See ''help QuatOfLaueGroup'' for a list of valid options.  Oh, and have a great day  -) you filthy animal')
    
        
        #
    tol = 1*10**-6
    rots = []
    rots2 = []
    for R in angleAxis:
        w,p1,p2,p3 = R
        #rots2.append(quaternion_map(Unit_Quaternion(R)))
        p = numpy.array([p1,p2,p3])
        q = Cross(p,e1func())
        if Mag(q)<tol:
            #print 'new'
            q = Cross(p,e2func())
        r = Cross(p,q)
        Rot = Rodrigues(w,p,q,r)
        Rot2 = Rodrigues3(Unit(p),w)
        rots2.append(Rot2)
        rots.append(Rot)
    return rots,rots2

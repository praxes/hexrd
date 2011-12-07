
from Vector_funcs import *
from Symmetry_groups import Sym_Group
from Diffraction_Analysis_10 import *
from data_class import inv_dict
from Vector_Data_Structures import *
from Misc_funcs import *
"""
def myfunc(alist,x):
    print alist
    print x
    print x**2
    return x**3
class myclass:
    def __init__(self,x):
        self.x = x
    def myfunc(self,arg1,arg2):
        return myfunc(arg1,arg2)
"""
hkls = []
rots,rots2 = Sym_Group('Oh')
f = open('./bcc.hkls')
for line in f:
    c = line.replace('\t',',')
    hkls.append(eval(c))
f.close()
angstrom = 10**-10
a = 2.83*angstrom
b = 2.83*angstrom
c = 2.83*angstrom
alpha = rad(90)
beta = rad(0)
gamma = rad(90)
lam = keV2Angstrom(80.725)*angstrom
g1 = a*e1func()
g2 = b*erfunc(alpha)
g3 = c*erhofunc(beta,gamma)
g_vecs = numpy.array([g1,g2,g3])
gij = Get_Metric(g1,g2,g3)
g_Vecs = dot(inv(gij),g_vecs)
g_1,g_2,g_3 = g_Vecs
orig_hkls = inv_dict()
hash_hkls = {}
hash_hkls_base = {}
only_hkls = {}
#picks out a complete set of recip vectors taking out symmetries and -h,-k,-l, by using a dictionary with Hash_3Vector's __hash__ function to id duplicates 
signs = [-1,1]
for sgn in range(2):
    sign_ = signs[sgn]
    for h,k,l in hkls:
        r0 = sign_*(h*g_1+k*g_2+l*g_3)        
        for rotation in rots:
            r0_rot = dot(rotation,r0)                        
            r0_hash = Hash_3Vector(r0_rot,tol =10**-4)
            hash_hkls[r0_hash]=(h,k,l)
#minimal set of recip vectors        
for h,k,l in hkls:
    r0 = h*g_1+k*g_2+l*g_3
    r0_ = Hash_3Vector(r0,tol =10**-4)
    hash_hkls_base[r0_] = (h,k,l)

#hash_hkls= { Hash_3Vector : (h,k,l)} -> Flip this around
tmp = inv_dict()
tmp.Absorb(hash_hkls)
tmp.Flip()
#tmp = {(hkl):(hash_3Vector1, hash_3Vector2,...)}
#orig_hkls = tmp
tmp2 = inv_dict()
for key in tmp.keys():
    hash_vecs = tmp[key]
    for g_vec in hash_vecs:
        tmp2.Add(key,g_vec.array)
orig_hkls = tmp2.Copy()

#orig_hkls = {(hkl):(numpy.array(g1_1),numpy.array(g1_2),...)
#orig_hkls is a complete set of non duplicate {hkl : vector list} which characterize the crystal
a_ = Vector_Data_Structure(dx = 2.5*.8*10**8,dy = 2.5*.8*10**8,dz = 2.5*.8*10**8,da1 = rad(1),da2 = rad(1),da3 = rad(1), fmap = Get_Pixel_Polar_Coords,inv_fmap = From_Pixel_Polar_Coords)
b_ = Vector_Data_Structure(dx = 2.5*.8*10**8,dy = 2.5*.8*10**8,dz = 2.5*.8*10**8,da1 = 2.5*.8*10**8,da2 = .1,da3 = .1, fmap = Get_Spherical_Coords,inv_fmap = From_Spherical_Coords) #b_.a1_a2_a3_dict: {(r,theta,phi):id}
a_.Define_kwargs(lam = lam)
c_ = Vector_Data_Structure(dx = 2.5*.8*10**8,dy = 2.5*.8*10**8,dz = 2.5*.8*10**8,da1 = 2.5*.8*10**8,da2 = .1,da3 = .1, fmap = Get_Spherical_Coords,inv_fmap = From_Spherical_Coords)
b_.Define_kwargs(e3 = e3func(),e1 = e1func(), e2 = e2func())
c_.Define_kwargs(e3 = e3func(),e1 = e1func(), e2 = e2func())


#base_gIs = {}
#base_gIs_r = inv_dict()
for gi in hash_hkls_base.keys():
    h,k,l = hash_hkls_base[gi]
    b_.Add(gi.array,(h,k,l))

#set up a deformation/rotation
p = Unit(e1func()+0.3*e2func())
q = Unit(Cross(p,e3func()))
r = Unit(Cross(p,q))
k1,k2,k3 = 1.1,1,1
#Distortion = k1*Bond(p,p)+k2*Bond(q,q)+k3*Bond(r,r)
Distortion = Ifunc()
rot_axis = Unit(p+.2*q - .5*r)
w0 =3.*math.pi/4.
Rotation = Rodrigues3(rot_axis,w0)
F = dot(Rotation,Distortion)#F = R
invF_T = inv(F.T) #invF.T = R
ct = 0
for key in orig_hkls.keys():
    #h,k,l = key
    gi_vectors = orig_hkls[key]
    if orig_hkls.triggerdict[key]==1:
        for gi in gi_vectors:
            a_.Add(gi,ct)            
            c_.Add(dot(invF_T,gi),ct)
            ct+=1
    else:
        a_.Add(gi_vectors,ct)
        c_.Add(dot(invF_T,gi_vectors),ct)
        ct+=1
#a_ represents all basis vectors for this symmetry group and lattice parameters - uses symmetry group on the h,k,l->(gI) vector and the symmetry group on the (-h,-k,-l)-> -gI vector and removes redundant results 
#b_ represents a minimal set of basis vectors - just the hkl's form of gI without applying any symmetries or inversions
#c_ represents diffraction data which would be obtained from a complete sweep of the sample.


#hkl = b_.raw_data.keys()[4]
hkl = (8,4,2)
gI = b_.raw_data[(8,4,2)]
b_.Query_Point_1D(gI, 1,.2*b_.dx) #returns (8,4,2)
gI = a_.raw_data[5] #arbitrary input vector
b_.Query_Point_1D(gI, axis = 1, search_radius = .2*b_.dx) #returns [(8,3,1),(7,5,0),(7,4,3)] multiple gIs with this magnitude
b_.Query_Point_1D(gI, axis = 1, search_radius = 10*b_.dx) #here I open up the search radius... returns [(6, 5, 3), (6, 6, 0), (8, 2, 2), (8, 3, 1), (7, 5, 0), (7, 4, 3), (6, 6, 2), (7, 5, 2)]


gI = orig_hkls[(6,5,5)][3] #arbitrary selection from (6,5,5) equivalents
rI = dot(Rotation,gI) #equivalently, = dot(invF_T,gI)
c_.Query_Point(rI,search_radius = c_.dx) #returns [1699]
rI_data = c_.raw_data[1699] #matches the input data to c_ (which is the same as rI since it was constructed in the same way)
trialR = Rotation
import time
tinit = time.time()
#out = Query_Orientation(orig_hkls,c_,trialR)
t_python = time.time() - tinit
mags = []
for key in b_.raw_data.keys():
    mags.append(Int_Map(Mag(b_.raw_data[key])/c_.dx))
mags.sort()
t_ = numpy.array(mags)
t_ = numpy.unique(t_)


"""
run:
    bash-3.2$ . setup_paths2.5.txt
    bash-3.2$ python2.5
    >>>
    
"""
#raise Exception, "are you ready?"
import getting_started1 as g
import pyublas

"""tinit = time.time()
out2 = g.Query_Orientation(orig_hkls.dict,orig_hkls.triggerdict,c_.xyzdict.dict, c_.xyzdict.triggerdict, c_.dx, c_.dx, trialR)
t_boost = time.time() - tinit

print "t_python, t_boost",t_python,t_boost
"""
#raise Exception
tinit = time.time()
#out = Find_Friedel_Pairs(c_)
tpython = time.time() - tinit
v1 = c_.raw_data[5] #array([  2.91741385e+10,   8.34807319e+08,   8.49381058e+09])
#friedel_ = out[5]
#v2 = c_.raw_data[friedel_] #array([ -2.91741385e+10,  -8.34807319e+08,  -8.49381058e+09])
#Mag(v1+v2) #close to zero
#tinit = time.time()
#out = g.Find_Friedel_Pairs(c_.raw_data.dict, c_.xyzdict.dict, c_.xyzdict.triggerdict, c_.dx, c_.dx)

tc = time.time() - tinit

print "tc, tpython",tc,tpython
#raise Exception

base_gIs_r = b_.a1_dict
a_rots = numpy.array(rots)
indexing_hkls = inv_dict()
r_keys = base_gIs_r.keys()
r_keys.sort()

#hkls_ = orig_hkls.keys()
#hkls_.sort()
max_index = 10
for key in r_keys[0:max_index]:
    hkl_ = base_gIs_r[key]
    if base_gIs_r.triggerdict[key]==1:
        for hkl__ in hkl_:
            indexing_hkls[hkl__] = orig_hkls[hkl__]
            indexing_hkls.triggerdict[hkl__] = orig_hkls.triggerdict[hkl__]
    else:
        indexing_hkls[hkl_] = orig_hkls[hkl_]
        indexing_hkls.triggerdict[hkl_] = orig_hkls.triggerdict[hkl_]

indexing_struct = c_.Blank_Copy()
ct = 0
for key in indexing_hkls.keys():
    gIs = indexing_hkls[key]
    if indexing_hkls.triggerdict[key]==1:
        for gI in gIs:
            indexing_struct.Add(gI,ct)
            ct+=1
    else:
        indexing_struct.Add(gI,ct)
        ct+=1
raise Exception
"""time trials"""
timeit(c_.Query_Point,numpy.array([4,5,43]))/timeit(g.Query_Point,c_.xyzdict.dict,c_.xyzdict.triggerdict, numpy.array([4,5,43]), c_.dx,c_.dx) #~120x speed up over python implementation
timeit(Query_Orientation,indexing_struct.raw_data,c_,Rotation)/timeit(g.Query_Orientation,indexing_struct.raw_data.dict,indexing_struct.raw_data.triggerdict,c_.xyzdict.dict,c_.xyzdict.triggerdict,c_.dx,c_.dx,Rotation) #~230X speed up over python implementation

#timeit(

out2 = g.Query_Orientation(indexing_hkls.dict,indexing_hkls.triggerdict,c_.xyzdict.dict, c_.xyzdict.triggerdict, c_.dx, c_.dx, trialR)

out = g.Indexer_Fibers_Query_Along_Fiber(indexing_hkls, base_gIs_r, c_.raw_data,c_.xyzdict, c_.dx,c_.dx,a_rots,2)
#g.Query_Point_1D(base_gIs_r.dict,base_gIs_r.triggerdict,orig_hkls[(7,5,0)][4],b_.dx,b_.dx)
#out = g.Make_Fibers_Nathan_1()
out = g.Make_Fibers_Nathan_1(b_.raw_data.dict,base_gIs_r.dict, base_gIs_r.triggerdict,c_.raw_data.dict,c_.dx,c_.dx,orig_hkls.dict)

out2 = g.Indexer_Dist_Fiber_Nathan_1(indexing_hkls.dict, base_gIs_r.dict, base_gIs_r.triggerdict, c_.raw_data.dict, c_.dx,c_.dx,a_rots,0,20,indexing_hkls.dict)
out = g.Make_Fibers_Nathan_1(indexing_hkls.dict,base_gIs_r.dict, base_gIs_r.triggerdict,c_.raw_data.dict,c_.dx,c_.dx,indexing_hkls.dict)

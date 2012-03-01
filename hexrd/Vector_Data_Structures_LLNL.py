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
from math import asin,acos,cos,sin

from hexrd.data_class import inv_dict
from hexrd.data_class import layered_data
from hexrd.Vector_funcs import * 
from hexrd.Misc_funcs import *

def Get_Cartesian(vector,b1,b2,b3):
    x = dot(vector,b1)
    y = dot(vector,b2)
    z = dot(vector,b3)
    return x,y,z
def From_Cartesian(x,y,z,b1,b2,b3):
    return x*b1+y*b2+z*b3
     
def CCW(x_in_rad):
    if x_in_rad > 2*math.pi:
        revs = Int_Map(x_in_rad/(2*math.pi),rounding = 'ROUND_DOWN')
        return CCW(x_in_rad - revs*2*math.pi)
    if x_in_rad<0:
        return 2*math.pi - abs(x_in_rad)
    else:
        return x_in_rad    
def Construct_s_s0_Polar_Coords(a,b):
    return e3func()*cos(b)+sin(b)*(cos(a)*e1func()+sin(a)*e2func())-e1func()
def mapOme(ome,omemin,omemax):
    assert omemin<=omemax, 'check angle map range, is omemin < omemax, e.g. -60,60?'
    #omemin = rad(omemin)
    #omemax = rad(omemax)
    ome = CCW(ome)
    mapped_ome = None
    if ome>CCW(omemin):
        mapped_ome = -1*(2*math.pi - ome)
    if ome<CCW(omemax):
        mapped_ome = ome
    
    
    return mapped_ome
        
    
def Get_LLNL_Angles(v, lam = 1., tol = 1e-3):
    my_e1 = -e3func()
    my_e2 = -e1func()
    my_e3 = e2func()
    R__= Bond(my_e1,e1func())+Bond(my_e2,e2func())+Bond(my_e3,e3func())
    my_v = dot(R__.T,v)
    sols = Get_Pixel_Polar_Coords(my_v, lam, tol)
    sol_eta_theta = inv_dict()
    for sol_num in sols.keys():
        w,a,b = sols[sol_num]
        eta,theta = alphabeta_polar_to_etatheta(a,b)        
        sol_eta_theta[sol_num] = [2*theta,eta,w]
    return sol_eta_theta
def alphabeta_polar_to_etatheta(a,b):
    theta = .5*acos(sin(b)*cos(a))
    s = erhofunc(b,a,p=e3func(),q=e1func(),r=e2func())
    Proj = Proj_func(e1func())
    s1 = Unit(dot(Proj,s))
    quadrant = Quadrant2D(s1,e1 = -e2func(),e2 = e3func())
    #print s1, quadrant
    if quadrant == 1 or quadrant == 2:
        eta = Compute_Angle(s1,-e2func())
    else:
        eta = math.pi + Compute_Angle(s1,e2func())
    return eta,theta

def From_Pixel_Polar_Coords(w,a,b,lam=1.):
    s_s0 = Construct_s_s0_Polar_Coords(a,b)
    R = Rodrigues3(e3func(),w)
    return dot(R.T,s_s0/lam)
def Get_Pixel_Polar_Coords(v,lam=1.,tol = 10**-3):
    #print 'trying',v,lam,tol
    sols = inv_dict()
    g1,g2,g3 = v
    d = Mag(v)**-1
    try:
        theta = asin(lam/(2*d))    
        beta = acos(lam*g3)
        alpha1 = acos(cos(2*theta)/sin(beta))
    except ValueError:
        return sols
    
    alpha2 = -alpha1
    alphas = [alpha1,alpha2]
    A = -g2
    B = g1
    if A>=0:
        phi = asin(B/sqrt(A**2+B**2))
    else:
        phi = math.pi - asin(B/sqrt(A**2+B**2))
    
    C = (cos(2*theta) - 1)/lam
    w1 = asin(C/(sqrt(A**2+B**2)))-phi
    #what quadrant is R.v
    w_max = math.pi/2. - phi
    w2 = w1 + 2*abs((w_max - w1))
    Rv = dot(Rodrigues3(e3func(),w1),v)
    quadrant = Quadrant2D(Rv,e1 = e1func(),e2 = e2func())
    w1 = CCW(w1)
    w2 = CCW(w2)
    if quadrant == 2:
        pair1 = [w1,alpha1,beta]
        pair2 = [w2,alpha2,beta]
    elif quadrant == 3:
        #alpha should be negative here
        pair1 = [w1,alpha2,beta]
        pair2 = [w2,alpha1,beta]
    tol = Mag(v*lam)*tol
    
    ct = 0
    for w,a,b in [pair1,pair2]:        
        g_c = From_Pixel_Polar_Coords(w,a,b,lam)
        if Mag(g_c-v)<tol:
            #raise Exception, "does not match"
            sols.Add(ct,(w,a,b))
            ct+=1
    if len(sols.keys())==1:
        print 'only one solution'
        sols.Add(sols.keys()[0]+1,sols[sols.keys()[0]])
    
    return sols #pair1,pair2
def Construct_s_s0_eta_theta(eta,theta):
    """cos(2*theta)*e1func() + sin(2*theta)*(cos(eta)*-e2func()+sin(eta)*e3func())-e1func()"""
    return erhofunc(2*theta,eta,p = e1func(),q = -e2func(),r = e3func())-e1func()
def From_Pixel_Polar_Coords_2(w,eta,theta,lam = 1.):
    R = Rodrigues3(e3func(),w)
    s_s0 = Construct_s_s0_eta_theta(eta,theta)
    return dot(R.T,s_s0/lam)
def Get_Pixel(beam_center_px,L,s,mm_per_px=.2):
    """uses convention for 0,0 in upper left corner, beam_center given in row,col format facing the detector """
    row,col = beam_center
    x_c = row*(-e3func())+col*(-e2func())+L*e1func()
    Hypotenuse = dot(s,e1func())/L
    pos = Hypotenuse*s
def Get_Recip_Vector_From_Pixel((w,prow,pcol),center_px,procession,L0,lam,tol = 10**-3,mm_per_px = .2):
    """add a graphic file here will show each vector definition"""
    R = Rodrigues3(e3func(),w)
    p = dot(R,procession)
    L_act = L0 + dot(p,-e1func())
    center = Pixel(center_px)
    pos = Pixel((prow,pcol))
    pos_x,pos_y = center.to_x_y(pos.pos,center.pos)
    pos_x = pos_x*mm_per_px
    pos_y = pos_y*mm_per_px
    R_ = pos_x*(-e2func())+pos_y*e3func()+L0*e1func()
    r = R_ - p
    s = Unit(r)
    s_s0 = s - e1func()
    v = dot(R.T,s_s0/lam)
    return v
def FromLLNL(v):
    my_e1 = -e3func()
    my_e2 = -e1func()
    my_e3 = e2func()
    R__= Bond(my_e1,e1func())+Bond(my_e2,e2func())+Bond(my_e3,e3func())
    my_v = dot(R__.T,v)
    return my_v
def Get_LLNL_Pixel_Sols_From_Recip_Vector(v,detectorGeom,precession,lam,tol = 10**-3,mm_per_px = .2):
    center_px = numpy.array([detectorGeom.xc,detectorGeom.yc])
    center_px = center_px/mm_per_px
    a1 = Pixel((0,0))
    cent = numpy.array(a1.to_row_col(a1.from_x_y(center_px)))
    my_v = FromLLNL(v)
    return Get_Pixel_Coords_From_Recip_Vector(my_v,cent,precession,detectorGeom.workDist,lam,tol,mm_per_px)
    
def Get_Pixel_Coords_From_Recip_Vector(v,center_px,procession,L0,lam=1.,tol = 10**-3,mm_per_px = .2):
    """must be given center_px in row,col array (0,0) in upper left facing the detector
    returns w,row,col"""
    sols = Get_Pixel_Polar_Coords(v,lam,tol)
    #pixel bases
    p1 = -e3func()
    p2 = -e2func()
    #detector_basis
    d1 = -e2func()
    d2 = e3func()
    #detector origin
    #center_px
    tmp_center = Pixel(center_px)
    tmp = Pixel((0,0))#just for function access later    
    final_pxs = inv_dict()
    for key in sols.keys():
        w,a,b = sols[key]
        #eta,theta = alphabeta_polar_to_etatheta(a,b)
        R = Rodrigues3(e3func(),w)
        G_pos = dot(R,procession)
        proc_x = dot(G_pos,e1func())
        L_act = L0 - proc_x
        s = Construct_s_s0_Polar_Coords(a,b)+e1func()
        #or, s = erhofunc(b,a,p=e3func(),q=e1func(),r=e2func()) = erhofunc(b,a)
        Hyp = L_act/dot(s,e1func())
        r_act = G_pos + Hyp*s
        Proj = Proj_func(e1func())
        rel_pixel_pos = r_act/mm_per_px 
        Pr = dot(Proj,rel_pixel_pos)
        x_rel_pos = dot(Pr,-e2func())
        y_rel_pos = dot(Pr,e3func())
        loc = tmp.from_x_y((x_rel_pos,y_rel_pos),origin = tmp_center.pos)
        prow,pcol = tmp.to_row_col(loc)
        final_pxs.Add(key,(w,prow,pcol))
                              
        
    return final_pxs
    
def Get_Pixel_From_Recip_Vector(v,center_px,procession,L0,lam=1.,tol=10**-3,mm_per_px = .2):
    sols = Get_Pixel_Polar_Coords(v,lam,tol)
  
    tmp_center = Pixel(center_px)
    tmp = Pixel((0,0))
    row0,col0 = center_px
    center_pos = (row0*(-e3func())+col0*(-e2func()))*mm_per_px + L0*e1func()
    final_pxs = inv_dict()
    for key in sols.keys():
        w,a,b = sols[key]
        #eta,theta = alphabeta_polar_to_etatheta(a,b)
        R = Rodrigues3(e3func(),w)
        G_pos = dot(R,procession)
        proc_x = dot(G_pos,e1func())
        L_act = L0 - proc_x
        s = Construct_s_s0_Polar_Coords(a,b)+e1func()
        #or, s = erhofunc(b,a,p=e3func(),q=e1func(),r=e2func()) = erhofunc(b,a)
        Hyp = L_act/dot(s,e1func())
        r_act = G_pos + Hyp*s
        Proj = Proj_func(e1func())
        rel_pixel_pos = r_act/mm_per_px 
        Pr = dot(Proj,rel_pixel_pos)
        x_rel_pos = dot(Pr,-e2func())
        y_rel_pos = dot(Pr,e3func())
        loc = tmp.from_x_y((x_rel_pos,y_rel_pos),origin = tmp_center.pos)
        prow,pcol = tmp.to_row_col(loc)
        final_pxs.Add(key,Pixel((prow,pcol)))
                              
        
    return final_pxs
        
        
    
    

    
    
def Get_Pixel_Polar_Coords_2(v,lam=1.,tol = 10**-3):
    g1,g2,g3 = v
    d = Mag(v)**-1
    theta = asin(lam/(2*d))    
    beta = acos(lam*g3)
    eta1 = asin(cos(beta)/sin(2*theta))
    difference = math.pi/2 - abs(eta1) #always positive
    eta2 = eta1 + 2*difference*numpy.sign(eta1)
    alpha1 = acos(cos(2*theta)/sin(beta))
    alpha2 = -alpha1
    alphas = [alpha1,alpha2]
    etas = [eta1,eta2]
    A = -g2
    B = g1
    if A>=0:
        phi = asin(B/sqrt(A**2+B**2))
    else:
        phi = math.pi - asin(B/sqrt(A**2+B**2))
    C = (cos(2*theta) - 1)/lam
    w1 = asin(C/(sqrt(A**2+B**2)))-phi
    #what quadrant is R.v
    w_max = math.pi/2. - phi
    w2 = w1 + 2*abs((w_max - w1))
    Rv = dot(Rodrigues3(e3func(),w1),v)
    quadrant = Quadrant2D(Rv,e1 = e1func(),e2 = e2func())
    w1 = CCW(w1)
    w2 = CCW(w2)
    if quadrant == 2:
        pair1 = [w1,eta2,theta]
        pair2 = [w2,eta1,theta]
    elif quadrant == 3:
        #alpha should be negative here
        pair1 = [w1,eta1,theta]
        pair2 = [w2,eta2,theta]
    tol = Mag(v*lam)*tol
    sols = inv_dict()
    ct = 0
    for w,e,t in [pair1,pair2]:        
        g_c = From_Pixel_Polar_Coords_2(w,e,t,lam)
        if Mag(g_c-v)<tol:
            #raise Exception, "does not match"
            sols.Add(ct,(w,e,t))
            ct+=1
    if len(sols.keys())==1:
        print 'only one solution'
        sols.Add(ct+1,sols[sols.keys()[0]])
    
    return sols#pair1,pair2

    
#need to add:
#Spherical Coords
#recip_vector to w,eta,2theta
#recip_vector to w,pix_x,pix_y


class Vector_Data:
    """
    Vector_Data(dx = ,dy = ,dz = , b1 = [e1],b2 = [e2],b3 = [e3], base_map = [Cartesian], inv_base_map = [From_Cartesian]
    """
    def __init__(self,dx,dy,dz,b1 = e1func(),b2 = e2func(),b3 = e3func(),base_map = Get_Cartesian, inv_base_map = From_Cartesian):
        self.base_map = base_map
        self.inv_base_map = inv_base_map
        self.raw_data = inv_dict()
        self.b1 = b1
        self.b2 = b2
        self.b3 = b3
        self.xyzdict = inv_dict()
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.raw_data = inv_dict()
        self.core_kwargs = {}
        self.core_kwargs['b1']=b1
        self.core_kwargs['b2']=b2
        self.core_kwargs['b3']=b3
        self.allow_unique = False
    def keys(self):
        return self.raw_data.keys()
    
    def __getitem__(self,key):
        return self.raw_data[key]
    def Query_Point_Spherical(self,input_vector,search_radius = 'default'):
        if search_radius == 'default':
            search_radius = self.dx
        
        hitlist = self.Query_Point(input_vector,search_radius)
        final_hits = []
        for pt in hitlist:
            found_vec = self.raw_data[pt]
            radial_distance = Mag(found_vec - input_vector)
            if radial_distance<=search_radius:
                final_hits.append(pt)
        return final_hits
    def Query_Point_2DSpherical_Coords(self,input_vector,search_params = 'default'):
        """search by specifying radial search dist (~2theta) and by angular distance (solid angle between input vector and possible hits (~eta,~w) """
        r = Mag(input_vector)
        if search_params == 'default':
            radial_dist = self.dx
            angular_dist = asin(self.dx/r)
            search_radius = 'default'
        else:
            radial_dist,angular_dist,search_radius = search_params
        hits = self.Query_Point(input_vector,search_radius)
        final_hits = []
        for hit in hits:
            vec = self.raw_data[hit]
            solid_angle = Compute_Angle(vec,input_vector)
            r_ = Mag(vec)
            if abs(r_-r)<radial_dist and abs(solid_angle - angular_dist)<angular_dist:
                final_hits.append(hit)
        return final_hits

        
        

    def Query_Point_Spherical_Coords(self,input_vector,angular_search_dist,search_radius = 'default'):
        """compute the solid angle between input vector and all pts meeting Query_Point raw, return indices for pts within the solid angle magnitude"""
        hits = Query_Point(input_vector,search_radius)
        r,theta,phi = Get_Spherical_Coords(input_vector)
        #sols_ = Get_Pixel_Polar_Coords_2(
        #r = Mag(input_vector)
        final_hits = []
        for hit in hits:
            vec = self.raw_data[hit]
            projected_vec = Unit(vec)*r
            angular_distance = Compute_Angle(input_vector,projected_vec)
            if angular_distance<angular_search_dist:
                final_hits.append(hit)
        return final_hits
    
            
            

    def Query_Point(self,input_vector,search_radius = 'default'):
        if search_radius == 'default':
            search_radius = self.dx
        x_seed,y_seed,z_seed = self.base_map(input_vector,**self.core_kwargs)
        dx,dy,dz = self.dx,self.dy,self.dz
        #ix = Int_Map(float(x)/dx)
        #iy = Int_Map(float(y)/dy)
        #iz = Int_Map(float(z)/dz)
        #x_seed, y_seed, z_seed = seed_pt        
        x_ = Int_Map(search_radius/self.dx)
        y_ = Int_Map(search_radius/self.dy)
        z_ = Int_Map(search_radius/self.dz)
        x_array = numpy.arange(2*x_+1)-x_
        y_array = numpy.arange(2*y_+1)-y_
        z_array = numpy.arange(2*z_+1)-z_
        keys = []
        keys2 = inv_dict()
        for i in range(len(x_array)):
            for j in range(len(y_array)):
                for k in range(len(z_array)):
                    x = x_seed + x_array[i]*self.dx
                    ix = Int_Map(x/self.dx)
                    y = y_seed+y_array[j]*self.dy
                    iy = Int_Map(y/self.dy)
                    z = z_seed + z_array[k]*self.dz
                    iz = Int_Map(z/self.dz)
                    key = (ix,iy,iz)
                    if self.xyzdict.has_key(key):
                        if self.xyzdict.triggerdict[key]==1:
                            #multiple data_ids at this query point, use extend
                            list_keys = list(self.xyzdict[key])
                            keys.extend(list_keys)
                            #for data_id in self.xyzdict[key]:
                            #    keys2.Add(data_id,key) #should only have one key per data_id
                            
                            
                        else:
                            #single value, use append
                            keys.append(self.xyzdict[key])
                            #data_id = self.xyzdict[key]
                            #keys2.Add(data_id,key)
                            
                    else:
                        pass
                    
        return keys
                         
        
        
    def Base_Reset(self):
        raw_data = self.raw_data.Copy()
        self.raw_data = inv_dict()
        self.xyzdict = inv_dict()
        for key in raw_data.keys():
            self.Add(key,raw_data[key])
    def has_key(self,key):
        return self.raw_data.has_key(key)
    def Add(self,vector,key,**kwargs):
        if len(kwargs.keys())==0:
            x,y,z = self.base_map(vector,b1 = self.b1,b2 = self.b2,b3 = self.b3)
        else:
            x,y,z = self.base_map(vector,**kwargs)
        if self.raw_data.has_key(key) and self.allow_unique == False:
            print "duplicate entry attempted",key
            return None;
        #raise Exception, "all points need unique ID (no appending)"
        self.raw_data.Add(key,vector)
        #alternate: self.raw_data[key]=vector
        #self.raw_data[key]=vector
        dx,dy,dz = self.dx,self.dy,self.dz
        ix = Int_Map(float(x)/dx)
        iy = Int_Map(float(y)/dy)
        iz = Int_Map(float(z)/dz)
        self.xyzdict.Add((ix,iy,iz),key)
    def ReMap_Base(self,map_func,**kwargs):
        """kwargs must contain dx = ,dy = ,dz = """
        dx = kwargs['dx'];        dy = kwargs['dy'];        dz = kwargs['dz']
        del kwargs['dx'];        del kwargs['dy'];        del kwargs['dz']
        self.dx = dx;        self.dy = dy;        self.dz = dz
        raw_data = self.raw_data.Copy()
        xyzdict = self.xyzdict.Copy()
        self.raw_data = inv_dict()
        self.xyzdict = inv_dict()        
        temp = inv_dict()
        for key in raw_data.keys():
            vector = raw_data[key]
            new_vec = map_func(vector,**kwargs)
            self.Add(new_vec,key)
            
            
"""
from Vector_Data_Structures import *
a = Vector_Data_Structure(dx = .1,dy = .1,dz = .1,da1 = .1,da2 = .1,da3 = .1, fmap = Get_Spherical_Coords,inv_fmap = From_Spherical_Coords)
"""
class Vector_Data_Structure(Vector_Data):
    """
    from Vector_Data_Structures import *
    a = Vector_Data_Structure(dx = .1,dy = .1,dz = .1,da1 = .1,da2 = .1,da3 = .1, fmap = Get_Spherical_Coords,inv_fmap = From_Spherical_Coords)
    """
    def __init__(self,dx,dy,dz,da1,da2,da3,fmap,inv_fmap,b1 = e1func(),b2 = e2func(), b3 = e3func(), base_map = Get_Cartesian, inv_base_map = From_Cartesian):
        Vector_Data.__init__(self,dx,dy,dz,b1 = b1, b2 = b2, b3 = b3, base_map = base_map, inv_base_map = inv_base_map)
        self.da1 = da1
        self.da2 = da2
        self.da3 = da3
        self.a1_dict = inv_dict()
        self.a2_dict = inv_dict()
        self.a3_dict   = inv_dict()
        self.a1_a2_dict = inv_dict()
        self.a2_a3_dict = inv_dict()
        self.a3_a1_dict = inv_dict()
        self.a1_a2_a3_dict = inv_dict()
        self.fmap = fmap
        self.inv_fmap = inv_fmap
        self.kwargs = {}
        self.allow_unique = True
    def Get_Params(self):
        """returns (dx,dy,dz,da1,da2,da3,fmap,inv_fmap,b1,b2,b3,base_map,inv_base_map)"""
        return (self.dx,self.dy,self.dz,self.da1,self.da2,self.da3,self.fmap,self.inv_fmap,self.kwargs)#,self.b1,self.b2,self.b3,self.base_map,self.inv_base_map)
    def Define_kwargs(self,**kwargs):
        #for key in kwargs.keys():
        #    print key,kwargs[key]
        self.kwargs = kwargs
    def Blank_Copy(self):
        dx,dy,dz,da1,da2,da3,fmap,inv_fmap,kwargs = self.Get_Params()
        tmp = Vector_Data_Structure(dx,dy,dz,da1,da2,da3,fmap,inv_fmap)
        tmp.kwargs = kwargs
        return tmp
    def Query_Point_1D(self,input_vector,axis = 1,search_radius = 'default'):
        if axis == 1:
            dx = self.da1
            dict_ = self.a1_dict
        if axis ==2:
            dx = self.da2
            dict_ = self.a2_dict
        if axis ==3:
            dx = self.da3
            dict_ = self.a3_dict
        if search_radius == 'default':
            search_radius = dx
        return self.Query_Point_1D_(dict_,input_vector, dx, search_radius)
    
    def Query_Point_1D_(self,search_dict, input_vector, dx, search_radius):
        """ could be made @staticmethod ?"""
        r0 = Mag(input_vector)
        r_ = Int_Map(search_radius/dx)
        a1_start = -r_
        a1_stop = r_+1
        step = 1
        a1_arr = numpy.arange(a1_start,a1_stop,step)
        ret_list = []
        for a1 in a1_arr:
            a1_ = r0 + a1*dx
            a1key = Int_Map(a1_/dx)
            #print a1key
            if search_dict.has_key(a1key):
                a1_vals = search_dict[a1key]
                if search_dict.triggerdict[a1key] == 1:
                    ret_list.extend(a1_vals)
                else:
                    ret_list.append(a1_vals)
                
        return ret_list
          
    def Add_Val(self,(a1,a2,a3),key):
        da1,da2,da3 = self.da1,self.da2,self.da3
        ia1 = Int_Map(float(a1)/da1)#,rounding = 'ROUND_DOWN')
        ia2 = Int_Map(float(a2)/da2)#,rounding = 'ROUND_DOWN')
        ia3 = Int_Map(float(a3)/da3)#,rounding = 'ROUND_DOWN')        
        self.a1_dict.Add(ia1,key)
        self.a2_dict.Add(ia2,key)
        self.a3_dict.Add(ia3,key)
        self.a1_a2_dict.Add((ia1,ia2),key)
        self.a2_a3_dict.Add((ia2,ia3),key)
        self.a3_a1_dict.Add((ia3,ia1),key)
        self.a1_a2_a3_dict.Add((ia1,ia2,ia3),key)
    def Map1(self,arg):
        return arg*self.da1
    def Map2(self,arg):
        return arg*self.da2
    def Map3(self,arg):
        return arg*self.da3
    def Refine(self,refine_factor,selector = 1):
        dx,dy,dz,da1,da2,da3,fmap,inv_fmap,kwargs = self.Get_Params()
        if selector == 1:
            dx,dy,dz,da1 = numpy.array([dx,dy,dz,da1])*refine_factor
        
        #if selector ==2:
        #    dx,dy,dz,da2 = numpy.array([dx,dy,dz,da1])*refine_factor
        return self.New(dx,dy,dz,da1,da2,da3,fmap,inv_fmap, **kwargs)
    
    def New(self,dx,dy,dz,da1,da2,da3,fmap,inv_fmap,**kwargs):
        #self.kwargs = kwargs
        #print 'kwargs',kwargs,self.kwargs
        tmp = Vector_Data_Structure(dx = dx,dy = dy,dz = dz, da1 = da1, da2 = da2, da3 = da3, fmap = fmap, inv_fmap = inv_fmap)
        tmp.Define_kwargs(**kwargs)
        raw_data = self.raw_data
        try:
            for key in raw_data.keys():
                tmp.Add(raw_data[key],key)
        except ValueError:
            print "error in New with:"
            print key,raw_data[key]
            raise ValueError
        return tmp
    def ReMesh(self,dx,dy,dz,da1,da2,da3,**kwargs):
        b = self.New(dx,dy,dz,da1,da2,da3,self.fmap,self.inv_fmap,**kwargs)
        self.da1 = b.da1;self.da2 = b.da2;self.da3 = b.da3;self.a1_dict = b.a1_dict
        self.a2_dict = b.a2_dict;self.a3_dict = b.a3_dict;self.a1_a2_dict = b.a1_a2_dict
        self.a2_a3_dict = b.a2_a3_dict;self.a3_a1_dict = b.a3_a1_dict;self.a1_a2_a3_dict = b.a1_a2_a3_dict
        
    def Add(self,ndata,key):
        """map the input vector to a1,a2,a3, then use Int Map to put in on integer grid, stor up different slices of the dictionary"""
        kwargs = self.kwargs
        if self.has_key(key):
            return -1
            #print "duplicate"
        out = self.fmap(ndata,**kwargs)
        if isinstance(out,inv_dict):
            for akey in out.keys():
                self.Add_Val(out[akey],key)
        else:
            self.Add_Val(out,key)
        Vector_Data.Add(self,ndata,key)
    def ReMap(self,map_func,**kwargs):
        """kwargs must contain dx = ,dy = ,dz = """
        da1 = kwargs['da1'];        da2 = kwargs['da2'];        da3 = kwargs['da3']
        del kwargs['da1'];        del kwargs['da2'];        del kwargs['da3']
        self.da1 = da1;        self.da2 = da2;        self.da3 = da3
        raw_data = self.raw_data    
        temp = inv_dict()
        kwargs = self.kwargs
        for key in raw_data.keys():
            vector = raw_data[key]
            new_vec = map_func(vector,**kwargs)
            self.Add(new_vec,key)
    def __repr__(self):
        return str(self)
        
    def __str__(self):
        ret = 'da1, \t da2, \t da3 \n'
        ret = ret+str(self.da1)+'\t'+str(self.da2)+'\t'+str(self.da3) + '\n'
        ret = ret + 'dx, \t dy, \t dz \n'
        ret = ret + str(self.dx)+'\t'+str(self.dy)+'\t'+str(self.dz) + '\n'
        ret = ret+self.fmap.__str__()+'\n'+self.inv_fmap.__str__()
        return ret
    #def Size(self):
        
    def Reset(self):
        kwargs = self.kwargs
        raw_data = self.raw_data
        self.a1_dict = inv_dict();self.a2_dict = inv_dict();self.a3_dict = inv_dict(); self.a1_a2_dict = inv_dict();self.a2_a3_dict = inv_dict(); self.a3_a1_dict = inv_dict(); self.a1_a2_a3_dict = inv_dict()
        for key in raw_data.keys():
            self.Add(raw_data[key],key,**kwargs)
        
    
    def GetSlice(self,(a1_min,a1_max),(a2_min,a2_max),(a3_min,a3_max),scales = 'default'):
        """slices the dictionary using numpy.nonzero(arr>arrmin) to find a subset of points which meet certain criteria
        @scales: is a tuple of mappings into datastructure coordinates. Usually, scales = (self.da1,self.da2,self.da3)
        
        b = a.New(a.dx,a.dy,a.dz,a.da1,a.da2,a.da3,Get_Pixel_Coords_From_Recip_Vector,Get_Pixel_From_Recip_Vector,center_px = center_px,procession = p0,L0=L0,lam = lam)
        a1_min,a1_max,a1_scale = (0,.75*2*math.pi,b.da1)
        a2_min,a2_max,a2_scale = (0,2048,b.da2)
        a3_min,a3_max,a3_scale = (0,2048,b.da3)
        returns data_structure of same nature as self with only the filtered points, also a dictionary with {data_id: w} where w is the index (0 or 1) at which fmap()[w] gives the correct map result. 
        data_struct,out =    b.GetSlice((a1_min,a1_max),(a2_min,a2_max),(a3_min,a3_max),scales = (a.da1,a.da2,a.da3))
        
        """
        if scales=='default':
            a1_scale,a2_scale,a3_scale = self.da1,self.da2,self.da3
        else:
            a1_scale,a2_scale,a3_scale = scales
        #map all limits into data structure coordinates
        a1_min = float(a1_min)/a1_scale
        a2_min = float(a2_min)/a2_scale
        a3_min = float(a3_min)/a3_scale
        a1_max = float(a1_max)/a1_scale
        a2_max = float(a2_max)/a2_scale
        a3_max = float(a3_max)/a3_scale
        

        a1a2a3keys = numpy.array(self.a1_a2_a3_dict.keys()) #so slicing works
        a1slice = a1a2a3keys[:,0] #column 1 key = (a1,a2,a3)
        a2slice = a1a2a3keys[:,1] #column 2 key = (a1,a2,a3)
        a3slice = a1a2a3keys[:,2] #column 3 key = (a1,a2,a3)
        #now find the values within the sliced range
        a1s = self.GetArray(a1slice,a1_min,a1_max);a1s = numpy.unique(a1s)
        a2s = self.GetArray(a2slice,a2_min,a2_max);a2s = numpy.unique(a2s)
        a3s = self.GetArray(a3slice,a3_min,a3_max);a3s = numpy.unique(a3s)
        
        tmp = self.a1_a2_a3_dict
        #tmp : {fmap(data_id) : data_id}
        tmp2 = tmp.Get_Flipped_Copy()    
        #tmp2 : {data_id : fmap(data_id)}
        out = inv_dict()
        ct=0
        for a1key in a1s:
            for a2key in a2s:
                for a3key in a3s:
                    key = (a1key,a2key,a3key) #form all possible a1,a2,a3 in the valid range
                    if tmp.has_key(key):
                        data_ids = tmp[key] #tmp_out is the unique identifier(s) (if trigger = 1, multiples)
                        if tmp.triggerdict[key]==1:
                            for data_id in data_ids:
                                #the following snip is the same for if triggerdict = 0
                                stor = []
                                if tmp2.triggerdict[data_id]==1: #if there are multiple associations to this vector (true for most) e.g. if a vector input has more than one fmap result
                                    for m in range(len(tmp2[data_id])): #usually range(2)
                                        a1,a2,a3 = tmp2[data_id][m]
                                        if (a1_min<=a1<=a1_max) and (a2_min<=a2<=a2_max) and (a3_min<=a3<=a3_max):
                                            stor.append(m)
                                        else:
                                            pass
                            #for a1,a2,a3 in tmp2[tmp_out]:
                                else:
                                    #just store up first index
                                    stor.append(0)
                            #do
                                out.Add(data_id,stor)
                                ct+=1
                                
                        else:
                        #if tmp.triggerdict[key]!=1:
                            stor = []
                            if tmp2.triggerdict[data_ids]==1: #if there are multiple associations to this vector (true for most)
                                for m in range(len(tmp2[data_ids])):
                                    a1,a2,a3 = tmp2[data_ids][m]
                                    if (a1_min<=a1<=a1_max) and (a2_min<=a2<=a2_max) and (a3_min<=a3<=a3_max):
                                        stor.append(m)
                                    else:
                                        pass
                            #for a1,a2,a3 in tmp2[tmp_out]:
                            else:
                                stor.append(0)
                            out.Add(data_ids,stor)
                            ct+=1
        
        dx,dy,dz,da1,da2,da3,fmap,inv_fmap = self.Get_Params()
        blank_struct = Vector_Data_Structure(dx,dy,dz,da1,da2,da3,fmap,inv_fmap)
        blank_struct.kwargs = self.kwargs
        tmp = self.raw_data 
        #tmp: {data_id : input_data}
        for data_id in out.keys():
            r_vec = tmp[data_id]
            blank_struct.Add(r_vec,data_id)
            
        return blank_struct,out

    def GetArray(self,arr,lower,upper):
        """internal use"""
        a1_max = arr<=upper
        a1_min = arr>=lower
        a1_isect = (a1_max==a1_min)
        indices = numpy.nonzero(a1_isect)
        return arr[indices]
        
            
class Pixel:
    """
    a = Pixel((1024,1000))
    x,y = a.to_x_y(a.pos) #relative to lower left corner (Default center @ Pixel((2048,0))
    center = Pixel((1024,1024))
    x,y = a.to_x_y(a.pos,center = center.pos) #relative to 'center'
    
    """
    def __init__(self,(row,column)):
        self.b1 = numpy.array([1,0])
        self.b2 = numpy.array([0,1])
        self.origin = numpy.array([0,0])
        self.xy_origin = self.from_row_col((2048,0))
        self.pos = self.from_row_col((row,column),self.origin)
    @staticmethod
    def testme(arg1,arg2):
        print arg1,arg2
    def to_matlab(self,arr,origin = 'default'):
        if origin=='default':
            origin = self.origin
        rel_pos = arr-origin
        x1 = dot(rel_pos,self.b1)
        x2 = dot(rel_pos,-self.b2)
        return x1,x2
    def from_matlab(self,(x1,x2),origin = 'default'):
        if origin=='default':
            origin = self.origin
        return x1*self.b1 + x2*(-self.b2) + origin
    
    def from_row_col(self,(row,col),origin = 'default'):
        if origin == 'default':
            origin = self.origin
        return row*(-self.b2)+col*self.b1 + origin
    def to_row_col(self,arr,origin = 'default'):
        if origin == 'default':
            origin = self.origin
        rel_pos = arr-origin
        row = dot(rel_pos,-self.b2)
        col = dot(rel_pos,self.b1)
        return row,col
    def to_x_y(self,arr,origin = 'default'):
        if origin == 'default':
            origin = self.xy_origin
        rel_pos = arr-origin
        x = dot(rel_pos,self.b1)
        y = dot(rel_pos,self.b2)
        return x,y
    def from_x_y(self,(x,y),origin = 'default'):
        if origin == 'default':
            origin = self.xy_origin
        return x*self.b1+y*self.b2+origin
    def __str__(self):
        astr = 'Pixel \n row \t col \t x \t y \n'
        row,col = self.to_row_col(self.pos)
        x,y = self.to_x_y(self.pos)
        astr = astr+str(row)+'\t'+str(col)+'\t'+str(x)+'\t'+str(y)+'\n'
        return astr
    def __repr__(self):
        return self.__str__()

#get friedel pairs, then try index these
def Extract_Friedel_Pairs(rIs, friedel_pairs,search_radius = 'default'):
    """rIs - Vector_Data_Structure w/ raw_data{id:vector}"""
    #friedel_pairs = Find_Friedel_Pairs(rIs)
    friedel_structure = rIs.Blank_Copy()
                    
    for key in friedel_pairs.keys():
        id1 = key;
        id2 = friedel_pairs[id1]
        if friedel_pairs.triggerdict[id1]==1:
            for id_ in id2:
                v = rIs.raw_data[id_]
                #if friedel_structure.has_key(id_)!=True:
                #    friedel_structure.Add(v,id_)
                friedel_structure.Add(v,id_)
        else:
            id_ = id2;
            v = rIs.raw_data[id_]
            #if friedel_structure.has_key(id_)!=True:
            #    friedel_structure.Add(v,id_)
            friedel_structure.Add(v,id_)
    return friedel_structure

        
    
        
        

def Find_Friedel_Pairs(rIs,search_radius = 'default'):
    friedel_pairs = inv_dict()
    for rI_index in rIs.raw_data.keys():
        rI = rIs.raw_data[rI_index]
        friedel_pair = -rI
        pairs_found = rIs.Query_Point(friedel_pair,search_radius)
        for pair_found in pairs_found:
            friedel_pairs.Add(rI_index,pair_found)
    #friedel_pairs:{dataid1 : fridelpair_of_dataid1}
    #friedel_pairs:{dataid2 : (fridelpair1_of_dataid2, friedelpair2_of_dataid2,...)
        
    return friedel_pairs

def Query_Orientation(gIs, rIs, Rtrial,search_radius = 'default'):
    """gIs: inv_dict
    rIs: Vector_Data_Structure
    """
    num_hkls = len(gIs.keys())
    slots = inv_dict()
    for key in gIs.keys():
        gI_vals = gIs[key]
        if gIs.triggerdict[key]==1:
            for gI in gI_vals:
                r_test = dot(Rtrial,gI)
                hits = rIs.Query_Point(r_test,search_radius)
                for hit in hits:
                    slots.Add(key,hit)
        else:
            r_test = dot(Rtrial,gI_vals)
            hits = rIs.Query_Point(r_test,search_radius)
            for hit in hits:
                slots.Add(key,hit)
    return slots

def Query_Orientation_Unique(gIs,rIs,Rtrial,search_radius = 'default'):
    """gIs: indexing_hkls_struct
    rIs: Vector_Data_Structure
    """    
    num_hkls = len(gIs.keys())
    slots = inv_dict()
    allids = {}
    ct = 0
    crystal_struct = layered_data({})    
    for key in gIs.keys():
        gI = gIs[key]
        r_test = dot(Rtrial,gI)
        hits = rIs.Query_Point(r_test,search_radius)
        if len(hits)>1:
            print "grain splitting"
            for hit in hits:
                crystal_struct.SplitLayer
        else:
            assoc = tuple(key,hits[0])
            crystal_struct.AddItemToAll((key,hits[0]))
            
        
    return ct
                
                
"""class layered_data:
    def __init__(self,data_elem):
        self.dataelement = data_elem
        self.data = {}
        self.index = 0
    def keys(self):
        return self.data.keys()
    def __getitem__(self,key):
        return self.data[key]
    def GetNewElement(self):
        return self.dataelement.copy()
    def AddItem(self,(key,value)):
        for structkey in self.keys():
            #self[structkey].Add(key,value)
            #or..
            self[structkey][key]=value
    def AddLayer(self):
        self.data[self.index] = self.GetNewElement()
        self.index+=1
    
        
"""

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
import decimal
import math
import time

def Int_Map(arg,rounding = 'ROUND_HALF_DOWN'):
    tmp = decimal.Decimal(arg.__str__())
    int_arg = int(tmp.to_integral(rounding = rounding))
    return int_arg
def deg(x):
    return x*180./math.pi
def rad(x):
    return x*math.pi/180.
def Find_List_Location(a_sorted_list,arg):
    tmp = a_sorted_list[:]
    if arg<tmp[0]:
        tmp.insert(arg,0)
        out = 0
    if arg>tmp[-1]:
        tmp.insert(arg,len(tmp))
        out = len(tmp)+1
    else:
        for i in range(len(tmp)-1):
            if tmp[i]<arg and tmp[i+1]>arg:
                
                tmp.insert(i+1,arg)
                out = i+1
                break
    return tmp[:],out


    

def list_map(alist,index):
    return alist[index]
def list_map_inv(alist,val):
    return alist.index(val)


def timeit(afunc,*args):
    #print args
    tinit = time.time()
    afunc(*args)
    tfinal = time.time()
    return tfinal-tinit
            
def CCW(x_in_rad):
    """convert angle to counterclockwise angle in radians - restricted to between 0 and 2pi
    e.g.
    CCW(3*pi) = pi
    CCW(-pi/2) = 3*pi/2
    CCW(1.5pi) = 1.5 pi
    """
    if x_in_rad > 2*math.pi:
        revs = Int_Map(x_in_rad/(2*math.pi),rounding = 'ROUND_DOWN')
        return CCW(x_in_rad - revs*2*math.pi)
    if x_in_rad<0:
        return 2*math.pi - abs(x_in_rad)
    else:
        return x_in_rad

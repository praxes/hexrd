import decimal
import math
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
    import time
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

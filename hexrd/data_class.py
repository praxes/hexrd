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
 
import numpy
    
            
class inv_dict:
    '''(extendable) dictionary contains two data structures:
1. a dictionary of keys,values (access in python via inv_dict.dict)
2. a dictionary of keys,trigger, where trigger is either 0 or 1 depending on if the dictionary was `extended` (1) or not (0). 
These two structures allow for extending of the dictionary and dictionary inversion (turn keys->values)
example:
tmp = inv_dict()
tmp.Add(4,6)
tmp.Add(4,7)
tmp.Print()
tmp.Flip()
tmp.Print()
import numpy
tmp.Add(5,numpy.array(4,5,6))
tmp.Add(5,numpy.array(3,4,5))
tmp.Print()
adict = {}
a[4]=4
a[5]=10
tmp = inv_dict()
tmp.Absorb(a)
tmp.Print()

'''
    def __init__(self):
        self.dict={}
        self.triggerdict={}
        self.status=1
    def __getitem__(self,key):
        return self.dict[key]
    def __setitem__(self,key,value):
        self.dict[key]=value
    def keys(self):
        return self.dict.keys()
    def values(self):
        return self.dict.values()
    def has_key(self,key):
        return self.dict.has_key(key)
    def Add(self,key,value):
        if type(key) == numpy.ndarray:
            hash_key = tuple(key)
            self.Add(hash_key,value)
        else:
            #if(hasattr(key,"__hash__")):
            if self.has_key(key):
                self.Extend(key,value)
            else:
                self[key]=value
                self.triggerdict[key]=0
        
    def Extend(self,key,value):
        value_ = tuple([value])
        if self.triggerdict[key]==1:
            self[key] = self[key]+value_
        else:
            self.triggerdict[key]=1
            old_val = tuple([self[key]])
            new_val = old_val+value_
            self[key]=new_val
    def Flip(self,returncopy = False):
        tmp = inv_dict()
        for key in self.keys():
            value = self[key]
            if self.triggerdict[key]==1:
                for val in value:
                    tmp.Add(val,key)
            else:
                tmp.Add(value,key)
        if returncopy==False:
            self.dict = tmp.dict
            self.triggerdict = tmp.triggerdict
            self.status = self.status*-1
            return 
        else:
            return tmp
    def Copy(self):
        tmp = inv_dict()
        for key in self.keys():
            tmp[key]=self[key]
            tmp.triggerdict[key]=self.triggerdict[key]
        return tmp
    def Absorb(self,adict):
        for key in adict.keys():
            self.Add(key,adict[key])
        
    def __str__(self):
        return self.Print()
    def __repr__(self):
        print 'Inv_dictionary'
        return self.__str__()
    def Print(self):
        printstring = ''
        printstring+='Keys, Values \n'  
        for key in self.keys():
            printstring += key.__str__() + '\t'
            printstring += self[key].__str__() + '\n'
        return printstring
        
          
    
        
        
    
   
class layered_data:
    """test = layered_data({})
    test.add_layer()
    test[3]=6
    test[5]=8
    values = [99,100]
    key = 9
    test.split_and_add(key,values)
    test
    test[14]=14.6
    test
    test.baseids #current values
    test.split_and_add(45,(45,19))
    """
    def __init__(self,data):
        self.data = {}
        self.data_element = data
        self.globalct = 0
        self.baseids = [0]
    def __getitem__(self,key):
        return self.data[key]
    def add_layer(self,data = 'empty_element'):
        if (data == 'empty_element'):
            data = self.get_new_element()
        self.data[self.globalct] = data
        self.globalct+=1
    def get_new_element(self):
        """blank type"""
        return self.data_element.__new__(type(self.data_element))
    def __repr__(self):
        astr = ''
        for key in self.keys():
            astr+=str(key)
            astr+="\n"
            astr+=self[key].__repr__()
            astr+="\n"
        return astr
    def keys(self):
        return self.data.keys()
    def __setitem__(self,key,value):
        for globalkey in self.keys():
            self.data[globalkey][key] = value
    def copy(self):
        """copy only base ids"""
        copies = {}
        for baseid in self.baseids:
            new_ = layered_data(self.data_element)
            base_data = self.data[baseid]
            new_.data = {0:base_data.copy()}
            copies[baseid] = new_
        return copies
    def split_and_add(self,key,values):
        num_splits = len(values)
        #want to make new layered data with the different copies
        new_baseids = []
        
        for value in values:
            new_ = self.copy()
            for baseid in self.baseids:
                (new_[baseid])[key]=value
                new_baseids.append(self.globalct)
                self.add_layer(new_[baseid])
                
        self.baseids = new_baseids

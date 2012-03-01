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
def write_lines(open_file,data):
    for elem in data:
        open_file.write(elem.__str__())
        open_file.write('\t')
    open_file.write('\n')
def count_lines(filename):
    f = open(filename,'r')
    ct = 0
    for line in f:
        ct+=1
    return ct
def writedata(filename, data, row_function = False):
    f = open(filename, 'w')
    ct = -1
    for elem in data:
        ct+=1
        if row_function:
            out = [ct,row_function(elem)]
            write_lines(f,out)
            continue
        
        if len(elem)==1:
            elem = [ct, elem]
        ct+=1
        write_lines(f,elem)
    print 'wrote',f.name
    f.close()
class fileDataClass:
    def __init__(self, name):
        self.name = name
        self.data = {}
    def addpoint(self,x,y):
        try:
            self.data[tuple(x)] = y
        except TypeError:
            self.data[tuple([x])] = y
    
        
    def writefile(self):
        keys = self.data.keys()
        keys.sort()
        f = open(self.name, 'w')
        for key in keys:
            ndata = []
            ndata.extend(key)
            data = self.data[key]
            try:
                ndata.extend(data)
            except TypeError:
                ndata.extend([data])
            write_lines(f,ndata)
        print 'wrote',f.name
        return

class gnuplotDataClass(fileDataClass):
    def __init__(self,fname):
        fileDataClass.__init__(self,fname)
        self.gplotname = fname + '.gpi'
        self.epsname = fname + '.eps'
        
    def writefile(self, gplot_command):
        fileDataClass.writefile(self)
        gplotfile = open(self.gplotname, 'w')
        gplotfile.write('plot "'+str(self.name)+'" '+gplot_command+'\n')
        gplotfile.close()
    def writeline(self, gplot_command):
        gplotfile = open(self.gplotname, 'a')
        gplotfile.write(gplot_command+'\n')
        gplotfile.close()
    def writeimage(self):
        gplotfile = open(self.gplotname, 'a')
        gplotfile.write('\n')
        gplotfile.write('set terminal postscript enhanced\n')
        gplotfile.write('set output "'+self.epsname+'"'+'\n')
        gplotfile.write('replot')
        gplotfile.close()
        import os
        os.system('gnuplot ' + gplotfile.name)
        
        
                        
        
    

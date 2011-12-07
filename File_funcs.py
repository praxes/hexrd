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
        
        
                        
        
    

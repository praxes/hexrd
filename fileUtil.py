# DO-NOT-DELETE revisionify.begin() 
#
#   Copyright (c) 2007-2009 Lawrence Livermore National Security,
#   LLC. Produced at the Lawrence Livermore National Laboratory (Nathan
#   Barton <barton22@llnl.gov>) CODE-OCEC-08-104.
#   
#   Please also read the file NOTICES.
#   
#   This file is part of the mdef package (version 0.2) and is
#   free software: you can redistribute it and/or modify it under the
#   terms of the GNU Lesser General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#   
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   Lesser General Public License for more details.
#   
#   A copy of the GNU Lesser General Public License may be found in the
#   file NOTICES. If this file is missing, see
#   <http://www.gnu.org/licenses/>.
#
# DO-NOT-DELETE revisionify.end()
import mdef

import re
import sys
import os
import time
import string
import time

# from math import *

class Log:
    """for logging"""
    def __init__(self,logFileName=None,toScreen=True):
        self.logFile  = None
        self.toScreen = toScreen
        if logFileName:
            self.logFile = open(logFileName, 'w')
            self('for more information, see log file: %s' % logFileName )
        return
    def __call__(self, toLog, suppressToScreen=False):
        if self.logFile:
            #if isinstance(toLog,str):
            #    self.logFile.write(toLog)
            #elif hasattr(toLog,'__repr__'):
            #    self.logFile.write(toLog.__repr__())
            #else:
            #    raise RuntimeError, 'do not know how to write '
            self.logFile.write('\n'+str(toLog))
            self.logFile.flush()
        if self.toScreen and not suppressToScreen:
            print >> sys.stdout, str(toLog)
    def close(self):
        if self.logFile: self.logFile.close()
        return
    
class FileDescr:
    """base class"""
    def __init__(self, *args):
        raise RuntimeError, "need implementation"
    def __call__(self, **args):
        raise RuntimeError, "need implementation"
    def __str__(self):
        return str(self.__dict__)
    def __repr__(self):
        'this is not really a good implementation of this method'
        return str(self.__class__) + "(" + str(self.__dict__) + ")"
    def filename(self):
        raise RuntimeError, "need implementation"
        
class FileLink(FileDescr):
    def __init__(self, fileName):
        self.fileName = fileName
        if os.path.isfile(os.path.join('/',fileName)):
            # assume is absolute path
            self.absFileName = os.path.join('/',fileName)
        else:
            # assume is relative path
            self.absFileName = os.path.join(os.getcwd(),fileName)
        return
    def __call__(self, **args):
        os.system('ln -s %s' % self.absFileName)
        return
    def filename(self):
        return self.fileName

class FileLinkWild(FileDescr):
    """must be given an absolute path"""
    def __init__(self, fileNameWild):
        # assume is absolute path
        self.absFileNameWild = fileNameWild
        return
    def __call__(self, **args):
        '''
        might be better to go to glob.glob() call here
        '''
        pLs = os.popen('ls %s' % self.absFileNameWild)
        files = pLs.readlines()
        pLs.close()
        for file in files:
            os.system('ln -s %s' % file)
        return

def resolveWild(fname):
    '''
    might be better to go to glob.glob() call here
    '''
    pipe = os.popen('ls %s' % (fname))
    tmp = pipe.read()
    pipe.close()
    return (string.split(tmp)[0])

def catList(lines,sep=''):
    val = ""
    for line in lines[:-1]:
        val = val + line + sep
    val = val + lines[-1]
    return val

def fileToForm(fname):
    f = open(fname,'r')
    fLines = f.readlines()
    f.close()
    form = catList(fLines)
    return form

class FileForm(FileDescr):
    def __init__(self, fileName, formFile=None, form=None, dictForDefs=None):
        self.fileName = fileName
        self.dictForDefs = dictForDefs
        if form:
            self.form = form
        elif formFile:
            self.form = fileToForm(formFile)
        else:
            self.form = fileToForm(fileName)
            # raise RuntimeError, 'need form or formFile'
        return
    def filename(self):
        return self.fileName
    def __call__(self, **args):
        f = open(self.fileName,'w')
        if self.dictForDefs is not None:
            print >> f, dictToDefs(self.dictForDefs)
        if args.has_key('dict'):
            valDict = args['dict']
            print >> f, self.form % valDict,
        else:
            valDict = None
            print >> f, self.form, # essentially just a copy
        f.close()
        return
    def __str__(self):
        dictTemp = {}
        dictTemp.update(self.__dict__)
        form = dictTemp.pop('form') # pop so that do not print what may be a long form
        return str(self.__class__) + "(" + str(self.__dict__) + ")"
    def __repr__(self):
        return str(self)

def readFloatDataA(fname):
    '''
    This is definitely faster than readFloatData when it is appropriate,
    which is pretty much any time the data can be interpreted as an array and
    numpy is available
    '''
    import numpy as num
    if hasattr(num, 'loadtxt'):
        data = num.loadtxt(fname, dtype=float, comments='#')
    else:
        data = readFloatData(fname)
    return data
#
def readFloatData(fname=None):
    data = []
    if fname:
      f = open(fname)
    else:
      f = sys.stdin
    line = f.readline()
    while line:

        # strip off end of line
        temp = re.split("\n",line)
        line = temp[0]

        # strip off trailing comment
        temp = re.split("#",line)
        line = temp[0]

        fields = re.split("[ \t]+",line)
        if len(fields[0]) == 0:
            # had some leading white space
            iFieldFirst = 1
        else:
            iFieldFirst = 0
        if len(fields[len(fields)-1]) == 0:
            # had some trailing white space
            iFieldLast = len(fields)-1
        else:
            iFieldLast = len(fields)

        if (iFieldLast >= iFieldFirst):
            dataLine = []
            for field in fields[iFieldFirst:iFieldLast] :
                dataLine.append(float(field))
            data.append(dataLine)

        line = f.readline()

    if fname:
      f.close()
    return data

def readDataFlat(fname):
    data = []
    f = open(fname)
    line = f.readline()
    while line:

        # strip off end of line
        temp = re.split('\n',line)
        line = temp[0]

        # strip off trailing comment
        temp = re.split('#',line)
        line = temp[0]

        fields = re.split('[ \t]+',line)
        if len(fields[0]) == 0:
            # had some leading white space
            iFieldFirst = 1
        else:
            iFieldFirst = 0
        if len(fields[len(fields)-1]) == 0:
            # had some trailing white space
            iFieldLast = len(fields)-1
        else:
            iFieldLast = len(fields)

        if (iFieldLast >= iFieldFirst):
            #dataLine = []
            for field in fields[iFieldFirst:iFieldLast] :
                #dataLine.append(eval(field))
                data.append(eval(field))
            #data.append(dataLine)

        line = f.readline()

    f.close()
    return data

def archiveDir(dirName):
    timeStamp = time.strftime('%Y_%m_%d__%H_%M_%S')
    newDir = dirName+"-"+timeStamp
    os.rename(dirName,newDir)
    # logger not necessarily set up yet
    #logger('moved directory %s to %s' % (dirName,newDir))
    print 'moved directory %s to %s' % (dirName,newDir)
    return

def archiveFile(fileName):
    timeStamp = time.strftime('%Y_%m_%d__%H_%M_%S')
    newFile = fileName+"-"+timeStamp
    os.rename(fileName,newFile)
    # logger not necessarily set up yet
    #logger('moved fileectory %s to %s' % (fileName,newFile))
    print 'moved file %s to %s' % (fileName,newFile)
    return

def getFromPipe(command, werr=False):
    import os
    if werr:
        import subprocess
        # pipei, pipeo = os.popen4(command)
        p = subprocess.Popen(command,
                             shell=True,
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             close_fds=True)
        pipei, pipeo = (p.stdin, p.stdout)
        val = pipeo.readline()[:-1] # [:-1] drops end-of-line character
        pipeo.close()
        pipei.close()
    else:
        pipe = os.popen(command)
        val = pipe.readline()[:-1] # [:-1] drops end-of-line character
        pipe.close()
    return val

def getSysType():
    import os
    if (os.environ.has_key('SYS_TYPE')):
        sysType = os.environ['SYS_TYPE']
    else:
        # raise RuntimeError, 'need SYS_TYPE environment variable'
        sysType = getFromPipe('babel-config --query-var=target')
    return sysType

def getNCorePerNode():
    import os, string
    sysType = getSysType()
    if sysType == 'aix_5_64_fed':
       return 8
    else: 
       pipe = os.popen('grep processor /proc/cpuinfo | wc --lines')
       tmp = pipe.readline()
       nCores = int((string.split(tmp))[0])
       pipe.close()
       return nCores

def getHostName():
    import re, os
    return re.split('[0-9]', os.environ['HOSTNAME'])[0]

def getBankNames():
    pMdiag = os.popen('mdiag -u %s' % (os.environ['USER']))
    line = pMdiag.readline()
    bankNames = []
    while line:
        splitLine = line.split('ALIST=')
        if len(splitLine) > 1:
            bankNames += splitLine[1].split(',')
        line = pMdiag.readline()
    pMdiag.close()
    return bankNames

def rmSafe(fileName):
    import os
    if os.access(fileName, os.F_OK) : os.remove(fileName)
    return

def rmDirF(dirName):
    import os
    for root, dirs, files in os.walk(dirName, topdown=False):
        for file in files: os.remove(os.path.join(root,file)) # os.remove(file)
        for dir  in dirs:  os.rmdir(os.path.join(root,dir))
    os.removedirs(dirName)
    return

def rmWorkDir(workdir):
    #os.rmdir(workdir)
    os.system('/bin/rm -rf %s' % (workdir))
    return

def listFiles(filename):
    'most useful if filename contains a wildcard'
    import os
    pLs = os.popen('ls '+filename)
    files = pLs.readlines()
    pLs.close()
    return [fTemp.split()[0] for fTemp in files]
    
def rmWild(filename):
    import os
    p = os.popen('ls %s*' % (filename))
    line = p.readline()
    while line:
      os.remove(line.split()[0])
      line = p.readline()
    return

def getScratchBaseDir():
    
    sysType = getSysType()
    resourceName = None
    #
    if sysType == 'chaos_3_x86_elan3':
        hostName = getHostName()
        if hostName == 'ingot':
            scratchBaseDir = '/nfs/petasim/' + os.environ['USER']
        elif hostName == 'alc':
            scratchBaseDir = '/p/lscratchb/' + os.environ['USER']
            resourceName = 'lscratchb'
        else:
            raise RuntimeError, 'unknown HOSTNAME environment variable  '
    elif sysType == 'chaos_3_x86_64_ib' or sysType == 'chaos_4_x86_64_ib':
        scratchBaseDir = '/p/lscratchb/' + os.environ['USER']
        resourceName = 'lscratchb'
    elif sysType == 'aix_5_64_fed':
        scratchBaseDir = '/p/gscratcha/' + os.environ['USER']
        resourceName = 'gscratcha'
    elif sysType.count('i386-apple-darwin'):
        scratchBaseDir = '/tmp'
    else:
        raise RuntimeError, 'getting scratch base directory: unknown SYS_TYPE environment variable  '
    return scratchBaseDir, resourceName

def argListToStr(argv):
    if len(argv) > 0:
        argsString = reduce(lambda a,b : str(a)+' '+str(b), argv)
    else:
        argsString = ''
    return argsString


def dictToDefs(d):
    s = ''
    for key, val in d.iteritems():
        if isinstance(val,float) or isinstance(val,int):
            s += 'def %s {%s}\n' % (key,val)
    return s


def readFLT(fname, structured=False):
    """
    Read a Fable-style .flt file

    JVB 2011/03/24
    """
    import numpy as num
    fid = open(fname, 'r')
    hdr = fid.readline().strip().split()
    nFields = len(hdr) - 1
    tmp = num.fromfile(fid, dtype=float, sep=" ")
    nRecords = len(tmp) / nFields
    tmp = tmp.reshape(nRecords, nFields)
    if structured:
        flt_dtype = [[hdr[i+1], 'f8'] for i in range(nFields)]
        intCols = [3,13,14,15,16,17,19,20,21,22,27,28,29]
        for i in intCols:
            flt_dtype[i][1] = 'i4'
        flt_dtype = [tuple(flt_dtype[i]) for i in range(nFields)]
        retval = num.empty(nRecords, dtype=flt_dtype)
        for i in range(nRecords):
            thisRow = list(tmp[i, :])
            thisRowInts = num.array(tmp[i, intCols], dtype='i4')
            jj = 0
            for j in intCols:
                thisRow[j] = thisRowInts[jj]
                jj+=1
            retval[i] = tuple(thisRow)
    else:
        retval = tmp
    fid.close()
    return retval

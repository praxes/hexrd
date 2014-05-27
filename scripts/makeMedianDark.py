import sys, os

from hexrd.xrd import detector
from hexrd.xrd import xrdutil

fileInfo = [('/home/bernier2/Documents/data/ruby_00025.ge3', 2), 
            ('/home/bernier2/Documents/data/ruby_00026.ge3', 2)]

nChunk      = 5
medianSize  = 3
medianRange = (-50,50)

def medianDarkGE(fileInfo, nChunk, medianSize, medianRange):
    reader    = detector.ReadGE(fileInfo, subtractDark=False, doFlip=False)    
    numFrames = reader.getNFrames()
    mdresults = xrdutil.darkFromStack(reader, nFrames=numFrames, nChunk=nChunk, 
                                      medianSize=medianSize, medianRange=medianRange, 
                                      checkIntensityResolution=False)
    darkFrame, deadPixels, std = mdresults
    return darkFrame, deadPixels, std

if __name__ == "__main__":
    """
    Usage:

    python makeMedianDark.py LT3MEW ge2 341 341 5 1 ./ dark
    """
    argv = sys.argv[1:]
    
    fileRoot   = argv[0]
    fileSuffix = argv[1]
    fileStart  = argv[2]
    fileStop   = argv[3]
    zpad       = argv[4]
    nEmtpy     = argv[5]
    outputDir  = argv[6]
    outputName = argv[7]
    
    numImages = (int(fileStop) - int(fileStart)) + 1
    fileInfo = []
    zpad_str = "_%"+"0%dd" % (int(zpad))
    
    for i in range(numImages):
        if fileSuffix.strip() == '':
            fname = fileRoot+zpad_str
        else:
            fname = fileRoot+zpad_str+"."+fileSuffix 
        fileInfo.append( ( fname % (int(fileStart) + i), int(nEmtpy) ) )
        pass
    
    print "processing fileInfo: "
    print fileInfo
    
    darkFrame, deadPixels, std = medianDarkGE(fileInfo, nChunk, medianSize, medianRange)
    fw = detector.FrameWriter(2048, 2048, 
                              filename=os.path.join(outputDir, outputName+'.'+fileSuffix), 
                              nbytesHeader=8192)
    fw.write(darkFrame)
    fw.close()

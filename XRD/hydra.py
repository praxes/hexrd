#! /usr/bin/env python
#
#  Copyright-Info-Goes-Here
#
"""Hydra detector tools.
"""
#
# ---------------------------------------------------CLASS:  Hydra
#
class Hydra(object):
    """Hydra image processing"""
    def __init__(self):
        """Constructor for Hydra."""
        #
        #  These arrays need four entries
        #
        self.readers = 4*[None]
	self.images  = 4*[None]
        #
        return
    #
    # ============================== API
    #
    #                     ========== Properties
    #
    #                     ========== Public Methods
    #
    def loadImages(self):
        """Load the four hydra images"""
        print 'loading images'
        reader1 = self.readers[0]
        aggMode = reader1.aggModeOp
        nrFrame = reader1.getNumberOfFrames() # number of reader frames
        if aggMode:
            rdFrames = nrFrame
        else:
            rdFrames = 1
            pass
        
        for i in range(4):
            ri = self.readers[i].makeReader()
            self.images[i] = ri.read(nframes= rdFrames, 
                                     sumImg = aggMode)
            print 'done %d' % i
            pass

        return
    #
    pass  # end class
#
# -----------------------------------------------END CLASS:  Hydra


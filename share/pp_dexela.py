import time

from hexrd import imageseries

PIS = imageseries.process.ProcessedImageSeries

class PP_Dexela(object):
    """PP_Dexela"""
    PROCFMT = 'frame-cache'
    RAWFMT = 'hdf5'
    RAWPATH = '/imageseries'
    DARKPCTILE = 50

    def __init__(self, fname, omw, frame_start=0):
        """Constructor for PP_Dexela"""
        #
	self.fname = fname
        self.omwedges = omw
        self.frame_start = frame_start
        self.use_frame_list = (self.frame_start > 0)
        self.raw = imageseries.open(self.fname, self.RAWFMT, path=self.RAWPATH)
        self._dark = None

        return

    @property
    def oplist(self):
        return [('dark', self.dark ), ('flip', 't'), ('flip', 'hv') ]

    @property
    def framelist(self):
        return range(self.frame_start, self.nframes+self.frame_start)
    #
    # ============================== API
    #
    @property
    def nframes(self):
        return self.omwedges.nframes

    def omegas(self):
        return self.omwedges.omegas

    def save_omegas(self, fname):
        self.omwedges.save_omegas(fname)

    def processed(self):
        if self.use_frame_list:
            kw = dict(frame_list=self.framelist)

        return PIS(self.raw, self.oplist, **kw)

    @property
    def dark(self, nframes=50):
        """build and return dark image"""
        if self._dark is None:
            print "building dark images using %s frames (may take a while) ... " % nframes
            start = time.clock()
            self._dark = imageseries.stats.percentile(
                    self.raw, self.DARKPCTILE, nframes=nframes
            )
            elapsed = (time.clock() - start)
            print "done building background (dakr) image: elapsed time is %f seconds" \
              % elapsed

        return self._dark

    def save_processed(self, name, threshold):
        fcname = '%s-fcache.yml' % name
        cache = '%s-cachefile.npz' % name
        omname = '%s-omegas.npy' % name
        imageseries.write(self.processed(), fcname, self.PROCFMT,
                  threshold=threshold,
                  cache_file=cache)
        self.save_omegas(omname)

    pass  # end class

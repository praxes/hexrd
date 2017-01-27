import time
import os

from hexrd import imageseries

PIS = imageseries.process.ProcessedImageSeries

class PP_Dexela(object):
    """PP_Dexela"""
    PROCFMT = 'frame-cache'
    RAWFMT = 'hdf5'
    RAWPATH = '/imageseries'
    DARKPCTILE = 50

    def __init__(self, fname, omw, flips, frame_start=0):
        """Constructor for PP_Dexela"""
        #
	self.fname = fname
        self.omwedges = omw
        self.flips = flips
        self.frame_start = frame_start
        self.use_frame_list = (self.frame_start > 0)
        self.raw = imageseries.open(self.fname, self.RAWFMT, path=self.RAWPATH)
        self._dark = None

        print 'On Init: ', self.nframes, self.fname, self.omwedges.nframes,\
          len(self.raw)
        return

    @property
    def oplist(self):
        return [('dark', self.dark)] + self.flips

    @property
    def framelist(self):
        return range(self.frame_start, self.nframes)
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
            usenframes = min(nframes, self.nframes)
            print "building dark images using %s frames (may take a while)"\
            " ... " % usenframes
            start = time.clock()
            self._dark = imageseries.stats.percentile(
                    self.raw, self.DARKPCTILE, nframes=usenframes
            )
            elapsed = (time.clock() - start)
            print "done building background (dakr) image: elapsed time is %f seconds" \
              % elapsed

        return self._dark

    def save_processed(self, name, threshold):
        dname = '%s-fcache-dir' % name
        tcname = '%s-fcache-tmp.yml' % name
        fcname = '%s-fcache.yml' % name
        cache = '%s-cachefile.npz' % name
        omname = '%s-omegas.npy' % name

        pname = lambda s: os.path.join(dname, s) # prepend fc directory

        os.mkdir(dname)

        # Steps:
        # * write frame cache with no omegas to temporary file
        # * write omegas to file
        # * modify temporary file to include omegas
        imageseries.write(self.processed(), pname(tcname), self.PROCFMT,
                  threshold=threshold,
                  cache_file=cache)
        self.save_omegas(pname(omname))
        # modify yaml
        with open(pname(tcname), 'r') as f:
            s = f.read()
        m0 = 'meta: {}'
        m1 = 'meta:\n  omega: ! load-numpy-array %s' % omname
        with open(pname(fcname), 'w') as f:
            f.write(s.replace(m0, m1))
        os.remove(pname(tcname))

    pass  # end class

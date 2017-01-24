"""Handle omega (specimen rotation) metadata

* OmegaWedges class specifies omega metadata in wedges
"""
import numpy as np

from .baseclass import ImageSeries

OMEGA_KEY = 'omega'

class OmegaImageSeries(ImageSeries):
    """ImageSeries with omega metadata"""
    DFLT_TOL = 1.0e-6

    def __init__(self, ims):
        """This class is initialized with an existing imageseries"""
        # check for omega metadata
        if OMEGA_KEY in ims.metadata:
            self._omega = ims.metadata[OMEGA_KEY]
            if len(ims) != self._omega.shape[0]:
                msg = 'omega array mismatch: array has %s frames, expecting %s'
                msg = msg % (self._omega.shape[0], len(ims))
                raise OmegaSeriesError(msg)
        else:
            raise OmegaSeriesError('Imageseries has no omega metadata')

        super(OmegaImageSeries, self).__init__(ims)
        self._make_wedges()

    def _make_wedges(self, tol=DFLT_TOL):
        nf = len(self)
        om = self.omega

        # find the frames where the wedges break
        starts = [0]
        delta = om[0, 1] - om[0, 0]
        omlast = om[0, 1]
        for f in range(1, nf):
            if delta <= 0:
                raise OmegaSeriesError('omega array must be increasing')
            # check whether delta changes or ranges not contiguous
            d = om[f,1] - om[f,0]
            if (np.abs(d - delta) > tol) or (np.abs(om[f,0] - omlast) > tol):
                starts.append(f)
                delta = d
            omlast = om[f, 1]
        starts.append(nf)

        self._omegawedges = OmegaWedges(nf)
        for s in range(len(starts) - 1):
            ostart = om[starts[s], 0]
            ostop = om[starts[s + 1] - 1, 1]
            steps = starts[s+1] - starts[s]
            self._omegawedges.addwedge(ostart, ostop, steps)

    @property
    def omega(self):
        """return omega range array (nframes, 2)"""
        return self._omega

    @property
    def omegawedges(self):
        return self._omegawedges

    @property
    def nwedges(self):
        return self.omegawedges.nwedges

    def wedge(self, i):
        """return i'th wedge as a dictionary"""
        d = self.omegawedges.wedges[i]
        delta = (d['ostop'] - d['ostart'])/d['nsteps']
        d.update(delta=delta)
        return d


class OmegaWedges(object):
    """piecewise linear omega ranges"""
    def __init__(self, nframes):
        """Constructor for OmegaWedge"""
	self.nframes = nframes
        self._wedges = []
    #
    # ============================== API
    #
    @property
    def omegas(self):
        """n x 2 array of omega values, one per frame"""
        if self.nframes != self.wframes:
            msg = "number of frames (%s) does not match "\
                  "number of wedge frames (%s)" %(self.nframes, self.wframes)
            raise OmegaSeriesError(msg)

        oa = np.zeros((self.nframes, 2))
        wstart = 0
        for w in self.wedges:
            ns = w['nsteps']
            wr = range(wstart, wstart + ns)
            wa0 = np.linspace(w['ostart'], w['ostop'], ns + 1)
            oa[wr, 0] = wa0[:-1]
            oa[wr, 1] = wa0[1:]
            wstart += ns

        return oa

    @property
    def nwedges(self):
        """number of wedges"""
        return len(self._wedges)

    @property
    def wedges(self):
        """list of wedges (dictionaries)"""
        return self._wedges

    def addwedge(self, ostart, ostop, nsteps, loc=None):
        """add wedge to list"""
        d = dict(ostart=ostart, ostop=ostop, nsteps=nsteps)
        if loc is None:
            loc = self.nwedges

        self.wedges.insert(loc, d)

    def delwedge(self, i):
        """delete wedge number i"""
        self.wedges.pop(i)

    @property
    def wframes(self):
        """number of frames in wedges"""
        wf = [w['nsteps'] for w in self.wedges]
        return np.int(np.sum(wf))

    def save_omegas(self, fname):
        """save omegas to text file"""
        np.save(fname, self.omegas)

    pass  # end class


class OmegaSeriesError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

"""Handle omega (specimen rotation) metadata

* OmegaWedges class specifies omega metadata in wedges
"""
import numpy as np

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
            raise OmegaWedgesError(msg)

        oa = np.zeros((self.nframes, 2))
        wstart = 0
        for w in self.wedges:
            ns = w['nsteps']
            wr = range(wstart, wstart + ns)
            wa0 = np.linspace(w['ostart'], w['ostop'], ns + 1)
            oa[wr, 0] = wa0[:-1]
            oa[wr, 1] = wa0[1:]

        return oa

    @property
    def nwedges(self):
        """number of wedges"""
        return len(self._wedges)

    @property
    def wedges(self):
        """list of wedges (dictionary)"""
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


class OmegaWedgesError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

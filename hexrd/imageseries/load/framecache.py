"""Adapter class for frame caches
"""
import os

import numpy as np
from scipy.sparse import csr_matrix
import yaml

from . import ImageSeriesAdapter
from ..imageseriesiter import ImageSeriesIterator
from .metadata import yamlmeta

class FrameCacheImageSeriesAdapter(ImageSeriesAdapter):
    """collection of images in HDF5 format"""

    format = 'frame-cache'

    def __init__(self, fname, **kwargs):
        """Constructor for frame cache image series

        *fname* - filename of the yml file
        *kwargs* - keyword arguments (none required)
        """
        self._fname = fname
        self._load_yml()
        self._load_cache()

    def _load_yml(self):
        with open(self._fname, "r") as f:
            d = yaml.load(f)
        datad = d['data']
        self._cache = datad['file']
        self._nframes = datad['nframes']
        self._shape = tuple(datad['shape'])
        self._dtype = np.dtype(datad['dtype'])
        self._meta = yamlmeta(d['meta'], path=self._cache)

    def _load_cache(self):
        """load into list of csr sparse matrices"""
        bpath = os.path.dirname(self._fname)
        if os.path.isabs(self._cache):
            cachepath = self._cache
        else:
            cachepath = os.path.join(bpath, self._cache)
        arrs = np.load(cachepath)

        self._framelist = []
        for i in range(self._nframes):
            row = arrs["%d_row" % i]
            col = arrs["%d_col" % i]
            data = arrs["%d_data" % i]
            frame = csr_matrix((data, (row, col)),
                               shape=self._shape, dtype=self._dtype)
            self._framelist.append(frame)

    @property
    def metadata(self):
        """(read-only) Image sequence metadata
        """
        return self._meta

    def load_metadata(self, indict):
        """(read-only) Image sequence metadata

        Currently returns none
        """
        #### Currently not used: saved temporarily for np.array trigger
        metad = {}
        for k, v in indict.items():
            if v == '++np.array':
                newk = k + '-array'
                metad[k] = np.array(indict.pop(newk))
                metad.pop(newk, None)
            else:
                metad[k] = v
        return metad

    @property
    def dtype(self):
        return self._dtype

    @property
    def shape(self):
        return self._shape

    def __getitem__(self, key):
        return self._framelist[key].toarray()

    def __iter__(self):
        return ImageSeriesIterator(self)

    #@memoize
    def __len__(self):
        return self._nframes

    pass  # end class

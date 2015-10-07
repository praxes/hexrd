"""Adapter class for frame caches
"""
from . import ImageSeriesAdapter
from ..imageseriesiter import ImageSeriesIterator

import numpy as np
from scipy.sparse import csr_matrix
import yaml

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
        metad = d['meta']
        self._cache = datad['file']
        self._nframes = datad['nframes']
        self._shape = tuple(datad['shape'])
        self._dtype = np.dtype(datad['dtype'])
        self._meta = metad

    def _load_cache(self):
        """load into list of csr sparse matrices"""
        arrs = np.load(self._cache)

        self._framelist = []
        for i in range(self._nframes):
            row = arrs["%d_row" % i]
            col = arrs["%d_col" % i]
            data = arrs["%d_data" % i]
            frame = csr_matrix((data, (row, col)), shape=self._shape)
            self._framelist.append(frame)

    @property
    def metadata(self):
        """(read-only) Image sequence metadata

        Currently returns none
        """
        return self._meta

    @property
    def dtype(self):
        return self._dtype

    @property
    def shape(self):
        return self._shape

    def __getitem__(self, key):
        return self._framelist[key]

    def __iter__(self):
        return ImageSeriesIterator(self)

    #@memoize
    def __len__(self):
        return self._nframes

    pass  # end class

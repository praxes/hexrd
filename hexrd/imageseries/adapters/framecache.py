"""Adapter class for frame caches
"""
from . import ImageSeriesAdapter
from ..imageseriesiter import ImageSeriesIterator

import numpy as np
from scipy.sparse import csr_matrix

class FrameCacheImageSeriesAdapter(ImageSeriesAdapter):
    """collection of images in HDF5 format"""

    format = 'frame-cache'

    def __init__(self, fname, **kwargs):
        """Constructor for frame cache image series

        *fname* - filename of the yml file 
        *kwargs* - keyword arguments (none required 
        """
        self._fname = fname
        self._load()

    def _load(self):
        """load into list of csr sparse matrices"""
        arrs = np.load(self._fname)
        self._shape = tuple(arrs['shape'].tolist())
        self._dtype = None

        arrsh = tuple(arrs['shape'])
        nk = len(arrs.files) - 1
        self._nframes = nk/3
        self._framelist = []
        for i in range(self._nframes):
            row = arrs["%d_row" % i]
            col = arrs["%d_col" % i]
            data = arrs["%d_data" % i]
            frame = csr_matrix((data, (row, col)), shape=arrsh)
            self._framelist.append(frame)
            if self._dtype is None:
                self._dtype = data.dtype
    
    def metadata(self):
        """(read-only) Image sequence metadata

        Currently returns none
        """
        return None
    
    def shape(self):
        return self._dtype


    def __getitem__(self, key):
        with self._dset as dset:
            return dset.__getitem__(key)

    def __iter__(self):
        return ImageSeriesIterator(self)

    #@memoize
    def __len__(self):
        return self._nframes

    pass  # end class

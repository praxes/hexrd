"""Adapter class for numpy array (3D)
"""
from . import ImageSeriesAdapter
from ..imageseriesiter import ImageSeriesIterator

class ArrayImageSeriesAdapter(ImageSeriesAdapter):
    """collection of images in numpy array"""

    format = 'array'

    def __init__(self, fname, **kwargs):
        """Constructor for frame cache image series

        *fname* - should be None 
        *kwargs* - keyword arguments
                 . 'data' = a 3D array (double/float)
                 . 'metadata' = a dictionary
        """
        self._data = kwargs['data']
        self._meta = kwargs.pop('metadata', dict())
        self._shape = self._data.shape
        self._nframes = self._shape[0]
        self._nxny = self._shape[1:3]

    @property
    def metadata(self):
        """(read-only) Image sequence metadata

        Currently returns none
        """
        return self._meta

    @property
    def shape(self):
        return self._nxny

    @property
    def dtype(self):
        return self._data.dtype
    
    def __getitem__(self, key):
        return self._data[key]

    def __iter__(self):
        return ImageSeriesIterator(self)

    #@memoize
    def __len__(self):
        return self._nframes

    pass  # end class

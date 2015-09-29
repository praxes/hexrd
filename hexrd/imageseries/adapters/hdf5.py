"""HDF5 adapter class
"""
import h5py
from . import ImageSeriesAdapter
from ..imageseriesiter import ImageSeriesIterator

class HDF5ImageSeriesAdapter(ImageSeriesAdapter):
    """collection of images in HDF5 format"""

    format = 'hdf5'
    #The code below failed with: "Error when calling the metaclass bases"
    #                  "'property' object is not callable"
    #@property
    #def format(self):
    #    return 'hdf5'

    @property
    def _dset(self):
        # return a context manager to ensure proper file handling
        # always use like: "with self._dset as dset:"
        return H5ContextManager(self.__h5name, self.__path)

    @property
    #@memoize
    def metadata(self):
        """(read-only) Image sequence metadata

        Currently returns any dimension scales in a dictionary
        """
        mdict = {}
        with self._dset as dset:
            for k in dset.dims[0].keys():
                mdict[k] = dset.dims[0][k][...]

        return mdict

    @property
    def dtype(self):
        with self._dset as dset:
            return dset.dtype
        
    @property
    #@memoize so you only need to do this once
    def shape(self):
        with self._dset as dset:
            return dset.shape[1:]

    def __init__(self, fname, **kwargs):
        """Constructor for H5FrameSeries

        *fname* - filename of the HDF5 file
        *kwargs* - keyword arguments, choices are:
           path - (required) path of dataset in HDF5 file
        """
        self.__h5name = fname
        self.__path = kwargs['path']

    def __getitem__(self, key):
        with self._dset as dset:
            return dset.__getitem__(key)

    def __iter__(self):
        return ImageSeriesIterator(self)

    #@memoize
    def __len__(self):
        with self._dset as dset:
            return len(dset)

    pass  # end class


class H5ContextManager:

    def __init__(self, fname, path):
        self._fname = fname
        self._path = path
        self._f = None

    def __enter__(self):
        self._f = h5py.File(self._fname, 'r')
        return self._f[self._path]

    def __exit__(self, *args):
        self._f.close()


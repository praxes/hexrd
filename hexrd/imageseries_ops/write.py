"""Write imageseries to various formats"""

import abc

import numpy
import h5py

def write(ims, fname, fmt, opts):
    """write imageseries to file with options

    *ims* - an imageseries
    *fname* - name of file
    *fmt* - a format string
    *opts* - namespace of options
    """
    print('writing format ', fmt, ' to file: ', fname)
    wcls = _Registry.getwriter(fmt)
    print(wcls)
    w = wcls(ims, fname, opts)
    w.write()

# Registry

class _RegisterWriter(abc.ABCMeta):

    def __init__(cls, name, bases, attrs):
        abc.ABCMeta.__init__(cls, name, bases, attrs)
        _Registry.register(cls)

class _Registry(object):
    """Registry for imageseries writers"""
    writer_registry = dict()

    @classmethod
    def register(cls, wcls):
        """Register writer class"""
        if wcls.__name__ is not 'Writer':
            cls.writer_registry[wcls.fmt] = wcls

    @classmethod
    def getwriter(cls, name):
        """return instance associated with name"""
        return cls.writer_registry[name]
    #
    pass  # end class

class Writer(object):
    """Base class for writers"""
    __metaclass__ = _RegisterWriter
    fmt = None
    def __init__(self, ims, fname, opts):
        self._ims = ims
        self._shape = ims.shape
        self._dtype = ims.dtype
        self._nframes = len(ims)
        self._fname = fname
        self._opts = opts

    pass # end class
        
class WriteH5(Writer):
    fmt = 'hdf5'

    def __init__(self, ims, fname, opts):
        Writer.__init__(self, ims, fname, opts)
        self._path = self._opts['path']
    
    def _open_dset(self):
        """open HDF5 file and dataset"""
        f = h5py.File(self._fname, "a")
        s0, s1 = self._shape
        
        return f.create_dataset(self._path, (self._nframes, s0, s1), self._dtype,
                                compression="gzip")
    #
    # ======================================== API
    #
    def write(self):
        """Write imageseries to HDF5 file"""
        print('writing ', self.fmt)
        ds = self._open_dset()
        for i in range(self._nframes):
            ds[i, :, :] = self._ims[i]

        # next: add metadata
        
    pass # end class

class WriteFrameCache(Writer):

    fmt = 'frame-cache'
    def __init__(self, ims, fname, opts):
        Writer.__init__(self, ims, fname, opts)
        self._thresh = self._opts['threshold']

    def write(self):
        """writes frame cache for imageseries

        presumes sparse forms are small enough to contain all frames
        """
        print('writing ', self.fmt)
        arrd = dict()
        for i in range(self._nframes):
            frame = self._ims[i]
            mask = frame > self._thresh
            row, col = mask.nonzero()
            arrd['%d_data' % i] = frame[mask]
            arrd['%d_row' % i] = row
            arrd['%d_col' % i] = col
        
        numpy.savez_compressed(self._fname, **arrd)

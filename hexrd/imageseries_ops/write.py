"""Write imageseries to various formats"""

import abc

import h5py

def write(ims, fname, fmt, opts):
    """write imageseries to file with options

    *ims* - an imageseries
    *fname* - name of file
    *fmt* - a format string
    *opts* - namespace of options
    """
    print('writing format ', fmt, ' to file: ', fname)
    w = _Registry.getwriter(fmt)
    w.write(ims, fname, opts)
    

# Registry

class _RegisterWriter(abc.ABCMeta):

    def __init__(cls, name, bases, attrs):
        abc.ABCMeta.__init__(cls, name, bases, attrs)
        _Registry.register(cls)

class _Registry(object):
    """Registry for symmetry type instances"""
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
        
class WriteH5(Writer):
    fmt = 'hdf5'

    @classmethod
    def write(cls, ims, fname, opts):
        """Write imageseries to HDF5 file"""
        print('writing hdf5')
        pass
        
    pass # end class

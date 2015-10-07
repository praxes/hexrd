import abc
import pkgutil

from ..imageseriesabc import ImageSeriesABC
from .registry import Registry

# Metaclass for adapter registry

class _RegisterAdapterClass(abc.ABCMeta):

    def __init__(cls, name, bases, attrs):
        abc.ABCMeta.__init__(cls, name, bases, attrs)
        Registry.register(cls)

class ImageSeriesAdapter(ImageSeriesABC):

    __metaclass__ = _RegisterAdapterClass

    format = None

# import all adapter modules

for loader, name, ispkg in pkgutil.iter_modules(__path__):
    if name is not 'registry':
        __import__(name, globals=globals())
    #
    # couldn't get the following line to work due to relative import issue:
    #     loader.find_module(name).load_module(name)

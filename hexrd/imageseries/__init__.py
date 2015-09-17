"""Handles series of images

This file contains the generic ImageSeries class
and a function for loading. Adapters for particular
data formats are managed in a subpackage.
"""
from imageseriesabc import ImageSeriesABC
import adapters

class ImageSeries(ImageSeriesABC):
    """collection of images

    Basic sequence class with additional properties for image shape and
    metadata (possibly None).
    """

    def __init__(self, adapter):
        """Build FrameSeries from adapter instance

        *adapter* - object instance based on abstract Sequence class with
        properties for image shape and, optionally, metadata.
        """
        self.__adapter = adapter

        return

    def __getitem__(self, key):
        return self.__adapter[key]

    def __len__(self):
        return len(self.__adapter)

    def __getattr__(self, attrname):
        return getattr(self.__adapter, attrname)

    pass  # end class

def open(filename, format=None, **kwargs):
    # find the appropriate adapter based on format specified
    reg = adapters.Registry.adapter_registry
    adapter = reg[format](filename, **kwargs)
    return ImageSeries(adapter)

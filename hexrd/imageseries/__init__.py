"""Handles series of images

This file contains the generic ImageSeries class
and a function for loading. Adapters for particular
data formats are managed in the "load" subpackage.
"""
from .baseclass import ImageSeries
from . import load
from . import process
from . import save
from . import stats

def open(filename, format=None, **kwargs):
    # find the appropriate adapter based on format specified
    reg = load.Registry.adapter_registry
    adapter = reg[format](filename, **kwargs)
    return ImageSeries(adapter)

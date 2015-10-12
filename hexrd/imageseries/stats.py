"""Stats for imageseries"""
import numpy as np

def max(ims, nframes=0):
    nf = _nframes(ims, nframes)
    imgmax = ims[0]
    for i in range(1, nf):
        imgmax = np.maximum(imgmax, ims[i])
    return imgmax

def median(ims, nframes=0):
    """return image with median values over all frames"""
    # could be done by rectangle by rectangle if full series
    # too  big for memory
    nf = _nframes(ims, nframes)
    return np.median(_toarray(ims, nf), axis=0)

#
# ==================== Utilities
#
def _nframes(ims, nframes):
    """number of frames to use: len(ims) or specified number"""
    mynf = len(ims)
    return np.min((mynf, nframes)) if nframes > 0 else mynf

def _toarray(ims, nframes):
    ashp = (nframes,) + ims.shape
    a = np.zeros(ashp, dtype=ims.dtype)
    for i in range(nframes):
        a[i] = ims[i]

    return a

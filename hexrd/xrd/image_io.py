"""Image reading (mostly) and writing

Classes
-------
Framer2DRC: base class for reader/writers
ReadGeneric:
ReadGE:

ThreadReadFrame: class for using threads to read frames

Functions
---------
newGenericReader - returns a reader instance

"""
import copy
import os
import time
import warnings

import numpy as num

import imageseries

warnings.filterwarnings('always', '', DeprecationWarning)

class ReaderDeprecationWarning(DeprecationWarning):
    """Warnings on use of old reader features"""
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class OmegaImageSeries(object):
    """Facade for frame_series class, replacing other readers, primarily ReadGE"""
    OMEGA_TAG = 'omega'

    def __init__(self, fname, fmt='hdf5', **kwargs):
        """Initialize frame readerOmegaFrameReader

        *fileinfo* is a string
        *fmt* is the format to be passed to imageseries.open()
        *kwargs* is the option list to be passed to imageseries.open()

        NOTES:
        * The shape returned from imageseries is cast to int from numpy.uint64
          to allow for addition of indices with regular ints
        """
        self._imseries = imageseries.open(fname, fmt, **kwargs)
        self._shape = self._imseries.shape
        self._meta = self._imseries.metadata

        if self.OMEGA_TAG not in self._meta:
            raise RuntimeError('No omega data found in data file')
        return

    def __getitem__(self, k):
        return self._imseries[k]

    @property
    def nframes(self):
        """(get-only) number of frames"""
        return self._shape[0]

    @property
    def nrows(self):
        """(get-only) number of rows"""
        return self._shape[1]

    @property
    def ncols(self):
        """(get-only) number of columns"""
        return self._shape[2]

    @property
    def omega(self):
        """ (get-only) array of omega min/max per frame"""
        return self._meta[self.OMEGA_TAG]

    pass

class Framer2DRC(object):
    """Base class for readers.
    """
    def __init__(self, ncols, nrows,
                 dtypeDefault='int16', dtypeRead='uint16', dtypeFloat='float64'):
        self.__nrows = nrows
        self.__ncols = ncols
        self.__frame_dtype_dflt  = dtypeDefault
        self.__frame_dtype_read  = dtypeRead
        self.__frame_dtype_float = dtypeFloat

        self.__nbytes_frame  = num.nbytes[dtypeRead]*nrows*ncols

        return

    def get_nrows(self):
        return self.__nrows
    nrows = property(get_nrows, None, None)

    def get_ncols(self):
        return self.__ncols
    ncols = property(get_ncols, None, None)

    def get_nbytesFrame(self):
        return self.__nbytes_frame
    nbytesFrame = property(get_nbytesFrame, None, None)

    def get_dtypeDefault(self):
        return self.__frame_dtype_dflt
    dtypeDefault = property(get_dtypeDefault, None, None)

    def get_dtypeRead(self):
        return self.__frame_dtype_read
    dtypeRead = property(get_dtypeRead, None, None)

    def get_dtypeFloat(self):
        return self.__frame_dtype_float
    dtypeFloat = property(get_dtypeFloat, None, None)

    def getEmptyMask(self):
        """convenience method for getting an emtpy mask"""
        # this used to be a class method
        return num.zeros([self.nrows, self.ncols], dtype=bool)

class OmegaFramer(object):
    """Omega information associated with frame numbers"""
    def __init__(self, omegas):
        """Initialize omega ranges

        *omegas* is nframes x 2

        Could check for monotonicity.
        """
        self._omegas = omegas
        self._omin = omegas.min()
        self._omax = omegas.max()
        self._omean = omegas.mean(axis=1)
        self._odels = omegas[:, 1] - omegas[:, 0]
        self._delta = self._odels[0]
        self._orange = num.hstack((omegas[:, 0], omegas[-1, 1]))

        return

    # property:

    def _omin(self):
        return self._omin

    def _omax(self):
        return self._omax

    def getDeltaOmega(self, nframes=1):
        return self._omax - self._omin

    def getOmegaMinMax(self):
        return self._omin, self._omax

    def frameToOmega(self, frame):
        """can frame be nonintegral? round to int ... """
        return self._omean[frame]

    def omegaToFrame(self, omega):
        return  num.searchsorted(self._orange) - 1


    def omegaToFrameRange(self, omega):
        # note: old code assumed single delta omega
        return omeToFrameRange(omega, self._omean, self._delta)


class ReadGeneric(Framer2DRC, OmegaFramer):
    """Generic reader with omega information
"""
    def __init__(self, filename, ncols, nrows, *args, **kwargs):

        Framer2DRC.__init__(self, ncols, nrows, **kwargs)
        return

    def read(self, nskip=0, nframes=1, sumImg=False):
        """
        sumImg can be set to True or to something like numpy.maximum
        """
        raise RuntimeError("Generic reader not available for reading")

    def getNFrames(self):
        return 0


    def getWriter(self, filename):
        return None

class ReadGE(Framer2DRC,OmegaFramer):
    """General reader for omega scans

    Originally, this was for reading GE format images, but this is now
    a general reader accessing the OmegaFrameReader facade class. The main
    functionality to read a sequence of images with associated omega ranges.

    ORIGINAL DOCS
    =============

    *) In multiframe images where background subtraction is requested but no
       dark is specified, attempts to use the
       empty frame(s).  An error is returned if there are not any specified.
       If there are multiple empty frames, the average is used.

    """
    def __init__(self, file_info, *args, **kwargs):
        """Initialize the reader

        *file_info* is now just the filename
        *kwargs* is a dictionary
                 keys include: "path" path in hdf5 file

        Of original kwargs, only using "mask"
        """
        self._fname = file_info
        self._kwargs = kwargs
        try:
            self._omis = OmegaImageSeries(file_info, **kwargs)
            Framer2DRC.__init__(self, self._omis.nrows, self._omis.ncols)
            OmegaFramer.__init__(self, self._omis.omega)
        except:
            self._omis = None
        self.mask = None


        # counter for last global frame that was read
        self.iFrame = -1

        return


    def __call__(self, *args, **kwargs):
        return self.read(*args, **kwargs)

    @classmethod
    def makeNew(cls):
        """return another copy of this reader"""
        return cls(self._fname, **self._kwargs)

    def getWriter(self, filename):
        return None

    def getNFrames(self):
        """number of total frames with real data, not number remaining"""
        return self._omis.nframes

    def getFrameOmega(self, iFrame=None):
        """if iFrame is none, use internal counter"""
        if iFrame is None:
            iFrame = self.iFrame
        if hasattr(iFrame, '__len__'):
            # in case last read was multiframe
            oms = [self.frameToOmega(frm) for frm in iFrame]
            retval = num.mean(num.asarray(oms))
        else:
            retval = self.frameToOmega(iFrame)
        return retval


    def readBBox(self, bbox, raw=True, doFlip=None):
        """
        with raw=True, read more or less raw data, with bbox = [(iLo,iHi),(jLo,jHi),(fLo,fHi)]

        """
        # implement in OmegaFrameReader
        nskip = bbox[2][0]
        bBox = num.array(bbox)
        sl_i = slice(*bBox[0])
        sl_j = slice(*bBox[1])
        'plenty of performance optimization might be possible here'
        if raw:
            retval = num.empty( tuple(bBox[:,1] - bBox[:,0]), dtype=self.__frame_dtype_read )
        else:
            retval = num.empty( tuple(bBox[:,1] - bBox[:,0]), dtype=self.__frame_dtype_dflt )
        for iFrame in range(retval.shape[2]):
            thisframe = reader.read(nskip=nskip)
            nskip = 0
            retval[:,:,iFrame] = copy.deepcopy(thisframe[sl_i, sl_j])
        return retval

    def getDark(self):
        return 0

    def indicesToMask(self, indices):
      """Create mask from list of indices

      Indices can be a list of indices, as from makeIndicesTThRanges
      """
      mask = self.getEmptyMask()
      if hasattr(indices,'__len__'):
        for indThese in indices:
          mask[indThese] = True
      else:
        mask[indices] = True
      return mask

    def read(self, nskip=0, nframes=1, sumImg=False):
        """Read one or more frames, possibly operating on them

        This returns a single frame is nframes is 1, multiple
        frames if nframes > 1 with sumImg off, or a single frame
        resulting from some operation on the multiple frames if
        sumImg is true or a function.

        *sumImg* can be set to True or to a function of two frames like numpy.maximum
        *nskip* applies only to the first frame
        """
        self.iFrame = num.atleast_1d(self.iFrame)[-1] + nskip

        multiframe = nframes > 1
        sumimg_callable = hasattr(sumImg, '__call__')

        if not multiframe:
            self.iFrame += 1
            img = self._omis[self.iFrame]
            if self.mask is not None:
                img[self.mask] = 0
            return img

        # multiframe case
        self.iFrame = self.iFrame + 1 + range(nframes)

        if not sumImg:
            # return multiple frames
            imgs = self._omis[self.iFrame]
            for i in range(nframes):
                if self.mask is not None:
                    imgs[i, self.mask] = 0
            return imgs

        # Now, operate on frames consecutively
        op = sumImg if sumimg_callable else num.add

        ifrm = self.iFrame + 1

        img = self._omis[ifrm]
        for i in range(1, nframes):
            ifrm += 1
            img = op(img, self._omis[ifrm])
        if not sumimg_callable:
            img = img * (1.0/nframes)

        if self.mask is not None:
            img[self.mask] = 0

        return img

    def close(self):
        return

    @classmethod
    def display(cls,
                thisframe,
                roi = None,
                pw  = None,
                **kwargs
                ):
        warnings.warn('display method on readers no longer implemented',
                      ReaderDeprecationWarning)

#
# Module functions
#
def omeToFrameRange(omega, omegas, omegaDelta):
    """
    check omega range for the frames in
    stead of omega center;
    result can be a pair of frames if the specified omega is
    exactly on the border
    """
    retval = num.where(num.abs(omegas - omega) <= omegaDelta*0.5)[0]
    return retval

def newGenericReader(ncols, nrows, *args, **kwargs):
    """ Currently just returns a Framer2DRC
    """

    # retval = Framer2DRC(ncols, nrows, **kwargs)
    filename = kwargs.pop('filename', None)
    retval = ReadGeneric(filename, ncols, nrows, *args, **kwargs)

    return retval


"""Class for processing frames or frame groups"""
import numpy as np

from .baseclass import ImageSeries

class ProcessedImageSeries(ImageSeries):
    """Images series with mapping applied to frames"""
    FLIP = 'flip'
    DARK = 'dark'
    RECT = 'rectangle'

    def __init__(self, imser, **kwargs):
        """Instantiate imsageseries based on existing one with mapping options

        *imser* - an existing imageseries
        *kwargs* - dictionary for processing options
        """
        self._imser = imser
        self._opts = kwargs

    def __getitem__(self, key):
        return self._process_frame(key)

    def __len__(self):
        return len(self._imser)

    def _process_frame(self, key):
        # apply flip at end
        img = self._imser[key]
        img = self._subtract_dark(img)
        img = self._rectangle(img)
        img = self._flip(img)
        return img

    def _subtract_dark(self, img):
        # need to check for values below zero
        if self.DARK not in self._opts:
            return img

        dark = self._opts[self.DARK]
        return np.where(img > dark, img-dark, 0)

    def _rectangle(self, img):
        # restrict to rectangle
        if self.RECT in self._opts:
            r = self._opts[self.RECT]
            return img[r[0,0]:r[0,1], r[1,0]:r[1,1]]
        else:
            return img

    def _flip(self, img):
        if self.FLIP in self._opts:
            flip = self._opts['flip']
        else:
            return img

        if flip in ('y','v'): # about y-axis (vertical)
            pimg = img[:, ::-1]
        elif flip in ('x', 'h'): # about x-axis (horizontal)
            pimg = img[::-1, :]
        elif flip in ('vh', 'hv', 'r180'): # 180 degree rotation
            pimg = img[::-1, ::-1]
        elif flip in ('t', 'T'): # transpose (possible shape change)
            pimg = img.T
        elif flip in ('ccw90', 'r90'): # rotate 90 (possible shape change)
            pimg = img.T[:, ::-1]
        elif flip in ('cw90', 'r270'): # rotate 270 (possible shape change)
            pimg = img.T[::-1, :]
        else:
            pimg = img

        return pimg

    def _toarray(self, nframes=0):
        mynf = len(self)
        nf = np.min((mynf, nframes)) if nframes > 0 else mynf
        ashp = (nf,) + self.shape
        a = np.zeros(ashp, dtype=self.dtype)
        for i in range(nf):
            a[i] = self.__getitem__(i)

        return a

    @property
    def dtype(self):
        return self._imser.dtype

    @property
    def shape(self):
        return self._imser.shape

    def median(self, nframes=0):
        return np.median(self._toarray(nframes=nframes), axis=0)

    pass # end class

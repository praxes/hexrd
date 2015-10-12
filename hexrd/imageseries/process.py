"""Class for processing individual frames"""
import numpy as np

from .baseclass import ImageSeries

class ProcessedImageSeries(ImageSeries):
    """Images series with mapping applied to frames"""
    FLIP = 'flip'
    DARK = 'dark'
    RECT = 'rectangle'

    _opdict = {}

    def __init__(self, imser, oplist):
        """imsageseries based on existing one with image processing options

        *imser* - an existing imageseries
        *oplist* - list of processing operations;
                   a list of pairs (key, data) pairs, with key specifying the
                   operation to perform using specified data

        """
        self._imser = imser
        self._oplist = oplist

        self.addop(self.DARK, self._subtract_dark)
        self.addop(self.FLIP, self._flip)
        self.addop(self.RECT, self._rectangle)

    def __getitem__(self, key):
        return self._process_frame(key)

    def __len__(self):
        return len(self._imser)

    def _process_frame(self, key):
        img = np.copy(self._imser[key])
        for op in self.oplist:
            key, data = op
            func = self._opdict[key]
            img = func(img, data)

        return img

    def _subtract_dark(self, img, dark):
        # need to check for values below zero
        return np.where(img > dark, img-dark, 0)

    def _rectangle(self, img, r):
        # restrict to rectangle
        return img[r[0,0]:r[0,1], r[1,0]:r[1,1]]

    def _flip(self, img, flip):
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
    #
    # ==================== API
    #
    @property
    def dtype(self):
        return self._imser.dtype

    @property
    def shape(self):
        return self._imser.shape

    @classmethod
    def addop(cls, key, func):
        """Add operation to processing options

        *key* - string to use to specify this op
        *func* - function to call for this op: f(data)
        """
        cls._opdict[key] = func

    @property
    def oplist(self):
        """list of operations to apply"""
        return self._oplist

    pass # end class

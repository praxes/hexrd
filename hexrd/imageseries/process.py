"""Class for processing frames or frame groups"""
from hexrd.imageseries import ImageSeries

class ProcessedImageSeries(ImageSeries):
    """Images series with mapping applied to frames"""
    FLIP = 'flip'
    DARK = 'dark'
    
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
        img = self._imser[key]
        img = self._subtract_dark(img)
        img = self._flip(img)
        return img

    def _subtract_dark(self, img):
        if self.DARK in self._opts:
            return img - self._opts[self.DARK]

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
        

    @property 
    def dtype(self):
        return self._imser.dtype

    @property
    def shape(self):
        return self._imser.shape

    
    pass # end class

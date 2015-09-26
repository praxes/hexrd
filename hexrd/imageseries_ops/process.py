"""Class for processing frames or frame groups"""

class ModifiedImageSeries(Imageseries):
    """Images series with mapping applied to frames"""
    def __init__(self, imser, cfg):
        """Instantiate imsageseries based on existing one with mapping options

        *imser* - an existing imageseries
        *cfg* - configuration for processing options
        """
        self._imser = imser

    def __getitem__(self, key):
        return self._imser.__getitem__(key)
        
    def __len__(self):
        return len(self._imser)

    @property 
    def dtype(self):
        return self._imser.dtype

    @property
    def shape(self):
        return self._imser.shape

    def process_frame(self, key):
        img = self._imser[key]
        return self._process_frame(img)
    
    pass # end class

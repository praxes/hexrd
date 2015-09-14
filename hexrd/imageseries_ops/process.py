"""Class for processing frames or frame groups"""

class ModifiedImageSeries(Imageseries):
    """Images series with mapping applied to frames"""
    def __init__(self, imser, cfg):
        """Instantiate imsageseries based on existing one with mapping options

        *imser* - an existing imageseries
        *cfg* - configuration for processing options
        """
        self._imser = imser
    pass # end class

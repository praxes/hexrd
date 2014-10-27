import os
import logging
import multiprocessing as mp

from .config import Config
#from .detector import DetectorConfig
#from .findorientations import FindOrientationsConfig
#from .fitgrains import FitGrainsConfig
from .imageseries import ImageSeriesConfig
from .material import MaterialConfig


logger = logging.getLogger(__name__)


class RootConfig(Config):


    __config_map = {
        #'detector': DetectorConfig,
        #'find_orientations': FindOrientationsConfig,
        #'fit_grains': FitGrainsConfig,
        'image_series': ImageSeriesConfig,
        'material': MaterialConfig,
        }


    def _get_config(self, ctype):
        temp = self._cfg.get(ctype, None)
        temp = temp if temp is not None else {}
        return self.__config_map[ctype](temp, self._cfg)


    @property
    def analysis_name(self):
        return str(self._cfg.get('analysis_name', 'analysis'))


    @property
    def detector(self):
        return self._get_config('detector')


    @property
    def find_orientations(self):
        return self._get_config('find_orientations')


    @property
    def fit_grains(self):
        return self._get_config('fit_grains')


    @property
    def image_series(self):
        return self._get_config('image_series')


    @property
    def material(self):
        return self._get_config('material')


    @property
    def multiprocessing(self):
        # determine number of processes to run in parallel
        multiproc = self._cfg.get('multiprocessing', -1)
        ncpus = mp.cpu_count()
        if multiproc == 'all':
            return ncpus
        elif multiproc == 'half':
            temp = ncpus / 2
            return temp if temp else 1
        elif isinstance(multiproc, int):
            if multiproc >= 0:
                if multiproc > ncpus:
                    logger.warning(
                        'Resuested %s processes, %d available, using %d',
                        multiproc, ncpus, ncpus
                        )
                    return ncpus
                return multiproc if multiproc else 1
            else:
                temp = ncpus + multiproc
                if temp < 1:
                    logger.warning(
                        'Cannot use less than 1 process, requested %d of %d'
                        ' available. Defaulting to 1',
                        temp, ncpus
                        )
                    return 1
                return temp
        else:
            temp = ncpus - 1
            logger.warning(
                "Invalid value %s for multiprocessing, defaulting to %d"
                " of %d available",
                multiproc, temp, ncpus
                )
            return temp


    @property
    def working_dir(self):
        temp = self._cfg.get('working_dir', os.getcwd())
        if not os.path.isdir(temp):
            raise IOError(
                'working directory "%s" not found' % temp
                )
        return temp

import os
import logging
import multiprocessing as mp
import sys

from .config import Config
from .detector import DetectorConfig
from .findorientations import FindOrientationsConfig
#from .fitgrains import FitGrainsConfig
from .imageseries import ImageSeriesConfig
from .material import MaterialConfig
from utils import null


logger = logging.getLogger('hexrd.config')


class RootConfig(Config):


    __config_map = {
        'detector': DetectorConfig,
        'find_orientations': FindOrientationsConfig,
        #'fit_grains': FitGrainsConfig,
        'image_series': ImageSeriesConfig,
        'material': MaterialConfig,
        }


    def get(self, key, default=null, path_exists=False):
        args = key.split(':')
        args, item = args[:-1], args[-1]
        temp = self._cfg
        for arg in args:
            temp = temp.get(arg, {})
            # intermediate block may be None:
            temp = {} if temp is None else temp
        try:
            res = temp[item]
        except KeyError:
            if default is not null:
                logger.info(
                    '%s not specified, defaulting to %s', key, default
                    )
                res = temp.get(item, default)
            else:
                raise RuntimeError(
                    '%s must be specified in configuration file' % key
                    )
        if path_exists and res and not os.path.exists(res):
            raise IOError(
                '%s "%s" not found' % (key, res)
                )
        return res



    def _get_config(self, ctype):
        temp = self._cfg.get(ctype, None)
        temp = temp if temp is not None else {}
        return self.__config_map[ctype](self)


    @property
    def analysis_name(self):
        return str(self.get('analysis_name', default='analysis'))


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
        multiproc = self.get('multiprocessing', default=-1)
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
        return self.get(
            'working_dir', default=os.getcwd(), path_exists=True)

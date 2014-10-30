import os
import logging
import multiprocessing as mp
import sys

from .config import Config
from .detector import DetectorConfig
from .findorientations import FindOrientationsConfig
from .fitgrains import FitGrainsConfig
from .imageseries import ImageSeriesConfig
from .material import MaterialConfig
from utils import null


logger = logging.getLogger('hexrd.config')


class RootConfig(Config):


    def get(self, key, default=null):
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
        return res


    @property
    def analysis_name(self):
        return str(self.get('analysis_name', default='analysis'))


    @property
    def analysis_dir(self):
        return os.path.join(self.working_dir, self.analysis_name)


    @property
    def detector(self):
        return DetectorConfig(self)


    @property
    def find_orientations(self):
        return FindOrientationsConfig(self)


    @property
    def fit_grains(self):
        return FitGrainsConfig(self)


    @property
    def image_series(self):
        return ImageSeriesConfig(self)


    @property
    def material(self):
        return MaterialConfig(self)


    @property
    def multiprocessing(self):
        # determine number of processes to run in parallel
        multiproc = self.get('multiprocessing', default=-1)
        ncpus = mp.cpu_count()
        if multiproc == 'all':
            res = ncpus
        elif multiproc == 'half':
            temp = ncpus / 2
            res = temp if temp else 1
        elif isinstance(multiproc, int):
            if multiproc >= 0:
                if multiproc > ncpus:
                    logger.warning(
                        'Resuested %s processes, %d available',
                        multiproc, ncpus, ncpus
                        )
                    res = ncpus
                else:
                    res = multiproc if multiproc else 1
            else:
                temp = ncpus + multiproc
                if temp < 1:
                    logger.warning(
                        'Cannot use less than 1 process, requested %d of %d',
                        temp, ncpus
                        )
                    res = 1
                else:
                    res = temp
        else:
            temp = ncpus - 1
            logger.warning(
                "Invalid value %s for multiprocessing",
                multiproc
                )
            res = temp
        logger.info("Using %d of %d available processors", res, ncpus)
        return res


    @property
    def working_dir(self):
        temp = self.get('working_dir', default=os.getcwd())
        if os.path.exists(temp):
            return temp
        raise IOError(
            '"working_dir": "%s" does not exist'
            )

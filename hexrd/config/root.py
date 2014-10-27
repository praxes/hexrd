import os
import logging
import multiprocessing as mp

logger = logging.getLogger(__name__)


class ConfigRoot(object):


    @property
    def analysis_name(self):
        return self._cfg.get('analysis_name', '')


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
            raise IOError()
        return temp


    def __init__(self, cfg):
        self._cfg = cfg

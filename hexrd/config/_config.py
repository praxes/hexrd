import copy
import logging
import multiprocessing as mp
import os

import yaml

logger = logging.getLogger(__name__)


class Config(object):


    @property
    def analysis_name(self):
        return self._cfg.get('analysis_name', '')


    @property
    def multiprocessing(self):
        # determine number of processes to run in parallel
        multiproc = self._cfg.get('multiprocessing', -1)
        ncpus = mp.cpu_count()
        if multiproc == 'all':
            pass
        elif multiproc == 'half':
            ncpus /= 2
        elif isinstance(multiproc, int):
            ncpus = multiproc if multiproc >= 0 else ncpus+multiproc
        else:
            ncpus -= 1
            logger.warning(
                "Invalid value %s for multiprocessing, defaulting to -1"
                % multiproc
                )

        # if ncpus is 0, make it 1
        if ncpus < 1:
            ncpus = 1
            logger.warning(
                'Cannot use less than 1 process, requested %d. '
                'Defaulting to 1' % (ncpus)
                )
        elif ncpus > mp.cpu_count():
            ncpus = mp.cpu_count()
            logger.warning(
                'Resuested %s processes, using %d'
                % (multiproc, ncpus)
                )
        return ncpus


    @property
    def working_dir(self):
        return self._cfg.get('working_dir', os.getcwd())


    def __init__(self, cfg):
        self._cfg = cfg


def merge_dicts(a, b):
    "Returns a merged dict, updating values in `a` with values from `b`"
    a = copy.deepcopy(a)
    for k,v in b.iteritems():
        if isinstance(v, dict):
            merge_dicts(a[k], v)
        else:
            a[k] = v
    return a


def open(file_name):
    with file(file_name) as f:
        res = []
        for cfg in yaml.load_all(f):
            try:
                # take the previous config section and update with values
                # from the current one
                res.append(merge_dicts(res[-1], cfg))
            except IndexError:
                # this is the first config section
                res.append(cfg)
        return [Config(i) for i in res]

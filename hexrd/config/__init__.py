import yaml

from . import root
from . import utils


def open(file_name):
    """
    Reads configuration settings from a yaml file.

    Returns a list of configuration objects, one for each document section in
    the file.
    """
    with file(file_name) as f:
        res = []
        for cfg in yaml.load_all(f):
            try:
                # take the previous config section and update with values
                # from the current one
                res.append(utils.merge_dicts(res[0], cfg))
            except IndexError:
                # this is the first config section
                res.append(cfg)
        return [root.RootConfig(i, i) for i in res]

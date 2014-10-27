class Config(object):

    def __init__(self, cfg, rootcfg):
        self._cfg = cfg
        self._rootcfg = rootcfg


    def _get_nested_val(self, *args):
        args, item = args[:-1], args[-1]
        temp = self._cfg
        for arg in args:
            temp = temp.get(arg, {})
            temp = {} if temp is None else temp
        return temp[item]

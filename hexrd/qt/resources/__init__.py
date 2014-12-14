import os


class _ResourceManager(object):

    _root = os.path.split(__file__)[0]

    def __getitem__(self, key):
        return os.path.join(self._root, key)


resources = _ResourceManager()

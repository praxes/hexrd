import os


class _ResourceManager(object):

    _root = os.path.split(__file__)[0]

    def __init__(self):
        self._data = {}

    def __getitem__(self, key):
        return os.path.join(self._root, self._data[key])

    def __setitem__(self, key, val):
        self._data[key] = val


image_files = _ResourceManager()
image_files['splash'] = 'hexrd.png'

ui_files = _ResourceManager()
ui_files['main_window'] = 'mainwindow.ui'

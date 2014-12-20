from PyQt4 import QtGui, QtCore, uic

from .resources import resources


class ImageSeries(object):

    # this is just a placeholder
    pass



class ImageSeriesController(QtGui.QDialog):


    @property
    def image_series(self):
        return None


    def __init__(self, cfg, name=None):
        super(ImageSeriesController, self).__init__()
        uic.loadUi(resources['imageseries.ui'], self)
        self._configure_tool_buttons()

        self.rotateCheckBox.toggled.emit(False)

        self.load_settings()


    def _configure_tool_buttons(self):
        self.imagesToolButton.addAction(self.actionGetImageFiles)
        self.darkImageToolButton.addAction(self.actionGetDarkFile)
        self.pathToolButton.addAction(self.actionGetPath)


    def load_settings(self):
        pass


    @QtCore.pyqtSlot(name='on_buttonBox_accepted')
    def save_settings(self):
        pass


    @QtCore.pyqtSlot()
    def on_actionGetImageFiles_triggered(self):
        pass


    @QtCore.pyqtSlot(name='on_actionGetImageFiles_triggered')
    @QtCore.pyqtSlot(name='on_actionGetPath_triggered')
    @QtCore.pyqtSlot(name='on_actionGetDarkFile_triggered')
    def not_implemented(self):
        dlg = QtGui.QMessageBox.information(
            self, 'Not Implemented',
"""The requested feature has not been implemented.

Please consider filing an issue report at
http://github.com/praxes/hexrd/issues, if one does not already exist.""",
            buttons=QtGui.QMessageBox.Ok
            )



def get_image_series(cfg, name=None, ui=False):
    if ui:
        ctlr = ImageSeriesController(cfg, name)
        ctlr.exec_()
    return ctlr.image_series

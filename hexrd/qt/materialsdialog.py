from PyQt4 import QtGui, QtCore, uic

from .resources import resources


class Material(object):

    # this is just a placeholder
    pass



class MaterialDialogController(QtGui.QDialog):


    @property
    def material(self):
        return None


    def __init__(self, name=None):
        super(MaterialDialogController, self).__init__()
        uic.loadUi(resources['materialsdialog.ui'], self)

        self.load_settings()


    def load_settings(self):
        pass


    @QtCore.pyqtSlot(name='on_buttonBox_accepted')
    def save_settings(self):
        pass


    def not_implemented(self):
        dlg = QtGui.QMessageBox.information(
            self, 'Not Implemented',
"""The requested feature has not been implemented.

Please consider filing an issue report at
http://github.com/praxes/hexrd/issues, if one does not already exist.""",
            buttons=QtGui.QMessageBox.Ok
            )



def get_material(name=None, ui=False):
    if ui:
        ctlr = MaterialDialogController(name)
        ctlr.exec_()
    return ctlr.material

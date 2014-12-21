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
        self._config_ui()

        self.load_settings()


    def _config_ui(self):
        self.spaceGroupComboBox.lineEdit().setReadOnly(True)
        for i in range(1,231):
            self.spaceGroupComboBox.addItem(str(i))
        self.spaceGroupComboBox.lineEdit().setAlignment(QtCore.Qt.AlignCenter)
        self.hallSymbolComboBox.lineEdit().setReadOnly(True)
        self.hermannMauguinComboBox.lineEdit().setReadOnly(True)

        self.spaceGroupComboBox.setCurrentIndex(224)


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

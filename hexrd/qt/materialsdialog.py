from PyQt4 import QtGui, QtCore, uic

from hexrd.xrd import spacegroup
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
        self.hallSymbolComboBox.lineEdit().setReadOnly(True)
        self.hermannMauguinComboBox.lineEdit().setReadOnly(True)
        for k in spacegroup.sgid_to_hall:
            self.spaceGroupComboBox.addItem(k)
            self.hallSymbolComboBox.addItem(spacegroup.sgid_to_hall[k])
            self.hermannMauguinComboBox.addItem(spacegroup.sgid_to_hm[k])
        self.spaceGroupComboBox.lineEdit().setAlignment(QtCore.Qt.AlignCenter)
        self.hallSymbolComboBox.lineEdit().setAlignment(QtCore.Qt.AlignCenter)
        self.hermannMauguinComboBox.lineEdit().setAlignment(
            QtCore.Qt.AlignCenter
            )

        self.spaceGroupComboBox.setCurrentIndex(522)


    def load_settings(self):
        pass


    @QtCore.pyqtSlot(name='on_buttonBox_accepted')
    def save_settings(self):
        pass


    @QtCore.pyqtSlot(int, name='on_spaceGroupComboBox_currentIndexChanged')
    @QtCore.pyqtSlot(int, name='on_hallSymbolComboBox_currentIndexChanged')
    @QtCore.pyqtSlot(int, name='on_hermannMauguinComboBox_currentIndexChanged')
    def set_spacegroup(self, val):
        try:
            self.spaceGroupComboBox.blockSignals(True)
            self.hallSymbolComboBox.blockSignals(True)
            self.hermannMauguinComboBox.blockSignals(True)
            self.spaceGroupComboBox.setCurrentIndex(val)
            self.hallSymbolComboBox.setCurrentIndex(val)
            self.hermannMauguinComboBox.setCurrentIndex(val)
            sgid = int(self.spaceGroupComboBox.currentText().split(':')[0])
            for sgids, lg in spacegroup._pgDict.items():
                if sgid in sgids:
                    self.laueGroupLineEdit.setText(lg[0])
                    break
            self.latticeTypeLineEdit.setText(spacegroup._ltDict[lg[1]])
        finally:
            self.spaceGroupComboBox.blockSignals(False)
            self.hallSymbolComboBox.blockSignals(False)
            self.hermannMauguinComboBox.blockSignals(False)


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

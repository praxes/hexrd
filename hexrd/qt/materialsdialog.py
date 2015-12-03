from PyQt4 import QtGui, QtCore, uic

from hexrd.xrd import spacegroup
from hexrd.xrd.material import Material
from .resources import resources

class MaterialDialogController(QtGui.QDialog):


    def __init__(self, name=None):

        super(MaterialDialogController, self).__init__()
        uic.loadUi(resources['materialsdialog.ui'], self)

        self._mat = Material()
        # in order: a, b, c, alpha, beta, gamma
        self._lpboxes = [self.doubleSpinBox_2, self.doubleSpinBox_4,
                         self.doubleSpinBox_3, self.doubleSpinBox_7,
                         self.doubleSpinBox_5, self.doubleSpinBox_6 ]
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

    def _lp_blocksignals(self, tf):
        self.doubleSpinBox_2.blockSignals(tf)
        self.doubleSpinBox_3.blockSignals(tf)
        self.doubleSpinBox_4.blockSignals(tf)
        self.doubleSpinBox_5.blockSignals(tf)
        self.doubleSpinBox_6.blockSignals(tf)
        self.doubleSpinBox_7.blockSignals(tf)

    def load_settings(self):
        pass

    @property
    def material(self):
        return self._mat

    @property
    def lpboxes(self):
        return self._lpboxes

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

    @QtCore.pyqtSlot(int, name='on_spaceGroupComboBox_currentIndexChanged')
    def enable_latparams(self, val):
        """enable independent lattice parameters"""
        # lattice parameters are stored in the old "ValUnit" class
        self._lp_blocksignals(True)
        nlp = 6
        A  = 'angstrom'
        D  = 'degrees'
        m = self.material
        sgid = int(self.spaceGroupComboBox.currentText().split(':')[0])
        self.material.sgnum = sgid
        reqp = m.spaceGroup.reqParams
        lprm = m.latticeParameters
        for i in range(nlp):
            boxi = self.lpboxes[i]
            boxi.setEnabled(i in reqp)
            u = A if i < 3 else D
            boxi.setValue(lprm[i].getVal(u))
        self._lp_blocksignals(False)

    @QtCore.pyqtSlot(float, name='on_doubleSpinBox_2_valueChanged')
    @QtCore.pyqtSlot(float, name='on_doubleSpinBox_3_valueChanged')
    @QtCore.pyqtSlot(float, name='on_doubleSpinBox_4_valueChanged')
    @QtCore.pyqtSlot(float, name='on_doubleSpinBox_5_valueChanged')
    @QtCore.pyqtSlot(float, name='on_doubleSpinBox_6_valueChanged')
    @QtCore.pyqtSlot(float, name='on_doubleSpinBox_7_valueChanged')
    def set_latparams(self, val):
        """update all the lattice parameter boxes when one changes"""
        # note: material takes reduced set of lattice parameters but outputs
        #       all six
        self._lp_blocksignals(True)
        nlp = 6
        A  = 'angstrom'
        D  = 'degrees'
        m = self.material
        reqp = m.spaceGroup.reqParams
        nreq = len(reqp)
        lp_red = nreq*[0.0]
        for i in range(nreq):
            boxi = self.lpboxes[reqp[i]]
            lp_red[i] = boxi.value()
        m.latticeParameters = lp_red
        lprm = m.latticeParameters
        for i in range(nlp):
            u = A if i < 3 else D
            boxi = self.lpboxes[i]
            boxi.setValue(lprm[i].getVal(u))
        self._lp_blocksignals(False)

    @QtCore.pyqtSlot(name='on_lineEdit_6_editingFinished')
    def set_name(self):
        n = self.lineEdit_6.text()
        print n, type(n)
        self.material.name = self.lineEdit_6.text()

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

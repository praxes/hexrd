import logging
import multiprocessing
import os
import sys
import time

# use API v2 to prepare for migration to Py3/PyQt5:
import sip
sip.setapi('QString', 2)
sip.setapi('QTextStream', 2)
sip.setapi('QVariant', 2)

from PyQt4.QtCore import Qt, QObject, QSettings, pyqtSlot
from PyQt4.QtGui import (
    qApp, QApplication, QFileDialog, QMainWindow, QMessageBox, QPixmap,
    QSplashScreen
    )
from PyQt4.uic import loadUi

from matplotlib import cm

from hexrd import config
from hexrd.qt.resources import image_files, ui_files


class QLogStream(object):

    def __init__(self, dest):
        self._dest = dest

    def write(self, val):
        # append adds its own line endings, need to rstrip
        self._dest.append(val.rstrip())



def add_handler(log_level, stream=None):
    logger = logging.getLogger('hexrd')
    h = logging.StreamHandler(stream)
    h.setLevel(log_level)
    h.setFormatter(
        logging.Formatter('%(asctime)s - %(message)s', '%y-%m-%d %H:%M:%S')
        )
    logger.addHandler(h)


class GraphicsCanvasController(QObject):

    def __init__(self, ui):
        self.ui = ui

        cmaps = sorted(i[:-2] for i in dir(cm) if i.endswith('_r'))
        ui.colorMapComboBox.addItems(cmaps)
        ui.colorMapComboBox.setCurrentIndex(cmaps.index('gist_heat'))

        ui.minSpinBox.valueChanged[int].connect(
            ui.maxSpinBox.setMinimum
            )
        ui.maxSpinBox.valueChanged[int].connect(
            ui.minSpinBox.setMaximum
            )


class MainController(QMainWindow):


    def __init__(self, log_level, cfg=None):
        super(MainController, self).__init__()
        # Create and display the splash screen
        splash_pix = QPixmap(image_files['splash'])
        splash = QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)
        splash.setMask(splash_pix.mask())
        splash.show()
        # sleep seems necessary on linux to render the image
        time.sleep(0.01)
        qApp.processEvents()

        # give the splash screen a little time to breathe
        time.sleep(2)

        loadUi(ui_files['main_window'], self)

        # now that we have the ui, configure the logging widget
        add_handler(log_level, QLogStream(self.loggerTextEdit))

        self.gc_ctlr = GraphicsCanvasController(self)

        self._connect_signals()

        self.load_config(cfg)

        self.setWindowTitle('HEXRD')
        self.settings = QSettings('hexrd', 'hexrd')
        try:
            self.restoreGeometry(self.settings.value('geometry'))
            self.restoreState(self.settings.value('state'))
            self.cfgToolBox.setCurrentIndex(
                self.settings.value('currentTool', type=int)
                )
        except TypeError:
            raise
        self.show()
        splash.finish(self)


    def _connect_signals(self):
        self.actionExit.triggered.connect(self.close)


    @pyqtSlot()
    def on_actionDocumentation_triggered(self):
        import webbrowser
        webbrowser.open_new_tab('http://hexrd.readthedocs.org/en/latest')


    @pyqtSlot()
    def on_actionAbout_triggered(self):
        dlg = QMessageBox.about(
            self, 'About HEXRD',
"""HEXRD provides a collection of resources for analysis of x-ray diffraction
data, especially high-energy x-ray diffraction. HEXRD is comprised of a
library and API for writing scripts, a command line interface, and an
interactive graphical user interface.

HEXRD is an open-source project originally conceived by Joel Bernier, and
developed by Joel Bernier, Darren Dale, and Donald Boyce, et.al.
"""
            )


    def on_analysisNameLineEdit_editingFinished(self):
        self.cfg.analysis_name = str(self.analysisNameLineEdit.text())


    @pyqtSlot()
    def on_changeWorkingDirButton_clicked(self):
        temp = QFileDialog.getExistingDirectory(
            parent=self, caption='booya', directory=self.cfg.working_dir
            )
        if temp:
            self.workingDirLineEdit.setText(temp)
            self.cfg.working_dir = temp


    @pyqtSlot(int)
    def on_multiprocessingSpinBox_valueChanged(self, val):
        if self.cfg.multiprocessing != val:
            self.cfg.multiprocessing = val


    def closeEvent(self, event):
        if self.cfg.dirty:
            confirm = QMessageBox()
            confirm.setText('Configuration has been modified.')
            confirm.setInformativeText('Do you want to save your changes?')
            confirm.setStandardButtons(
                QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel
                )
            confirm.setDefaultButton(QMessageBox.Cancel)
            ret = confirm.exec_()
            if ret == QMessageBox.Save:
                if self.save_config() is False:
                    event.ignore()
                    return
            elif ret == QMessageBox.Cancel:
                return
        self.settings.setValue('geometry', self.saveGeometry())
        self.settings.setValue('state', self.saveState())
        self.settings.setValue('currentTool', self.cfgToolBox.currentIndex())
        event.accept()


    def load_config(self, cfg):
        self.cfg = cfg

        # general
        self.analysisNameLineEdit.setText(self.cfg.analysis_name)

        self.workingDirLineEdit.setText(self.cfg.working_dir)

        self.multiprocessingSpinBox.setMaximum(multiprocessing.cpu_count())
        self.multiprocessingSpinBox.setValue(self.cfg.multiprocessing)


    @pyqtSlot(name='on_actionSaveCalibration_triggered')
    @pyqtSlot(name='on_actionSaveMaterials_triggered')
    @pyqtSlot(name='on_actionPowderBinnedFit_triggered')
    @pyqtSlot(name='on_actionPowderDirectFit_triggered')
    @pyqtSlot(name='on_actionSingleCrystalFit_triggered')
    @pyqtSlot(name='on_actionCake_triggered')
    @pyqtSlot(name='on_actionPolarRebin_triggered')
    @pyqtSlot(name='on_actionFindOrientations_triggered')
    @pyqtSlot(name='on_actionFitGrains_triggered')
    def not_implemented(self):
        dlg = QMessageBox.information(
            self, 'Not Implemented',
"""The requested feature has not been implemented.

Please consider filing an issue report at
http://github.com/praxes/hexrd/issues, if one does not already exist.""",
            buttons=QMessageBox.Ok
            )


    @pyqtSlot(name='on_actionSaveConfiguration_triggered')
    def save_config(self):
        temp = QFileDialog.getSaveFileName(
            parent=self, caption='Save Configuration',
            directory=self.cfg.working_dir, filter='YAML files (*.yml)'
            )
        if temp:
            self.cfg.dump(temp)
            return True
        return False



def execute(args):
    app = QApplication(sys.argv)
    app.setApplicationName('HEXRD')

    # configure logging
    if args.debug:
        log_level = logging.DEBUG
        add_handler(log_level)
    else:
        log_level = logging.CRITICAL if args.quiet else logging.INFO

    cfg = config.open(args.config)[0]

    ctlr = MainController(log_level, cfg)

    sys.exit(app.exec_())
